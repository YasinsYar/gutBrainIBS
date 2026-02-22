from __future__ import annotations

import argparse
import json
import queue
import re
import shlex
import signal
import subprocess
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

import psutil


RULE_OUTPUTS = {
    "validate_inputs": ["data_processed/qc/input_validation.json"],
    "ingest_eqtl_colon_local": ["data_raw/qtl/colon_manifest.tsv"],
    "ingest_eqtl_brain": ["data_raw/qtl/brain_manifest.tsv"],
    "parse_eqtl_to_parquet": ["data_interim/eqtl_parsed.parquet", "data_interim/eqtl_parse_manifest.tsv"],
    "harmonize_variants": ["data_interim/eqtl_harmonized.parquet", "data_processed/qc/harmonization_qc.tsv"],
    "filter_significant_eqtl": ["data_interim/eqtl_significant.parquet"],
    "gwas_harmonization": ["data_interim/gwas_harmonized.parquet"],
    "shared_specific_colon": ["data_processed/tables/shared_specific_colon.tsv"],
    "heterogeneity_colon": ["data_processed/tables/heterogeneity_colon.tsv"],
    "coloc_colon_and_brain": ["data_processed/tables/coloc_colon.tsv", "data_processed/tables/coloc_brain.tsv"],
    "genetic_correlation": ["data_processed/tables/genetic_correlation.tsv"],
    "bidirectional_mr": ["data_processed/tables/mr_results.tsv"],
    "pathway_enrichment": ["data_processed/tables/pathway_enrichment.tsv"],
    "scrna_support_colon": ["data_processed/tables/scrna_support.tsv"],
    "build_figures": ["figures/pathway_enrichment_dotplot.png"],
    "assemble_results": ["data_processed/final/hypothesis_scorecard.json", "reports/run_qc.html"],
    "coloc_sensitivity": [
        "data_processed/sensitivity/coloc_sensitivity_summary.tsv",
        "data_processed/sensitivity/coloc_sensitivity_brain_detail.tsv",
    ],
    "mr_sensitivity": ["data_processed/sensitivity/mr_robustness.tsv"],
    "harmonization_sensitivity": [
        "data_processed/sensitivity/harmonization_strict_vs_nonstrict.tsv",
        "data_processed/sensitivity/coloc_brain_nonstrict.tsv",
    ],
    "run_sensitivity": ["data_processed/sensitivity/sensitivity_summary.tsv"],
    "build_presentation_assets": [
        "figures/presentation/10_mr_forest_robustness.png",
        "data_processed/tables/presentation/presentation_manifest.json",
        "reports/presentation_tables.html",
    ],
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def fmt_seconds(value: float | None) -> str:
    if value is None or value < 0 or not (value < 10**9):
        return "n/a"
    seconds = int(value)
    h = seconds // 3600
    m = (seconds % 3600) // 60
    s = seconds % 60
    return f"{h:02d}:{m:02d}:{s:02d}"


def progress_bar(done: int, total: int, width: int = 32) -> str:
    if total <= 0:
        return "[" + ("-" * width) + "]"
    frac = min(max(done / total, 0.0), 1.0)
    filled = int(round(frac * width))
    return "[" + ("#" * filled) + ("-" * (width - filled)) + "]"


def tree_processes(proc: psutil.Process) -> list[psutil.Process]:
    try:
        children = proc.children(recursive=True)
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        children = []
    return [proc] + children


def tree_metrics(proc: psutil.Process) -> tuple[float, float, float]:
    rss = 0
    cpu = 0.0
    io_write = 0
    for p in tree_processes(proc):
        try:
            rss += p.memory_info().rss
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
        try:
            t = p.cpu_times()
            cpu += float(t.user) + float(t.system)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
        try:
            io = p.io_counters()
            io_write += int(io.write_bytes)
        except (psutil.NoSuchProcess, psutil.AccessDenied, AttributeError):
            pass
    return rss / (1024 * 1024), cpu, io_write / (1024 * 1024)


def collect_output_sizes(paths: list[str]) -> tuple[dict[str, int], int]:
    sizes: dict[str, int] = {}
    total = 0
    for rel in paths:
        p = Path(rel)
        if p.exists():
            s = p.stat().st_size
            sizes[rel] = s
            total += s
        else:
            sizes[rel] = 0
    return sizes, total


def dir_size_mb(path: Path) -> float:
    if not path.exists() or not path.is_dir():
        return 0.0
    total = 0
    for p in path.rglob("*"):
        if p.is_file():
            try:
                total += p.stat().st_size
            except OSError:
                continue
    return total / (1024 * 1024)


def terminate_tree(pid: int, timeout: float = 20.0) -> None:
    try:
        root = psutil.Process(pid)
    except psutil.NoSuchProcess:
        return
    procs = tree_processes(root)
    for p in reversed(procs):
        try:
            p.terminate()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    _, alive = psutil.wait_procs(procs, timeout=timeout)
    for p in alive:
        try:
            p.kill()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass


def write_status(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="all")
    parser.add_argument("--profile", default="profiles/windows-default")
    parser.add_argument("--configfile", default="config/config.yaml")
    parser.add_argument("--forceall", action="store_true")
    parser.add_argument("--heartbeat-sec", type=int, default=15)
    parser.add_argument("--log-file", default="reports/run_managed.log")
    parser.add_argument("--status-file", default="reports/run_status.json")
    parser.add_argument("--pid-file", default="reports/run_status.pid")
    args, passthrough = parser.parse_known_args()

    log_path = Path(args.log_file)
    status_path = Path(args.status_file)
    pid_path = Path(args.pid_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    if pid_path.exists():
        raw = pid_path.read_text(encoding="utf-8").strip()
        if raw.isdigit():
            old_pid = int(raw)
            if psutil.pid_exists(old_pid):
                payload = {
                    "state": "already_running",
                    "message": f"existing managed run detected pid={old_pid}",
                    "updated_at_utc": utc_now(),
                }
                write_status(status_path, payload)
                print(f"[managed] existing run pid={old_pid}. use run_control stop/cleanup first.")
                return 2
        pid_path.unlink(missing_ok=True)

    cmd = ["snakemake", "--profile", args.profile, f"--configfile={args.configfile}"]
    if args.forceall:
        cmd.append("--forceall")
    cmd.append(args.target)
    cmd.extend(passthrough)

    print(f"[managed] command: {' '.join(shlex.quote(x) for x in cmd)}")
    print(f"[managed] status: {status_path}")
    print(f"[managed] log: {log_path}")

    start = time.time()
    last_hb = start
    last_output_ts = start
    last_done_ts = start
    total_steps: int | None = None
    done_steps = 0
    active_rules: set[str] = set()
    stopped_by_signal = False

    re_total = re.compile(r"^\s*total\s+(\d+)\s*$")
    re_done = re.compile(r"(\d+) of (\d+) steps \((\d+)%\) done")
    re_rule_start = re.compile(r"\b(?:localrule|rule)\s+([A-Za-z0-9_]+):")
    re_rule_finish = re.compile(r"Finished jobid:\s*\d+\s*\(Rule:\s*([A-Za-z0-9_]+)\)")
    re_harmonize = re.compile(r"\[harmonize_stream\] row_groups=(\d+)/(\d+) n_output=(\d+) elapsed=([0-9.]+)s")
    re_filter = re.compile(r"\[filter_eqtl_stream\] row_groups=(\d+)/(\d+) n_output=(\d+) elapsed=([0-9.]+)s")
    re_parse = re.compile(r"\[parse_eqtl_duckdb\] combine part=(\d+)/(\d+) rows_copied=(\d+)")
    re_harm_sens = re.compile(r"\[harmonization_sensitivity\] row_groups=(\d+)/(\d+) matched_chunks=(\d+) elapsed=([0-9.]+)s")
    re_coloc_sens = re.compile(r"\[coloc_sensitivity\] scenario=(\d+)/(\d+) id=([^\s]+) rows=(\d+)")

    with open(log_path, "w", encoding="utf-8") as log:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
        )
        root = psutil.Process(proc.pid)
        pid_path.write_text(str(proc.pid), encoding="utf-8")

        def _signal_handler(_signum, _frame) -> None:
            nonlocal stopped_by_signal
            stopped_by_signal = True
            terminate_tree(proc.pid)

        signal.signal(signal.SIGINT, _signal_handler)
        signal.signal(signal.SIGTERM, _signal_handler)

        out_queue: queue.Queue[str] = queue.Queue()

        assert proc.stdout is not None

        def _reader() -> None:
            for raw in proc.stdout:
                out_queue.put(raw)
            out_queue.put("")

        reader = threading.Thread(target=_reader, daemon=True)
        reader.start()

        startup = f"[progress] {progress_bar(0, 1)}   0% 0/? elapsed=00:00:00 eta=n/a active=n/a health=startup"
        print(startup)
        log.write(startup + "\n")

        tracked_paths = sorted({p for vals in RULE_OUTPUTS.values() for p in vals})
        _, prev_out_total = collect_output_sizes(tracked_paths)
        parts_dir = Path("data_interim/eqtl_parts")
        duckdb_tmp_dir = Path("data_interim/duckdb_tmp")
        prev_parts_mb = dir_size_mb(parts_dir)
        prev_tmp_mb = dir_size_mb(duckdb_tmp_dir)
        prev_cpu_total = 0.0
        prev_io_write = 0.0
        stream_done = False
        rc: int | None = None
        subprogress: dict | None = None

        while True:
            now = time.time()
            got_line = False
            try:
                line = out_queue.get(timeout=0.2)
                got_line = True
            except queue.Empty:
                line = None

            if got_line and line == "":
                stream_done = True
            elif got_line and line is not None:
                last_output_ts = now
                sys.stdout.write(line)
                log.write(line)
                log.flush()
                stripped = line.strip()

                mt = re_total.match(stripped)
                if mt:
                    total_steps = int(mt.group(1))
                md = re_done.search(stripped)
                if md:
                    done_steps = int(md.group(1))
                    total_steps = int(md.group(2))
                    last_done_ts = now
                ms = re_rule_start.search(stripped)
                if ms:
                    active_rules.add(ms.group(1))
                mf = re_rule_finish.search(stripped)
                if mf:
                    active_rules.discard(mf.group(1))

                mh = re_harmonize.search(stripped)
                if mh:
                    done = int(mh.group(1))
                    total = int(mh.group(2))
                    n_output = int(mh.group(3))
                    elapsed_sub = float(mh.group(4))
                    eta_sub = ((elapsed_sub / done) * (total - done)) if done > 0 and total > done else None
                    subprogress = {
                        "kind": "row_groups",
                        "step": "harmonize_variants",
                        "done": done,
                        "total": total,
                        "pct": round((done / total) * 100, 2) if total > 0 else 0.0,
                        "n_output": n_output,
                        "elapsed_sec": round(elapsed_sub, 1),
                        "eta_sec": round(eta_sub, 1) if eta_sub is not None else None,
                    }
                mfilt = re_filter.search(stripped)
                if mfilt:
                    done = int(mfilt.group(1))
                    total = int(mfilt.group(2))
                    n_output = int(mfilt.group(3))
                    elapsed_sub = float(mfilt.group(4))
                    eta_sub = ((elapsed_sub / done) * (total - done)) if done > 0 and total > done else None
                    subprogress = {
                        "kind": "row_groups",
                        "step": "filter_significant_eqtl",
                        "done": done,
                        "total": total,
                        "pct": round((done / total) * 100, 2) if total > 0 else 0.0,
                        "n_output": n_output,
                        "elapsed_sec": round(elapsed_sub, 1),
                        "eta_sec": round(eta_sub, 1) if eta_sub is not None else None,
                    }
                mp = re_parse.search(stripped)
                if mp:
                    done = int(mp.group(1))
                    total = int(mp.group(2))
                    rows_copied = int(mp.group(3))
                    subprogress = {
                        "kind": "parts",
                        "step": "parse_eqtl_to_parquet",
                        "done": done,
                        "total": total,
                        "pct": round((done / total) * 100, 2) if total > 0 else 0.0,
                        "rows_copied": rows_copied,
                    }
                mhs = re_harm_sens.search(stripped)
                if mhs:
                    done = int(mhs.group(1))
                    total = int(mhs.group(2))
                    matched = int(mhs.group(3))
                    elapsed_sub = float(mhs.group(4))
                    eta_sub = ((elapsed_sub / done) * (total - done)) if done > 0 and total > done else None
                    subprogress = {
                        "kind": "row_groups",
                        "step": "harmonization_sensitivity",
                        "done": done,
                        "total": total,
                        "pct": round((done / total) * 100, 2) if total > 0 else 0.0,
                        "matched_chunks": matched,
                        "elapsed_sec": round(elapsed_sub, 1),
                        "eta_sec": round(eta_sub, 1) if eta_sub is not None else None,
                    }
                mcs = re_coloc_sens.search(stripped)
                if mcs:
                    done = int(mcs.group(1))
                    total = int(mcs.group(2))
                    scenario_id = str(mcs.group(3))
                    rows = int(mcs.group(4))
                    subprogress = {
                        "kind": "scenarios",
                        "step": "coloc_sensitivity",
                        "done": done,
                        "total": total,
                        "pct": round((done / total) * 100, 2) if total > 0 else 0.0,
                        "scenario_id": scenario_id,
                        "rows": rows,
                    }

            if now - last_hb >= args.heartbeat_sec:
                elapsed = now - start
                if total_steps and done_steps > 0:
                    rate = done_steps / elapsed
                    remaining = max(total_steps - done_steps, 0)
                    eta_sec = remaining / rate if rate > 0 else None
                    pct = int(round((done_steps / total_steps) * 100))
                    bar = progress_bar(done_steps, total_steps)
                else:
                    eta_sec = None
                    pct = 0
                    bar = progress_bar(0, total_steps or 1)

                rss_mb, cpu_total, io_write_mb = tree_metrics(root)
                _, out_total = collect_output_sizes(tracked_paths)
                parts_mb = dir_size_mb(parts_dir)
                tmp_mb = dir_size_mb(duckdb_tmp_dir)
                delta_cpu = max(0.0, cpu_total - prev_cpu_total)
                delta_io = max(0.0, io_write_mb - prev_io_write)
                delta_out = max(0, out_total - prev_out_total) / (1024 * 1024)
                delta_parts = max(0.0, parts_mb - prev_parts_mb)
                delta_tmp = max(0.0, tmp_mb - prev_tmp_mb)
                prev_cpu_total = cpu_total
                prev_io_write = io_write_mb
                prev_out_total = out_total
                prev_parts_mb = parts_mb
                prev_tmp_mb = tmp_mb

                silent_for = now - last_output_ts
                no_done_for = now - last_done_ts
                if silent_for > args.heartbeat_sec * 6 and delta_cpu < 0.2 and delta_io < 1.0 and delta_out < 1.0 and delta_parts < 1.0 and delta_tmp < 1.0:
                    health = "likely_stalled"
                elif silent_for > args.heartbeat_sec * 2 and (delta_cpu >= 0.2 or delta_io >= 1.0 or delta_out >= 1.0 or delta_parts >= 1.0 or delta_tmp >= 1.0):
                    health = "working_no_logs"
                elif no_done_for > args.heartbeat_sec * 8:
                    health = "long_step"
                else:
                    health = "stable"

                active = ",".join(sorted(active_rules)) if active_rules else "n/a"
                msg = (
                    f"[progress] {bar} {pct:3d}% {done_steps}/{total_steps or '?'} "
                    f"elapsed={fmt_seconds(elapsed)} eta={fmt_seconds(eta_sec)} "
                    f"rss_mb={rss_mb:.1f} cpu_d={delta_cpu:.2f}s io_d={delta_io:.1f}MB out_d={delta_out:.1f}MB "
                    f"parts_d={delta_parts:.1f}MB tmp_d={delta_tmp:.1f}MB "
                    f"active={active} health={health}"
                )
                if subprogress:
                    msg += f" sub={subprogress.get('step')}:{subprogress.get('done')}/{subprogress.get('total')}({subprogress.get('pct')}%)"
                print(msg)
                log.write(msg + "\n")
                log.flush()

                status_payload = {
                    "state": "running",
                    "target": args.target,
                    "pid": proc.pid,
                    "started_at_utc": datetime.fromtimestamp(start, timezone.utc).isoformat(),
                    "updated_at_utc": utc_now(),
                    "done_steps": done_steps,
                    "total_steps": total_steps,
                    "elapsed_sec": elapsed,
                    "eta_sec": eta_sec,
                    "active_rules": sorted(active_rules),
                    "rss_mb": rss_mb,
                    "cpu_delta_sec": delta_cpu,
                    "io_write_delta_mb": delta_io,
                    "outputs_delta_mb": delta_out,
                    "parts_delta_mb": delta_parts,
                    "duckdb_tmp_delta_mb": delta_tmp,
                    "silent_for_sec": silent_for,
                    "health": health,
                    "log_file": str(log_path),
                    "command": cmd,
                }
                if subprogress:
                    status_payload["subprogress"] = subprogress
                write_status(status_path, status_payload)
                last_hb = now

            if stream_done and proc.poll() is not None and out_queue.empty():
                rc = proc.wait()
                break

        elapsed = time.time() - start
        if rc is None:
            rc = proc.wait()
        if stopped_by_signal or rc in (15, -15):
            final_state = "stopped"
        elif rc == 0:
            final_state = "completed"
        else:
            final_state = "failed"
        final_payload = {
            "state": final_state,
            "target": args.target,
            "pid": proc.pid,
            "return_code": rc,
            "elapsed_sec": elapsed,
            "done_steps": done_steps,
            "total_steps": total_steps,
            "updated_at_utc": utc_now(),
            "log_file": str(log_path),
        }
        write_status(status_path, final_payload)
        print(f"[managed] finished rc={rc} elapsed={fmt_seconds(elapsed)}")

    pid_path.unlink(missing_ok=True)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
