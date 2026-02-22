from __future__ import annotations

import argparse
import json
import re
import subprocess
import time
from pathlib import Path

import psutil


def process_tree(proc: psutil.Process) -> list[psutil.Process]:
    try:
        children = proc.children(recursive=True)
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        children = []
    return [proc] + children


def terminate_tree(pid: int, timeout: float = 20.0) -> bool:
    try:
        root = psutil.Process(pid)
    except psutil.NoSuchProcess:
        return False
    procs = process_tree(root)
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
    return True


def read_tail_lines(path: Path, max_bytes: int = 2_000_000) -> list[str]:
    if not path.exists():
        return []
    size = path.stat().st_size
    start = max(0, size - max_bytes)
    with open(path, "rb") as fh:
        fh.seek(start)
        blob = fh.read()
    text = blob.decode("utf-8", errors="ignore")
    return text.splitlines()


def enrich_subprogress(data: dict) -> dict:
    log_file = data.get("log_file")
    if not log_file:
        return data
    lines = read_tail_lines(Path(log_file))
    if not lines:
        return data

    re_h = re.compile(r"\[harmonize_stream\] row_groups=(\d+)/(\d+) n_output=(\d+) elapsed=([0-9.]+)s")
    re_f = re.compile(r"\[filter_eqtl_stream\] row_groups=(\d+)/(\d+) n_output=(\d+) elapsed=([0-9.]+)s")
    re_p = re.compile(r"\[parse_eqtl_duckdb\] combine part=(\d+)/(\d+) rows_copied=(\d+)")
    re_hs = re.compile(r"\[harmonization_sensitivity\] row_groups=(\d+)/(\d+) matched_chunks=(\d+) elapsed=([0-9.]+)s")
    re_cs = re.compile(r"\[coloc_sensitivity\] scenario=(\d+)/(\d+) id=([^\s]+) rows=(\d+)")

    for line in reversed(lines):
        mh = re_h.search(line)
        if mh:
            done = int(mh.group(1))
            total = int(mh.group(2))
            n_output = int(mh.group(3))
            elapsed = float(mh.group(4))
            eta = None
            if done > 0 and total > done and elapsed > 0:
                eta = round((elapsed / done) * (total - done), 1)
            data["subprogress"] = {
                "kind": "row_groups",
                "step": "harmonize_variants",
                "done": done,
                "total": total,
                "pct": round((done / total) * 100, 2) if total > 0 else 0.0,
                "n_output": n_output,
                "elapsed_sec": round(elapsed, 1),
                "eta_sec": eta,
            }
            return data

        mf = re_f.search(line)
        if mf:
            done = int(mf.group(1))
            total = int(mf.group(2))
            n_output = int(mf.group(3))
            elapsed = float(mf.group(4))
            eta = None
            if done > 0 and total > done and elapsed > 0:
                eta = round((elapsed / done) * (total - done), 1)
            data["subprogress"] = {
                "kind": "row_groups",
                "step": "filter_significant_eqtl",
                "done": done,
                "total": total,
                "pct": round((done / total) * 100, 2) if total > 0 else 0.0,
                "n_output": n_output,
                "elapsed_sec": round(elapsed, 1),
                "eta_sec": eta,
            }
            return data

        mp = re_p.search(line)
        if mp:
            done = int(mp.group(1))
            total = int(mp.group(2))
            copied = int(mp.group(3))
            data["subprogress"] = {
                "kind": "parts",
                "step": "parse_eqtl_to_parquet",
                "done": done,
                "total": total,
                "pct": round((done / total) * 100, 2) if total > 0 else 0.0,
                "rows_copied": copied,
            }
            return data
        mhs = re_hs.search(line)
        if mhs:
            done = int(mhs.group(1))
            total = int(mhs.group(2))
            matched = int(mhs.group(3))
            elapsed = float(mhs.group(4))
            eta = None
            if done > 0 and total > done and elapsed > 0:
                eta = round((elapsed / done) * (total - done), 1)
            data["subprogress"] = {
                "kind": "row_groups",
                "step": "harmonization_sensitivity",
                "done": done,
                "total": total,
                "pct": round((done / total) * 100, 2) if total > 0 else 0.0,
                "matched_chunks": matched,
                "elapsed_sec": round(elapsed, 1),
                "eta_sec": eta,
            }
            return data
        mcs = re_cs.search(line)
        if mcs:
            done = int(mcs.group(1))
            total = int(mcs.group(2))
            scenario_id = str(mcs.group(3))
            rows = int(mcs.group(4))
            data["subprogress"] = {
                "kind": "scenarios",
                "step": "coloc_sensitivity",
                "done": done,
                "total": total,
                "pct": round((done / total) * 100, 2) if total > 0 else 0.0,
                "scenario_id": scenario_id,
                "rows": rows,
            }
            return data
    return data


def cmd_status(status_file: Path) -> int:
    if not status_file.exists():
        print(f"status file not found: {status_file}")
        return 1
    data = json.loads(status_file.read_text(encoding="utf-8"))
    data = enrich_subprogress(data)
    print(json.dumps(data, ensure_ascii=False, indent=2))
    return 0


def fmt_hms(value: float | None) -> str:
    if value is None:
        return "n/a"
    secs = int(max(0, value))
    h = secs // 3600
    m = (secs % 3600) // 60
    s = secs % 60
    return f"{h:02d}:{m:02d}:{s:02d}"


def cmd_watch(status_file: Path, interval_sec: float, max_iterations: int) -> int:
    iterations = 0
    while True:
        if not status_file.exists():
            print(f"status file not found: {status_file}")
            return 1
        data = json.loads(status_file.read_text(encoding="utf-8"))
        data = enrich_subprogress(data)
        state = data.get("state", "unknown")
        active_rules = ",".join(data.get("active_rules", [])) or "n/a"
        done_steps = data.get("done_steps")
        total_steps = data.get("total_steps")
        elapsed = fmt_hms(data.get("elapsed_sec"))
        rss = data.get("rss_mb", 0.0)
        health = data.get("health", "n/a")
        sub = data.get("subprogress")
        if sub:
            subtxt = (
                f"{sub.get('step')} {sub.get('done')}/{sub.get('total')} "
                f"({sub.get('pct')}%) eta={fmt_hms(sub.get('eta_sec'))}"
            )
        else:
            subtxt = "n/a"
        now = time.strftime("%H:%M:%S")
        print(
            f"[{now}] state={state} steps={done_steps}/{total_steps} elapsed={elapsed} "
            f"rss={rss:.1f}MB active={active_rules} health={health} sub={subtxt}",
            flush=True,
        )
        if state in {"completed", "failed", "stopped"}:
            return 0 if state == "completed" else 2
        iterations += 1
        if max_iterations > 0 and iterations >= max_iterations:
            return 0
        time.sleep(max(1.0, interval_sec))


def cmd_stop(pid_file: Path, workspace: Path, unlock: bool) -> int:
    if not pid_file.exists():
        print(f"pid file not found: {pid_file}")
        return 1
    raw = pid_file.read_text(encoding="utf-8").strip()
    if not raw.isdigit():
        print(f"invalid pid file: {pid_file}")
        return 1
    pid = int(raw)
    stopped = terminate_tree(pid)
    pid_file.unlink(missing_ok=True)
    if unlock:
        subprocess.run(["snakemake", "--unlock"], cwd=str(workspace), check=False)
    print(f"stop requested for pid={pid}; stopped={stopped}")
    return 0


def cmd_cleanup(workspace: Path, unlock: bool) -> int:
    ws = str(workspace.resolve()).lower()
    killed = set()
    roots = set()
    for p in psutil.process_iter(["pid", "name", "cmdline", "cwd"]):
        try:
            pid = int(p.info["pid"])
            if pid == psutil.Process().pid:
                continue
            name = str(p.info.get("name", "")).lower()
            cmdline = " ".join(p.info.get("cmdline") or []).lower()
            cwd = str(p.info.get("cwd") or "").lower()
            in_ws = (ws and ws in cwd) or (ws and ws in cmdline)
            is_snakemake = ("snakemake" in name) or ("snakemake" in cmdline)
            if in_ws and is_snakemake:
                roots.add(pid)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
    for pid in sorted(roots):
        terminate_tree(pid)
        killed.add(pid)
    if unlock:
        subprocess.run(["snakemake", "--unlock"], cwd=str(workspace), check=False)
    print(f"cleanup done; killed={len(killed)} pids={sorted(killed)}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="action", required=True)

    p_status = sub.add_parser("status")
    p_status.add_argument("--status-file", default="reports/run_status.json")

    p_watch = sub.add_parser("watch")
    p_watch.add_argument("--status-file", default="reports/run_status.json")
    p_watch.add_argument("--interval-sec", type=float, default=10.0)
    p_watch.add_argument("--max-iterations", type=int, default=0)

    p_stop = sub.add_parser("stop")
    p_stop.add_argument("--pid-file", default="reports/run_status.pid")
    p_stop.add_argument("--workspace", default=".")
    p_stop.add_argument("--unlock", action="store_true")

    p_cleanup = sub.add_parser("cleanup")
    p_cleanup.add_argument("--workspace", default=".")
    p_cleanup.add_argument("--unlock", action="store_true")

    args = parser.parse_args()

    if args.action == "status":
        return cmd_status(Path(args.status_file))
    if args.action == "watch":
        return cmd_watch(Path(args.status_file), float(args.interval_sec), int(args.max_iterations))
    if args.action == "stop":
        return cmd_stop(Path(args.pid_file), Path(args.workspace), bool(args.unlock))
    return cmd_cleanup(Path(args.workspace), bool(args.unlock))


if __name__ == "__main__":
    raise SystemExit(main())
