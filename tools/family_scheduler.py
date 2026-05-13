#!/usr/bin/env python3
"""
family_scheduler.py — Tiered parallel scheduler for family_generate (C++ region solver)

Mirrors ScheduledRegionSolver.wl: same cache/status/tier/sector-worker architecture,
but drives the C++ binary instead of LIERegions`regionsBySectors.

Usage:
  python3 tools/family_scheduler.py <family_name> [options]

  <family_name> must correspond to families/<family_name>.json

Pipeline:
  1. Enumerate all subsectors from family JSON's topSector
  2. Pre-generate common IBP data via `family_generate --prepare`
  3. For each tier (Quick → Fast → Retry → Long):
     - Launch parallel `family_generate --sector X` workers
     - Track progress via cache/status.json (same format as MMA scheduler)
     - On timeout/error, re-queue to next tier
  4. After all sectors done, `family_generate --assemble` produces final IBPMat + RingData

State file format (cache/status.json):
  { "<sectorBin>": {"State":"done|running|queued|failed|trivial|cached", "Time":N, "Regions":N, "Tier":"..."} }
  Same as MMA scheduler's status.json, for live dashboard compatibility.

Examples:
  # Default: Quick=300s, Fast=3600s, Retry=36000s, 4 workers
  python3 tools/family_scheduler.py bub00

  # Custom tiers and parallelism
  python3 tools/family_scheduler.py TB123 --tiers "Quick=120 Fast=1200 Retry=86400" --workers 6

  # Resume from existing cache
  python3 tools/family_scheduler.py NP322 --resume
"""

import json, os, sys, time, subprocess, signal, glob, argparse
from pathlib import Path
from collections import OrderedDict

# ── Default configuration (mirrors MMA's $DefaultScheduleConfig) ──────────

DEFAULT_TIERS = OrderedDict([
    ("Quick", 300),
    ("Fast",  3600),
    ("Retry", 36000),
])

DEFAULT_WORKERS = 4
POLL_INTERVAL  = 3  # seconds between worker reaps

# ── Scheduler ────────────────────────────────────────────────────────────

class FamilyScheduler:
    def __init__(self, family_name, tiers=None, max_workers=None, resume=False, cache_dir=None, output_dir=None):
        proj_root = self._find_project_root()
        self.family_name = family_name
        self.family_json = proj_root / "families" / f"{family_name}.json"
        self.family_generate = proj_root / "family_generate"

        if not self.family_json.exists():
            raise FileNotFoundError(f"Family config not found: {self.family_json}")
        if not self.family_generate.exists():
            raise FileNotFoundError(f"Binary not found: {self.family_generate}. Run 'make family_generate' first.")

        self.tiers = tiers or DEFAULT_TIERS
        self.max_workers = max_workers or DEFAULT_WORKERS

        self.cache_dir = Path(cache_dir) if cache_dir else proj_root / "cache" / f"{family_name}_scheduled"
        self.output_dir = Path(output_dir) if output_dir else self.cache_dir / "output"
        self.sector_cache_dir = self.cache_dir / "SectorCache"

        self.status_file = self.cache_dir / "status.json"
        self.lock_file = self.cache_dir / "instance.lock"
        self.log_dir = self.cache_dir / "logs"

        self.status = {}  # sectorKey -> {State, Time, Regions, Tier}
        self.resume = resume
        self.running_workers = {}  # pid -> {key, tier_idx, start_time, timeout}

    # ── helpers ──────────────────────────────────────────────────────

    @staticmethod
    def _find_project_root():
        """Find project root by walking up from CWD looking for families/."""
        p = Path.cwd().resolve()
        for _ in range(10):
            if (p / "families").is_dir() and (p / "family_generate").exists():
                return p
            if (p / ".git").exists() and (p / "families").is_dir():
                return p
            p = p.parent
        # fallback: just use CWD
        return Path.cwd().resolve()

    @staticmethod
    def _sector_to_str(sector_bits):
        return "".join(str(b) for b in sector_bits)

    def _log(self, msg):
        print(f"[SCHEDULER] {msg}", flush=True)

    # ── cache/status management ──────────────────────────────────────

    def _load_status(self):
        if self.status_file.exists():
            try:
                return json.loads(self.status_file.read_text())
            except (json.JSONDecodeError, OSError):
                return {}
        return {}

    def _save_status(self):
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        tmp = self.status_file.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(self.status, indent=2, ensure_ascii=False))
        tmp.replace(self.status_file)

    def _acquire_lock(self):
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        if self.lock_file.exists():
            pid_str = self.lock_file.read_text().strip()
            if pid_str:
                try:
                    pid = int(pid_str)
                    os.kill(pid, 0)
                    self._log(f"[LOCK] Another instance running (PID {pid}). Use --resume or kill it.")
                    return False
                except (ProcessLookupError, ValueError):
                    self._log(f"[LOCK] Stale lock (PID {pid_str}). Removing.")
                    self.lock_file.unlink(missing_ok=True)
        self.lock_file.write_text(str(os.getpid()))
        return True

    def _release_lock(self):
        self.lock_file.unlink(missing_ok=True)

    # ── sector enumeration ───────────────────────────────────────────

    def _enumerate_sectors(self):
        """Enumerate all subsectors from family JSON's topSector."""
        config = json.loads(self.family_json.read_text())
        top_sector = config.get("topSector", [])
        ne = len(top_sector)
        # Count active propagator positions
        active = [i for i, v in enumerate(top_sector) if v == 1]
        # Minimum subsector size = max(1, nLoop) (always at least 1)
        # We don't know nLoop from JSON alone, get it from config
        n_loop = len(config.get("loopMomenta", []))
        min_size = max(1, n_loop)

        # Generate all subsets of active positions with size >= min_size
        sectors = []
        from itertools import combinations
        for size in range(min_size, len(active) + 1):
            for combo in combinations(active, size):
                sector = [0] * ne
                for idx in combo:
                    sector[idx] = 1
                sectors.append(sector)

        # Sort: by weight ascending, then by binary value descending (MMA convention)
        sectors.sort(key=lambda s: (sum(s), -int("".join(str(b) for b in s), 2)))
        return config, sectors

    # ── worker management ────────────────────────────────────────────

    def _launch_worker(self, sector_key, sector, tier_idx, timeout):
        """Launch family_generate --sector for one sector."""
        self.sector_cache_dir.mkdir(parents=True, exist_ok=True)
        self.log_dir.mkdir(parents=True, exist_ok=True)
        log_file = self.log_dir / f"worker_{sector_key}.log"

        cmd = [
            str(self.family_generate), str(self.family_json),
            "--sector", sector_key,
            "--cache-dir", str(self.sector_cache_dir),
        ]
        self._log(f"  CMD: {' '.join(cmd)}")
        with open(log_file, "w") as lf:
            proc = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT)

        self.running_workers[proc.pid] = {
            "key": sector_key,
            "tier_idx": tier_idx,
            "start_time": time.time(),
            "timeout": timeout,
            "proc": proc,
        }

    def _reap_workers(self):
        """Check all running workers; return lists of completed/failed/timeout."""
        done = []
        failed = []
        timeout = []
        to_kill = []

        for pid, info in list(self.running_workers.items()):
            proc = info["proc"]
            ret = proc.poll()

            elapsed = time.time() - info["start_time"]

            if ret is not None:
                # Process exited
                sector_key = info["key"]
                # Check for result (.bin file should exist)
                result_file = self.sector_cache_dir / f"IBPMat_{self.family_name}_sector_{sector_key}.bin"
                has_result = result_file.exists()
                # Also check for RingData
                ring_file = self.sector_cache_dir / f"RingData_{self.family_name}_sector_{sector_key}.bin"
                has_ring = ring_file.exists()
                if has_result and has_ring:
                    done.append((sector_key, info["tier_idx"]))
                    self._log(f"  [DONE] {sector_key} in {elapsed:.1f}s (exit={ret})")
                else:
                    # Even with exit=0, if no result file, it was a trivial sector (0 regions)
                    # Mark as done (trivial) since empty output was written
                    trivial_file = self.sector_cache_dir / f"IBPMat_{self.family_name}_sector_{sector_key}.bin"
                    if ret == 0 and trivial_file.exists():
                        done.append((sector_key, info["tier_idx"]))
                        self._log(f"  [TRIVIAL] {sector_key} in {elapsed:.1f}s (0 regions)")
                    else:
                        failed.append((sector_key, info["tier_idx"]))
                        self._log(f"  [FAIL] {sector_key} exit={ret} in {elapsed:.1f}s (no result)")
                del self.running_workers[pid]

            elif elapsed > info["timeout"]:
                # Timeout: kill the process
                sector_key = info["key"]
                self._log(f"  [TIMEOUT] {sector_key} after {elapsed:.1f}s (limit={info['timeout']}s)")
                proc.kill()
                timeout.append((sector_key, info["tier_idx"]))
                to_kill.append(pid)

        for pid in to_kill:
            del self.running_workers[pid]

        # Wait briefly for killed processes
        for pid in to_kill:
            try:
                os.waitpid(pid, 0)
            except OSError:
                pass

        return done, failed, timeout

    # ── main scheduling loop ─────────────────────────────────────────

    def run(self):
        tier_names = list(self.tiers.keys())
        n_tiers = len(tier_names)

        if not self._acquire_lock():
            return 1

        try:
            # ── 1. Enumerate sectors ────────────────────────────────
            config, sectors = self._enumerate_sectors()
            total = len(sectors)
            # Use the canonical name from JSON, not CLI arg
            self.family_name = config.get("name", self.family_name)
            self._log(f"Family: {self.family_name}, {total} sectors, {n_tiers} tiers: {' → '.join(tier_names)}")

            # ── 2. Initialize status ─────────────────────────────────
            if self.resume and self.status_file.exists():
                self.status = self._load_status()
                self._log(f"Resumed: {sum(1 for v in self.status.values() if v.get('State') in ('done','trivial','cached'))}/{total} already completed")
            else:
                self.status = {}
                for sector in sectors:
                    key = self._sector_to_str(sector)
                    self.status[key] = {"State": "pending", "Time": 0, "Regions": 0}
                self._save_status()

            # ── 3. Build tier queues ─────────────────────────────────
            queues = [[] for _ in range(n_tiers)]
            for sector in sectors:
                key = self._sector_to_str(sector)
                st = self.status.get(key, {}).get("State", "pending")

                if st in ("done", "trivial", "cached"):
                    continue

                # State-aware routing
                routed = False
                for ti in range(n_tiers - 1):
                    if st == f"{tier_names[ti]}-timeout" or st == f"{tier_names[ti]}-error":
                        queues[ti + 1].append(sector)
                        routed = True
                        break
                if not routed:
                    if st in ("error", "failed"):
                        queues[n_tiers - 1].append(sector)
                    elif st == "running":
                        queues[0].append(sector)
                    else:
                        queues[0].append(sector)

            pending_total = sum(len(q) for q in queues)
            self._log(f"Pending: {pending_total}/{total}")
            for ti in range(n_tiers):
                if queues[ti]:
                    self._log(f"  {tier_names[ti]}: {len(queues[ti])} sectors")

            if pending_total == 0:
                self._log("All sectors already completed. Running assembly...")
                self._assemble(config)
                return 0

            # ── 4. Tiered parallel execution ────────────────────────
            self.sector_cache_dir.mkdir(parents=True, exist_ok=True)
            self.log_dir.mkdir(parents=True, exist_ok=True)

            # Track sectors remaining per queue
            active = {ti: 0 for ti in range(n_tiers)}

            while True:
                # Reap completed workers
                done, failed_workers, timeout_workers = self._reap_workers()

                for key, ti in done:
                    active[ti] -= 1
                    # Reload status (worker may have written its own state)
                    self.status = self._load_status()
                    st = self.status.get(key, {}).get("State", "")
                    if st in ("pending", "running", "queued", ""):
                        self.status[key] = {"State": "done", "Time": 0, "Regions": 0}
                        self._save_status()

                for key, ti in failed_workers:
                    active[ti] -= 1
                    self.status[key] = {"State": "failed", "Time": 0, "Regions": 0}
                    self._save_status()

                for key, ti in timeout_workers:
                    active[ti] -= 1
                    self.status[key] = {"State": f"{tier_names[ti]}-timeout", "Time": 0, "Regions": 0}
                    self._save_status()
                    if ti < n_tiers - 1:
                        # Find the sector bits for re-queue
                        for sector in sectors:
                            if self._sector_to_str(sector) == key:
                                queues[ti + 1].append(sector)
                                break
                    else:
                        # Last tier timeout → mark as failed
                        self.status[key] = {"State": "failed", "Time": 0, "Regions": 0}
                        self._save_status()

                # Launch new workers (fill up to max_workers)
                total_active = sum(active.values())
                can_launch = self.max_workers - total_active

                while can_launch > 0:
                    launched = False
                    # Priority: lower tiers first (Quick before Fast before Retry)
                    for ti in range(n_tiers):
                        if can_launch <= 0:
                            break
                        if queues[ti]:
                            sector = queues[ti].pop(0)
                            key = self._sector_to_str(sector)
                            timeout = list(self.tiers.values())[ti]
                            self._launch_worker(key, sector, ti, timeout)
                            self.status[key] = {
                                "State": "running",
                                "Time": 0,
                                "Regions": 0,
                                "Tier": tier_names[ti],
                            }
                            active[ti] += 1
                            total_active += 1
                            can_launch -= 1
                            launched = True
                            self._log(f"  [LAUNCH] {key} → {tier_names[ti]} tier ({timeout}s) [PID {list(self.running_workers.keys())[-1]}]")
                    if not launched:
                        break

                # Save status for live dashboard
                self._save_status()

                # Check termination
                all_queues_empty = all(len(q) == 0 for q in queues)
                if all_queues_empty and total_active == 0:
                    self._log("All workers finished.")
                    break

                time.sleep(POLL_INTERVAL)

            # ── 5. Mark remaining non-done sectors as failed ─────────
            self.status = self._load_status()
            for sector in sectors:
                key = self._sector_to_str(sector)
                st = self.status.get(key, {}).get("State", "pending")
                if st not in ("done", "trivial", "cached"):
                    self.status[key] = {"State": "failed", "Regions": 0, "Time": 0}
            self._save_status()

            failed_count = sum(1 for v in self.status.values() if v.get("State") == "failed")
            if failed_count:
                self._log(f"WARNING: {failed_count}/{total} sectors failed")
            else:
                self._log(f"All {total} sectors completed successfully!")

            # ── 6. Assemble final output ─────────────────────────────
            self._assemble(config)

        finally:
            # Kill any remaining workers
            for pid, info in list(self.running_workers.items()):
                try:
                    info["proc"].kill()
                except OSError:
                    pass
            self._release_lock()

        return 0

    # ── assembly ────────────────────────────────────────────────────

    def _assemble(self, config):
        """Call family_generate --merge-sectors to produce final IBPMat + RingData."""
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Check if we have partial sector files to merge
        partial_files = list(self.sector_cache_dir.glob("IBPMat_*_sector_*.bin"))
        if not partial_files:
            self._log(f"No partial sector files found in {self.sector_cache_dir}")
            # Fall back to normal mode
            cmd = [
                str(self.family_generate), str(self.family_json),
                "--output", str(self.output_dir),
            ]
        else:
            cmd = [
                str(self.family_generate), str(self.family_json),
                "--merge-sectors", str(self.sector_cache_dir),
                "--output", str(self.output_dir),
            ]

        self._log(f"Assembling final output...")
        self._log(f"  {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            self._log("Assembly complete!")
            for line in result.stdout.strip().split("\n"):
                if line.strip():
                    self._log(f"  | {line}")
        else:
            self._log(f"Assembly FAILED (exit={result.returncode})")
            if result.stderr:
                for line in result.stderr.strip().split("\n"):
                    self._log(f"  ERR: {line}")
        # Merge stderr too
        if result.stderr:
            for line in result.stderr.strip().split("\n"):
                if line.strip():
                    self._log(f"  | {line}")

        # Export summary for task-manager dashboard
        summary_path = self.cache_dir / "status.json"
        if summary_path.exists():
            try:
                status = json.loads(summary_path.read_text())
                total = len(status)
                done = sum(1 for v in status.values() if v.get("State") in ("done", "trivial", "cached"))
                failed = sum(1 for v in status.values() if v.get("State") == "failed")
                self._log(f"\nSummary: {done}/{total} sectors completed, {failed} failed")
            except (json.JSONDecodeError, OSError):
                pass


# ── CLI ────────────────────────────────────────────────────────────────

def parse_tiers(tier_str):
    """Parse 'Quick=300 Fast=3600 Retry=36000' into OrderedDict."""
    tiers = OrderedDict()
    for part in tier_str.split():
        if "=" in part:
            name, val = part.split("=", 1)
            tiers[name.strip()] = int(val)
    return tiers


def main():
    parser = argparse.ArgumentParser(
        description="Tiered parallel scheduler for family_generate (C++ region solver)"
    )
    parser.add_argument("family", help="Family name (matches families/<name>.json)")
    parser.add_argument("--tiers", default=None,
                        help='Tier spec: "Quick=300 Fast=3600 Retry=36000"')
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                        help=f"Max parallel workers (default: {DEFAULT_WORKERS})")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from existing cache/status.json")
    parser.add_argument("--cache-dir", default=None,
                        help="Override cache directory")
    parser.add_argument("--output-dir", default=None,
                        help="Override output directory")
    args = parser.parse_args()

    tiers = parse_tiers(args.tiers) if args.tiers else DEFAULT_TIERS

    try:
        scheduler = FamilyScheduler(
            family_name=args.family,
            tiers=tiers,
            max_workers=args.workers,
            resume=args.resume,
            cache_dir=args.cache_dir,
            output_dir=args.output_dir,
        )
        return scheduler.run()
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    exit(main())
