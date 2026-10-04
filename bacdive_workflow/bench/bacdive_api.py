"""BacDive API harvest with rate limiting, retries, and raw-JSON caching.

Command:
    python -m bacdive_workflow.bench.bacdive_api harvest
        --training-labels <training_data_labels.csv>
        --cache-dir <bacdive_cache>
        --manifest <harvest.json>
"""

import argparse
import csv
import datetime
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

from bacdive_workflow.common import WorkflowError

BASE = "https://api.bacdive.dsmz.de/v2"
REQUESTS_PER_SECOND = 2.0
EMPTY_STOP_STREAK = 20


def _get(url, attempts=5):
    last = None
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                return response.read().decode()
        except urllib.error.HTTPError as error:
            if error.code == 400:
                raise WorkflowError(f"BacDive rejected the request (HTTP 400): {url}") from error
            last = error
            time.sleep(2**attempt)
        except OSError as error:
            last = error
            time.sleep(2**attempt)
    raise WorkflowError(f"BacDive request failed for {url}: {last}")


def probe_max_batch(opener=None):
    """Confirm the real per-request ID cap (100 works, 101 fails as of 2026-10-04)."""
    get = opener or _get
    get(f"{BASE}/fetch/{';'.join(str(i) for i in range(1, 101))}")
    try:
        get(f"{BASE}/fetch/{';'.join(str(i) for i in range(1, 102))}")
    except (OSError, WorkflowError):
        return 100
    return 101


def fetch_ids(ids, cache_path):
    """Fetch one batch, reusing the cached raw JSON when present."""
    cache_path = Path(cache_path)
    if cache_path.is_file():
        return json.loads(cache_path.read_text())
    url = f"{BASE}/fetch/{';'.join(str(i) for i in ids)}"
    try:
        payload = _get(url)
    finally:
        time.sleep(1 / REQUESTS_PER_SECOND)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(payload)
    return json.loads(payload)


def harvest(training_labels, cache_dir, manifest_path):
    cache_dir = Path(cache_dir)
    with open(training_labels, newline="", encoding="utf-8") as handle:
        training_ids = sorted({int(r["ID_strains"]) for r in csv.DictReader(handle)})
    batch_size = probe_max_batch()
    harvested = {"training": {"found": [], "missing": []}, "candidates": []}
    empty_streak = 0

    for start in range(0, len(training_ids), batch_size):
        chunk = training_ids[start : start + batch_size]
        cached = cache_dir / "train" / f"batch-{start // batch_size:04d}.json"
        payload = fetch_ids(chunk, cached)
        results = payload.get("results", {})
        for bacdive_id in chunk:
            if str(bacdive_id) in results:
                harvested["training"]["found"].append(bacdive_id)
            else:
                harvested["training"]["missing"].append(bacdive_id)

    highest = 0
    start = 1
    while True:
        chunk = list(range(start, start + batch_size))
        cached = cache_dir / "candidates" / f"batch-{start:07d}.json"
        payload = fetch_ids(chunk, cached)
        results = payload.get("results", {})
        if not results:
            empty_streak += 1
            if empty_streak >= EMPTY_STOP_STREAK:
                break
        else:
            empty_streak = 0
            found = sorted(int(i) for i in results)
            harvested["candidates"].extend(found)
            highest = max(highest, max(found))
        start += batch_size

    manifest = {
        "harvest_date": datetime.date.today().isoformat(),
        "batch_size": batch_size,
        "n_training_ids": len(training_ids),
        "n_training_found": len(harvested["training"]["found"]),
        "training_missing": harvested["training"]["missing"],
        "n_candidate_ids": len(harvested["candidates"]),
        "highest_candidate_id": highest,
        "empty_stop_streak": EMPTY_STOP_STREAK,
    }
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    index_path = manifest_path.with_suffix(".ids.json")
    index_path.write_text(json.dumps(harvested, indent=2) + "\n")
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description="Harvest BacDive records with caching.")
    parser.add_argument("command", choices=["harvest"])
    parser.add_argument("--training-labels", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        manifest = harvest(args.training_labels, args.cache_dir, args.manifest)
        print(
            f"Harvested {manifest['n_training_found']} training and {manifest['n_candidate_ids']} candidate records."
        )
        return 0
    except (WorkflowError, OSError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
