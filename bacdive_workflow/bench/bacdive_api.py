"""BacDive API harvest with rate limiting, retries, and raw-JSON caching."""

import json
import time
import urllib.request
from pathlib import Path

from bacdive_workflow.common import WorkflowError

BASE = "https://api.bacdive.dsmz.de/v2"
MAX_IDS_PER_REQUEST = 100
REQUESTS_PER_SECOND = 0.5


def fetch_batch(ids, cache_dir, opener=None):
    """Fetch one batch of BacDive IDs; reuse the cached raw JSON when present."""
    key = "-".join(str(i) for i in ids)
    cached = Path(cache_dir) / f"{key}.json"
    if cached.is_file():
        return json.loads(cached.read_text())
    url = f"{BASE}/fetch/{';'.join(str(i) for i in ids)}"
    attempts = 0
    while True:
        try:
            if opener is not None:
                payload = opener(url)
            else:
                with urllib.request.urlopen(url, timeout=60) as response:
                    payload = response.read().decode()
            break
        except OSError as error:
            attempts += 1
            if attempts >= 5:
                raise WorkflowError(f"BacDive fetch failed for {key}: {error}") from error
            time.sleep(2**attempts)
    cached.parent.mkdir(parents=True, exist_ok=True)
    cached.write_text(payload)
    time.sleep(1 / REQUESTS_PER_SECOND)
    return json.loads(payload)


def harvest(ids, cache_dir, batch_size=MAX_IDS_PER_REQUEST, opener=None):
    """Yield parsed records for every ID, caching each raw batch response."""
    for start in range(0, len(ids), batch_size):
        batch = ids[start : start + batch_size]
        yield from fetch_batch(batch, cache_dir, opener=opener)
