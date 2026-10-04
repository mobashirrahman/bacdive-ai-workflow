"""Check actual model predictions against the published upstream example."""

import argparse
import json
from pathlib import Path

EXPECTED = {
    "Acidophilic": (False, 98.5),
    "Gram-positive": (True, 98.93),
    "Spore-forming": (False, 96.47),
    "Aerobic": (False, 93.91),
    "Anaerobic": (True, 63.54),
    "Thermophilic": (False, 99.77),
    "Psychrophilic": (False, 99.96),
    "Flagellated motility": (False, 99.85),
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("prediction", type=Path)
    args = parser.parse_args()
    data = json.loads(args.prediction.read_text())
    if set(data["predictions"]) != set(EXPECTED):
        raise SystemExit("Example does not contain all eight traits.")
    for trait, (value, confidence) in EXPECTED.items():
        actual = data["predictions"][trait]
        if actual["prediction"] is not value or abs(actual["confidence"] - confidence) > 0.011:
            raise SystemExit(f"Published example differs for {trait}: {actual}")
    print("All eight predictions agree with the published example.")


if __name__ == "__main__":
    main()
