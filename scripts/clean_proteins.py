#!/usr/bin/env python3
"""Strip stray asterisks from a Prodigal protein FASTA in place.

Replaces the GNU-only `sed -i 's/\\*//g'` this workflow previously used, which
fails on BSD and macOS sed where `-i` requires an explicit backup suffix.
"""

import argparse
import sys
from pathlib import Path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("proteins", type=Path)
    args = parser.parse_args(argv)
    path = Path(args.proteins)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as error:
        parser.exit(1, f"Error: cannot read {path}: {error}\n")
    if "*" not in text:
        return 0
    temporary = path.with_name(f".{path.name}.clean")
    try:
        temporary.write_text(text.replace("*", ""), encoding="utf-8")
        temporary.replace(path)
    except OSError as error:
        temporary.unlink(missing_ok=True)
        parser.exit(1, f"Error: cannot write {path}: {error}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
