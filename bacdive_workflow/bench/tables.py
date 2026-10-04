"""Final table assembly and BENCHMARK.md marker regeneration."""

import re
import shutil
from pathlib import Path

MARKER = "<!-- benchmark:{name} -->"


def copy_final(tables, destination):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    for table in tables:
        shutil.copy2(table, destination / Path(table).name)


def markdown_table(path, columns=None, limit=50):
    import csv

    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if columns is None:
        columns = rows[0].keys() if rows else []
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in rows[:limit]:
        lines.append("| " + " | ".join(str(row.get(c, "")) for c in columns) + " |")
    return "\n".join(lines) + "\n"


def regenerate_markdown(document, tables):
    """Replace marker blocks in BENCHMARK.md; return True when content changed."""
    text = Path(document).read_text()
    changed = False
    for name, table in tables.items():
        pattern = re.compile(
            rf"{re.escape(MARKER.format(name=name))}\n```tab\n.*?```",
            re.DOTALL,
        )
        replacement = f"{MARKER.format(name=name)}\n```tab\n{markdown_table(table)}```"
        updated, count = pattern.subn(replacement, text)
        if count:
            text, changed = updated, True
    if changed:
        Path(document).write_text(text)
    return changed
