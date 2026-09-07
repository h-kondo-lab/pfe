#!/usr/bin/env python3
"""Build the docs copy of a Paper I manuscript from its sealed source.

The sealed Markdown manuscript is written for conversion to LaTeX, not for a
Markdown renderer, so three things do not survive a straight copy:

  * Figures appear only as `[Figure N] *caption*` placeholders. The images
    themselves live beside the LaTeX sources as PDFs and PNGs, and the PDF
    typesets them; a Markdown reader sees a caption with nothing above it.
  * The 43 reference entries sit on consecutive lines. Markdown joins
    consecutive lines into one paragraph, so the bibliography renders as a
    single run-on block instead of one entry per line as in the PDF.
  * Two footnote definitions are likewise adjacent and would be joined.

This script applies exactly those three repairs and nothing else: no wording,
number, table or equation is touched, so the docs copy stays faithful to the
sealed manuscript it is built from.

    python3 tools/build-manuscript.py \
        --source ~/Documents/pfe/papers/paper1/published/v4.1.1/full/supplementary/physics_from_existence_v4_1.md \
        --out docs/01-papers/physics_from_existence.md

Re-running it on an already-built file is a no-op, so it is safe to apply to
either the sealed source or the generated copy.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
from pathlib import Path

# Figure number -> image file, in the order the LaTeX sources include them.
FIGURES = {
    "1": "fig_potential.png",
    "2": "fig_pg23_v2.png",
}

FIGURE_LINE = re.compile(r"^\[Figure (\d+)\]\s*\*(.+)\*\s*$")
REFERENCE_LINE = re.compile(r"^\[\d+\] ")
FOOTNOTE_LINE = re.compile(r"^\[\^[^\]]+\]: ")


def _alt_text(caption: str) -> str:
    """First sentence of the caption, for the image's alt attribute."""
    for stop in ("。", ". "):
        head, sep, _ = caption.partition(stop)
        if sep:
            return head.strip()
    return caption.strip()


def build(text: str, img_dir: str) -> tuple[str, dict[str, int]]:
    counts = {"figures": 0, "references": 0, "footnotes": 0}
    lines = text.split("\n")
    out: list[str] = []

    for i, line in enumerate(lines):
        figure = FIGURE_LINE.match(line)
        if figure and figure.group(1) in FIGURES:
            number, caption = figure.group(1), figure.group(2)
            image = f"{img_dir}/{FIGURES[number]}"
            out.append(f"![{_alt_text(caption)}]({image})")
            out.append("")
            out.append(f"*Figure {number}. {caption}*")
            counts["figures"] += 1
            continue

        out.append(line)

        # Keep each reference and each footnote definition its own paragraph.
        following = lines[i + 1] if i + 1 < len(lines) else ""
        for pattern, key in ((REFERENCE_LINE, "references"), (FOOTNOTE_LINE, "footnotes")):
            if pattern.match(line) and pattern.match(following):
                out.append("")
                counts[key] += 1
                break

    return "\n".join(out), counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path, required=True, help="sealed manuscript Markdown")
    parser.add_argument("--out", type=Path, required=True, help="docs copy to write")
    parser.add_argument("--img-dir", default="img", help="image directory relative to --out (default: img)")
    parser.add_argument("--check", action="store_true", help="do not write; fail if --out is not up to date")
    args = parser.parse_args()

    if not args.source.exists():
        print(f"ERROR: source not found: {args.source}", file=sys.stderr)
        return 2

    source = args.source.read_text(encoding="utf-8")
    built, counts = build(source, args.img_dir)

    source_hash = hashlib.sha256(source.encode("utf-8")).hexdigest()
    built_hash = hashlib.sha256(built.encode("utf-8")).hexdigest()

    print(f"source {args.source}")
    print(f"  sha256 {source_hash}")
    print(f"  figures placed        : {counts['figures']}")
    print(f"  references separated  : {counts['references']}")
    print(f"  footnotes separated   : {counts['footnotes']}")
    print(f"built  {args.out}")
    print(f"  sha256 {built_hash}")

    if args.check:
        current = args.out.read_text(encoding="utf-8") if args.out.exists() else None
        if current == built:
            print("check: up to date")
            return 0
        print("check: OUT OF DATE — re-run without --check", file=sys.stderr)
        return 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(built, encoding="utf-8")
    print("written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
