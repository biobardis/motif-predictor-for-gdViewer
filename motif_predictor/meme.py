from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass(frozen=True)
class MemeMetadata:
    motif_widths: tuple[int, ...]

    @property
    def motif_count(self) -> int:
        return len(self.motif_widths)


_WIDTH_RE = re.compile(r"(?:^|\s)w\s*=\s*(\d+)(?:\s|$)")


def parse_meme_metadata(path: str) -> MemeMetadata:
    """Validate the MEME file and read motif widths for a short input summary."""
    widths: list[int] = []
    awaiting_matrix = False

    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue

            if line.startswith("MOTIF "):
                awaiting_matrix = True
                continue

            if awaiting_matrix and "letter-probability matrix:" in line.lower():
                match = _WIDTH_RE.search(line)
                if not match:
                    raise ValueError(
                        "Could not parse motif width from MEME letter-probability matrix line: "
                        + line
                    )
                width = int(match.group(1))
                if width < 1:
                    raise ValueError(f"Invalid motif width in MEME file: {width}")
                widths.append(width)
                awaiting_matrix = False

    if not widths:
        raise ValueError("No motif matrices with a parseable width were found in the MEME file")

    return MemeMetadata(motif_widths=tuple(widths))