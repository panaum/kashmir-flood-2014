"""Draw the decorative divider band for the website (docs/img/band.svg).

Single-colour line art, 1200 x 44. On the left a khatamband-style lattice
(crossing diagonals with eight-pointed stars at their centres); in the middle a
medallion where the star sits inside a circle; on the right Cosmati-style
marble inlay (interlaced roundels with lozenges between them). Decoration only.

Usage:  python3 scripts/20_draw_band.py
"""

import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
W, H, CY = 1200, 44, 22
TOP, BOT = 3, 41
STROKE = "#8A3A1B"        # chinar / sienna
ACCENT = "#C9952E"        # saffron / Naples, decorative only
LEFT_END, MID, RIGHT_START = 552, 600, 648
TILE = 36


def star(cx, cy, r):
    """Eight-pointed star: two squares, one turned 45 degrees."""
    a = [(cx - r, cy - r), (cx + r, cy - r), (cx + r, cy + r), (cx - r, cy + r)]
    d = r * math.sqrt(2)
    b = [(cx, cy - d), (cx + d, cy), (cx, cy + d), (cx - d, cy)]
    return [poly(a), poly(b)]


def poly(pts):
    return "<path d='M" + " L".join(f"{x:.2f} {y:.2f}" for x, y in pts) + " Z'/>"


def main():
    parts = [f"<path d='M0 {TOP} H{W} M0 {BOT} H{W}'/>"]

    # Khatamband lattice: diagonals rail to rail, stars where they cross.
    x = -TILE / 2
    while x < LEFT_END - TILE / 2:
        x0, x1 = x, x + TILE
        parts.append(f"<path d='M{x0:.1f} {TOP} L{x1:.1f} {BOT} M{x0:.1f} {BOT} L{x1:.1f} {TOP}'/>")
        cx = (x0 + x1) / 2
        if cx > 4:
            parts += star(cx, CY, 5.2)
            parts.append(f"<circle cx='{cx:.1f}' cy='{CY}' r='1.6'/>")
        x += TILE

    # Medallion: star inside a circle, flanked by small roundels.
    parts.append(f"<circle cx='{MID}' cy='{CY}' r='16'/>")
    parts.append(f"<circle cx='{MID}' cy='{CY}' r='13'/>")
    parts += star(MID, CY, 7.2)
    for dx in (-30, 30):
        parts.append(f"<circle cx='{MID + dx}' cy='{CY}' r='5'/>")
    parts.append(f"<path d='M{LEFT_END} {CY} H{MID - 35} M{MID + 35} {CY} H{RIGHT_START}'/>")

    # Cosmati inlay: roundels with an interlacing band, lozenges between.
    step, r = 34, 12
    cx = RIGHT_START + step / 2
    centres = []
    while cx < W + step:
        centres.append(cx)
        cx += step
    for c in centres:
        parts.append(f"<circle cx='{c:.1f}' cy='{CY}' r='{r}'/>")
        parts.append(f"<circle cx='{c:.1f}' cy='{CY}' r='{r - 4.5}'/>")
    for a, b in zip(centres, centres[1:]):
        m = (a + b) / 2
        parts.append(f"<path d='M{a:.1f} {CY - r} C{m:.1f} {CY - r} {m:.1f} {CY + r} {b:.1f} {CY + r}'/>")
        parts.append(f"<path d='M{a:.1f} {CY + r} C{m:.1f} {CY + r} {m:.1f} {CY - r} {b:.1f} {CY - r}'/>")
        parts.append(poly([(m, TOP + 2.5), (m + 4.5, TOP + 7), (m, TOP + 11.5), (m - 4.5, TOP + 7)]))
        parts.append(poly([(m, BOT - 2.5), (m + 4.5, BOT - 7), (m, BOT - 11.5), (m - 4.5, BOT - 7)]))

    svg = (f"<svg xmlns='http://www.w3.org/2000/svg' width='{W}' height='{H}' viewBox='0 0 {W} {H}'>"
           f"<g fill='none' stroke='{STROKE}' stroke-width='0.8' stroke-linejoin='round'>"
           + "".join(parts) + "</g>"
           f"<circle cx='{MID}' cy='{CY}' r='1.8' fill='{ACCENT}'/></svg>")
    out = ROOT / "docs" / "img" / "band.svg"
    out.write_text(svg)
    print(f"wrote {out.relative_to(ROOT)} ({len(svg)} bytes)")


if __name__ == "__main__":
    main()
