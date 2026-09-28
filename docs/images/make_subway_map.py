"""Generate the nf-core style subway map SVGs (light + dark) used in the README.

Same drawing style as mcbla-bulkatac-pipe's docs/images/make_subway_map.py.

Usage: python docs/images/make_subway_map.py docs/images
"""
import sys
from html import escape
from pathlib import Path

W, H = 1520, 560
Y = 220  # trunk
Y_ENR, Y_QC = 110, 410

LINES = {
    "core":  ("#24B064", "Processing: FASTQ to normalized counts"),
    "qc":    ("#F2B138", "QC"),
    "diff":  ("#1F6FEB", "Differential expression (edgeR)"),
    "enr":   ("#E8702A", "Enrichment (GSEApy)"),
}
DASHED = set()

X0, X_BAM, X_COUNTS, X_DE = 90, 630, 1090, 1250
PATHS = {
    "qc":   [(X_BAM, Y), (X_BAM, Y_QC - 60), (X_BAM + 60, Y_QC), (1030, Y_QC)],
    "enr":  [(X_DE, Y), (X_DE + (Y - Y_ENR), Y_ENR), (1440, Y_ENR)],
    "diff": [(X_COUNTS, Y), (1440, Y)],
    "core": [(X0, Y), (X_COUNTS, Y)],
}
DRAW_ORDER = ["qc", "enr", "diff", "core"]
LEGEND = ["core", "qc", "diff", "enr"]
LEGEND_COLS, LEGEND_DX = 4, 350

# stations: (x, y, label, label side, kind); kind: stop | hub | start | end
S = [
    (X0, Y, "FASTQ", "below", "start"),
    (200, Y, "fastp\npoly-X + quality\ntrim", "below", "stop"),
    (310, Y, "STAR\nsplice-aware", "below", "stop"),
    (420, Y, "samtools\nsort + index", "below", "stop"),
    (530, Y, "UMICollapse\nUMI dedup", "below", "stop"),
    (X_BAM, Y, "Final BAMs", "above", "hub"),
    (750, Y, "featureCounts\ngene bodies, 1/NH", "above", "stop"),
    (870, Y, "edgeR TMM\nnormalization", "above", "stop"),
    (985, Y, "Correlation\n+ PCA", "above", "stop"),
    (X_COUNTS, Y, "Normalized\ncounts", "below", "hub"),
    # QC
    (740, Y_QC, "RSeQC\nstrand, read\ndistribution", "below", "stop"),
    (840, Y_QC, "Gene body\ncoverage\n(housekeeping)", "below", "stop"),
    (940, Y_QC, "Qualimap\nRNA-seq QC", "below", "stop"),
    (1030, Y_QC, "MultiQC", "above", "end"),
    # differential expression
    (X_DE, Y, "edgeR QL\nper contrast", "below", "stop"),
    (1440, Y, "Tables, MA +\nvolcano plots", "below", "end"),
    # enrichment
    (1440, Y_ENR, "GSEApy prerank\nMSigDB Hallmark", "above", "end"),
]

STAGES = [(40, X_BAM + 30, "1", "Preprocessing"), (X_BAM + 30, X_COUNTS + 30, "2", "Counts & QC"),
          (X_COUNTS + 30, W - 20, "3", "Differential expression")]


def svg(theme):
    fg = "#1F2328" if theme == "light" else "#E6EDF3"
    muted = "#59636E" if theme == "light" else "#9198A1"
    band = "#F6F8FA" if theme == "light" else "#161B22"
    stfill = "#FFFFFF"
    stroke = "#1F2328"
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
         f'font-family="Helvetica Neue, Helvetica, Arial, sans-serif">',
         '<title>mcbla-plasmidsaurus-pipe subway map</title>']
    # stage bands
    for x0, x1, n, name in STAGES:
        o.append(f'<rect x="{x0}" y="20" width="{x1 - x0 - 8}" height="{H - 100}" rx="14" fill="{band}"/>')
        o.append(f'<circle cx="{x0 + 22}" cy="44" r="13" fill="{fg}"/>'
                 f'<text x="{x0 + 22}" y="49" font-size="14" font-weight="700" fill="{band}" text-anchor="middle">{n}</text>'
                 f'<text x="{x0 + 42}" y="49" font-size="15" font-weight="700" fill="{fg}">{escape(name)}</text>')
    # lines
    for key in DRAW_ORDER:
        pts = " ".join(f"{x},{y}" for x, y in PATHS[key])
        cap = ' stroke-dasharray="16 8" stroke-linecap="butt"' if key in DASHED else ' stroke-linecap="round"'
        o.append(f'<polyline points="{pts}" fill="none" stroke="{LINES[key][0]}" stroke-width="9" '
                 f'stroke-linejoin="round"{cap}/>')
    # stations
    for x, y, label, side, kind in S:
        if kind == "hub":
            o.append(f'<rect x="{x - 13}" y="{y - 13}" width="26" height="26" rx="13" fill="{stfill}" stroke="{stroke}" stroke-width="3.5"/>')
        elif kind in ("start", "end"):
            o.append(f'<rect x="{x - 10}" y="{y - 12}" width="20" height="24" rx="3" fill="{stfill}" stroke="{stroke}" stroke-width="3"/>'
                     f'<line x1="{x - 5}" y1="{y - 4}" x2="{x + 5}" y2="{y - 4}" stroke="{stroke}" stroke-width="2"/>'
                     f'<line x1="{x - 5}" y1="{y + 2}" x2="{x + 5}" y2="{y + 2}" stroke="{stroke}" stroke-width="2"/>')
        else:
            o.append(f'<circle cx="{x}" cy="{y}" r="9" fill="{stfill}" stroke="{stroke}" stroke-width="3"/>')
        lines = label.split("\n")
        lh = 15
        if side == "below":
            ys, anchor, lx = [y + 32 + i * lh for i in range(len(lines))], "middle", x
        elif side == "above":
            ys, anchor, lx = [y - 22 - (len(lines) - 1 - i) * lh for i in range(len(lines))], "middle", x
        else:
            ys, anchor, lx = [y + 5 + (i - (len(lines) - 1) / 2) * lh for i in range(len(lines))], "start", x + 18
        for i, (t, ty) in enumerate(zip(lines, ys)):
            weight = "700" if i == 0 else "400"
            col = fg if i == 0 else muted
            size = 13 if i == 0 else 11.5
            o.append(f'<text x="{lx}" y="{ty:.0f}" font-size="{size}" font-weight="{weight}" fill="{col}" text-anchor="{anchor}">{escape(t)}</text>')
    # legend
    ly = H - 58
    for i, key in enumerate(LEGEND):
        col, name = LINES[key]
        cx = 40 + (i % LEGEND_COLS) * LEGEND_DX
        cy = ly + (i // LEGEND_COLS) * 28
        cap = ' stroke-dasharray="10 5" stroke-linecap="butt"' if key in DASHED else ' stroke-linecap="round"'
        o.append(f'<line x1="{cx}" y1="{cy}" x2="{cx + 36}" y2="{cy}" stroke="{col}" stroke-width="8"{cap}/>'
                 f'<text x="{cx + 50}" y="{cy + 5}" font-size="13" fill="{fg}">{escape(name)}</text>')
    o.append("</svg>")
    return "\n".join(o)


out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
for t in ("light", "dark"):
    (out / f"subway_map_{t}.svg").write_text(svg(t) + "\n")
