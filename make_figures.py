"""
make_figures.py — Hollenbach Lab Demo
Generates two SVG figures from real IMGT/HLA data.
"""

import csv, math, os

OUT = os.path.dirname(__file__)

# ── Load data ─────────────────────────────────────────────────────────────────
rows = []
with open(os.path.join(OUT, "allele_counts.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        rows.append((row["locus"], int(row["allele_count"]), row["class"]))

growth = []
with open(os.path.join(OUT, "allele_growth.tsv")) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        growth.append((int(row["year"]), int(row["total_alleles"])))

# ── Figure 1: Per-locus allele count bar chart ────────────────────────────────
# Show the main loci only (exclude low-count non-classical for clarity)
SHOW_LOCI = [
    ("B",    "I"),
    ("C",    "I"),
    ("A",    "I"),
    ("DPB1", "II"),
    ("DQB1", "II"),
    ("DRB1", "II"),
    ("DQA1", "II"),
    ("DPA1", "II"),
    ("MICA", "other"),
    ("MICB", "other"),
    ("DRB3", "II"),
    ("DRB4", "II"),
    ("E",    "I"),
    ("G",    "I"),
    ("F",    "I"),
]
count_map = {r[0]: r[1] for r in rows}
bar_data = [(locus, count_map.get(locus, 0), cls) for locus, cls in SHOW_LOCI]
# Sort by count descending
bar_data.sort(key=lambda x: -x[1])

CLASS_COLOR = {
    "I":     "#a93226",   # red — classical class I
    "II":    "#1a5c8a",   # blue — classical class II
    "other": "#7f8c8d",   # grey — non-classical / MIC
}
CLASS_LABEL = {"I": "HLA class I", "II": "HLA class II", "other": "Non-classical / MIC"}

FW, FH = 700, 520
PAD_L = 90    # locus labels on left
PAD_R = 30
PAD_T = 80
PAD_B = 70
AW = FW - PAD_L - PAD_R
AH = FH - PAD_T - PAD_B

n = len(bar_data)
bar_h = min(28, (AH - (n - 1) * 6) // n)
gap   = 6
max_count = bar_data[0][1]  # already sorted

bars_svg = ""
for i, (locus, count, cls) in enumerate(bar_data):
    y = PAD_T + i * (bar_h + gap)
    bw = count / max_count * AW
    col = CLASS_COLOR[cls]
    # bar
    bars_svg += (f'<rect x="{PAD_L}" y="{y}" width="{bw:.1f}" height="{bar_h}" '
                 f'fill="{col}" rx="2" opacity="0.85"/>')
    # locus label (left)
    bars_svg += (f'<text x="{PAD_L - 8}" y="{y + bar_h/2 + 4:.1f}" '
                 f'text-anchor="end" font-size="11" fill="#333" font-weight="600">'
                 f'HLA-{locus}</text>')
    # count label (right of bar)
    label_x = PAD_L + bw + 6
    bars_svg += (f'<text x="{label_x:.1f}" y="{y + bar_h/2 + 4:.1f}" '
                 f'font-size="10" fill="{col}" font-weight="600">'
                 f'{count:,}</text>')

# X-axis ticks
xtick_svg = ""
for val in [0, 2000, 4000, 6000, 8000, 10000]:
    if val > max_count: continue
    tx = PAD_L + val / max_count * AW
    xtick_svg += (f'<line x1="{tx:.1f}" y1="{PAD_T}" x2="{tx:.1f}" '
                  f'y2="{PAD_T + AH}" stroke="#eee" stroke-width="1"/>')
    xtick_svg += (f'<text x="{tx:.1f}" y="{PAD_T + AH + 18}" '
                  f'text-anchor="middle" font-size="9" fill="#888">{val:,}</text>')

# Axis line
axis_svg = (f'<line x1="{PAD_L}" y1="{PAD_T}" x2="{PAD_L}" '
            f'y2="{PAD_T+AH}" stroke="#ccc" stroke-width="1"/>')

# Legend
leg_items = [("I", "HLA class I (classical)"), ("II", "HLA class II (classical)"), ("other", "Non-classical / MIC")]
leg_svg = ""
lx = PAD_L + AW * 0.52
ly = PAD_T + 4
for j, (cls, label) in enumerate(leg_items):
    col = CLASS_COLOR[cls]
    leg_svg += (f'<rect x="{lx:.0f}" y="{ly + j*18}" width="12" height="12" '
                f'fill="{col}" rx="2" opacity="0.85"/>')
    leg_svg += (f'<text x="{lx+17:.0f}" y="{ly + j*18 + 10}" '
                f'font-size="10" fill="#444">{label}</text>')

# Highlight annotation for HLA-B (most polymorphic)
hb_y = PAD_T  # first bar (HLA-B)
ann_svg = (
    f'<rect x="{PAD_L + count_map["B"]/max_count*AW - 2:.1f}" y="{hb_y - 2}" '
    f'width="4" height="{bar_h + 4}" rx="1" fill="#c0392b" opacity="0.6"/>'
    f'<text x="{PAD_L + count_map["B"]/max_count*AW + 70:.1f}" y="{hb_y + bar_h/2 - 6:.1f}" '
    f'font-size="9" fill="#c0392b">HLA-B: most polymorphic locus</text>'
    f'<text x="{PAD_L + count_map["B"]/max_count*AW + 70:.1f}" y="{hb_y + bar_h/2 + 6:.1f}" '
    f'font-size="9" fill="#c0392b">any two unrelated people share</text>'
    f'<text x="{PAD_L + count_map["B"]/max_count*AW + 70:.1f}" y="{hb_y + bar_h/2 + 18:.1f}" '
    f'font-size="9" fill="#c0392b">on avg. only ~1 allele in 4</text>'
)

svg1 = f"""<svg viewBox="0 0 {FW} {FH}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{FW//2}" y="22" text-anchor="middle" font-size="13" font-weight="600" fill="#222">
    HLA Allele Diversity per Locus — IPD-IMGT/HLA Release 3.64.0
  </text>
  <text x="{FW//2}" y="40" text-anchor="middle" font-size="10" fill="#666">
    46,005 total known alleles across all HLA loci (real data, July 2026)
  </text>
  <text x="{FW//2}" y="56" text-anchor="middle" font-size="10" fill="#666">
    This is the polymorphism landscape that short-read MHC assemblers must navigate
  </text>
  <text x="{PAD_L + AW/2:.0f}" y="{PAD_T + AH + 38}" text-anchor="middle"
        font-size="10" fill="#555">Number of known alleles</text>
  {axis_svg}{xtick_svg}{bars_svg}{leg_svg}{ann_svg}
</svg>"""

with open(os.path.join(OUT, "hla_polymorphism.svg"), "w") as f:
    f.write(svg1)
print("Wrote hla_polymorphism.svg")


# ── Figure 2: Allele count growth over time ───────────────────────────────────
FW2, FH2 = 700, 400
PAD_L2 = 70
PAD_R2 = 40
PAD_T2 = 70
PAD_B2 = 60
AW2 = FW2 - PAD_L2 - PAD_R2
AH2 = FH2 - PAD_T2 - PAD_B2

years  = [g[0] for g in growth]
totals = [g[1] for g in growth]
min_y, max_y = years[0], years[-1]
min_t, max_t = 0, max(totals) * 1.08

def gx(year):
    return PAD_L2 + (year - min_y) / (max_y - min_y) * AW2

def gy(total):
    return PAD_T2 + AH2 - total / max_t * AH2

# Polyline
pts = " ".join(f"{gx(y):.1f},{gy(t):.1f}" for y, t in growth)
line_svg = (f'<polyline points="{pts}" fill="none" stroke="#1a5c8a" '
            f'stroke-width="2.5" stroke-linejoin="round"/>')

# Filled area under the line
area_pts = (f"{gx(years[0]):.1f},{gy(0):.1f} " + pts +
            f" {gx(years[-1]):.1f},{gy(0):.1f}")
area_svg = (f'<polygon points="{area_pts}" fill="#1a5c8a" opacity="0.12"/>')

# Dots
dots2 = ""
for year, total in growth:
    cx, cy = gx(year), gy(total)
    dots2 += (f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="4" '
              f'fill="#1a5c8a" stroke="white" stroke-width="1.5"/>')

# Axis lines
ax2 = (f'<line x1="{PAD_L2}" y1="{PAD_T2}" x2="{PAD_L2}" '
       f'y2="{PAD_T2+AH2}" stroke="#ccc" stroke-width="1"/>'
       f'<line x1="{PAD_L2}" y1="{PAD_T2+AH2}" x2="{PAD_L2+AW2}" '
       f'y2="{PAD_T2+AH2}" stroke="#ccc" stroke-width="1"/>')

# X ticks
xt2 = ""
for yr in range(2000, 2026, 4):
    tx = gx(yr)
    xt2 += (f'<line x1="{tx:.1f}" y1="{PAD_T2+AH2}" x2="{tx:.1f}" '
            f'y2="{PAD_T2+AH2+5}" stroke="#999" stroke-width="1"/>'
            f'<text x="{tx:.1f}" y="{PAD_T2+AH2+18}" '
            f'text-anchor="middle" font-size="9" fill="#666">{yr}</text>')

# Y ticks
yt2 = ""
for val in [0, 10000, 20000, 30000, 40000]:
    ty = gy(val)
    yt2 += (f'<line x1="{PAD_L2-5}" y1="{ty:.1f}" x2="{PAD_L2}" '
            f'y2="{ty:.1f}" stroke="#999" stroke-width="1"/>'
            f'<text x="{PAD_L2-8}" y="{ty+4:.1f}" text-anchor="end" '
            f'font-size="9" fill="#666">{val//1000}k</text>')

# Annotation: MHConstructor published 2024
mhc_x = gx(2024)
mhc_y = gy(43000)
mhc_ann = (
    f'<line x1="{mhc_x:.1f}" y1="{mhc_y - 8:.1f}" '
    f'x2="{mhc_x:.1f}" y2="{mhc_y - 40:.1f}" '
    f'stroke="#c0392b" stroke-width="1" stroke-dasharray="3,2"/>'
    f'<rect x="{mhc_x - 88:.1f}" y="{mhc_y - 72:.1f}" width="176" height="32" '
    f'rx="3" fill="#fff5f5" stroke="#c0392b" stroke-width="1"/>'
    f'<text x="{mhc_x:.1f}" y="{mhc_y - 55:.1f}" text-anchor="middle" '
    f'font-size="9" fill="#c0392b" font-weight="600">MHConstructor published</text>'
    f'<text x="{mhc_x:.1f}" y="{mhc_y - 43:.1f}" text-anchor="middle" '
    f'font-size="9" fill="#c0392b">Genome Biology, Oct 2024</text>'
)

# Axis labels
ax_labels2 = (
    f'<text x="{PAD_L2 + AW2/2:.0f}" y="{FH2 - 8}" '
    f'text-anchor="middle" font-size="10" fill="#555">Year</text>'
    f'<text transform="rotate(-90,18,{PAD_T2 + AH2/2:.0f})" x="18" '
    f'y="{PAD_T2 + AH2/2:.0f}" text-anchor="middle" font-size="10" fill="#555">'
    f'Total known alleles</text>'
)

svg2 = f"""<svg viewBox="0 0 {FW2} {FH2}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{FW2//2}" y="22" text-anchor="middle" font-size="13" font-weight="600" fill="#222">
    HLA Allele Discovery Is Ongoing — The Reference Is Always Growing
  </text>
  <text x="{FW2//2}" y="40" text-anchor="middle" font-size="10" fill="#666">
    Total alleles in IPD-IMGT/HLA database by year · source: Robinson et al. 2024, HLA
  </text>
  <text x="{FW2//2}" y="56" text-anchor="middle" font-size="10" fill="#666">
    Point-in-time assemblies become outdated; haplotype-informed scaffolding is more robust than direct reference mapping
  </text>
  {ax2}{xt2}{yt2}{area_svg}{line_svg}{dots2}{mhc_ann}{ax_labels2}
</svg>"""

with open(os.path.join(OUT, "allele_growth.svg"), "w") as f:
    f.write(svg2)
print("Wrote allele_growth.svg")
