#!/usr/bin/env python3
"""
Participant flow diagram (fetal and preterm samples), drawn in the CONSORT / STROBE style:
  * rectangular boxes, thin outlines, black text, "n =" counts
  * a vertical main flow with arrows; exclusion boxes branch off to the side by a horizontal arrow
    from the main line, each stating the REASON and the numbers excluded
  * counts reconcile at every step (checked below)

One layout spec is rendered to SVG (editable in Inkscape/Illustrator; Arial) and to a 600-dpi TIFF/PNG
(Liberation Sans, which has Arial's metrics, so line breaks match).

Edit the numbers / wording in the DATA section and rerun:  python3 make_flowchart.py
Requires: Pillow.  Fonts: Liberation Sans (falls back to DejaVu Sans).
"""
import os
from PIL import Image, ImageDraw, ImageFont

# --------------------------------------------------------------------------------------------
# DATA (edit here)
# --------------------------------------------------------------------------------------------
PRETERM = dict(
    dataset=dict(participants=213, scans=546),
    qa=dict(participants=25, scans=186),
    usable=dict(participants=188, scans=360),
    injury=dict(participants=9, scans=14),
    final=dict(participants=179, scans=346),
)
FETAL = dict(
    dataset=dict(participants=385, scans=535),
    qa=dict(participants=224, scans=343),
    final=dict(participants=161, scans=192),
)
NB = "\u00a0"   # non-breaking space: keeps "< 120 volumes" and "n = 25" from splitting across lines
QA_REASON = f"Insufficient usable data after preprocessing (<{NB}120{NB}volumes, <{NB}4{NB}min)"
INJURY_REASON = f"Kidokoro brain injury score \u2265{NB}9"

# Palette (matches the figures)
RED, RED_TINT = "#cc2936", "#f2c9cd"
BLUE, BLUE_TINT = "#2a78d6", "#cadff5"
INK = "#000000"

# Geometry, in points (1 pt = 1/72 in). Final size: 7.2 in wide.
W = 518.4
MARGIN = 10
EXCL_W, ARROW_GAP, MAIN_W = 118.0, 16.0, 104.0
CENTER_GAP = W - 2 * MARGIN - 2 * (EXCL_W + ARROW_GAP + MAIN_W)
PAD = 6.0                       # inner padding of boxes
FS, FS_TITLE = 7.0, 8.0         # body / title font sizes (pt)
LEAD = 1.28                     # line spacing (x font size)
V_GAP = 34.0                    # vertical distance between boxes (room for the arrow)
LW_BOX, LW_ARROW = 1.0, 1.0     # line widths (pt)
HEAD_L, HEAD_W = 5.5, 4.6       # arrowhead length / width (pt)

FONT_DIRS = ["/usr/share/fonts/truetype/liberation", "/usr/share/fonts/truetype/dejavu"]
def _font_path(bold):
    for d, names in [(FONT_DIRS[0], ("LiberationSans-Bold.ttf" if bold else "LiberationSans-Regular.ttf",)),
                     (FONT_DIRS[1], ("DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf",))]:
        for n in names:
            p = os.path.join(d, n)
            if os.path.exists(p):
                return p
    raise FileNotFoundError("No Liberation Sans / DejaVu Sans font found")

_MEASURE_DPI = 600
def text_width(s, size_pt, bold=False):
    f = ImageFont.truetype(_font_path(bold), round(size_pt * _MEASURE_DPI / 72 * 4))   # x4 for precision
    return f.getlength(s) / 4 * 72 / _MEASURE_DPI

def wrap(s, size_pt, width_pt, bold=False):
    words, lines, cur = s.split(" "), [], ""   # split on plain spaces only, so non-breaking spaces keep words together
    for w in words:
        t = (cur + " " + w).strip()
        if text_width(t, size_pt, bold) <= width_pt or not cur:
            cur = t
        else:
            lines.append(cur); cur = w
    if cur:
        lines.append(cur)
    return lines

# --------------------------------------------------------------------------------------------
# Checks: counts reconcile at every step
# --------------------------------------------------------------------------------------------
for k in ("participants", "scans"):
    assert PRETERM["dataset"][k] - PRETERM["qa"][k] == PRETERM["usable"][k], k
    assert PRETERM["usable"][k] - PRETERM["injury"][k] == PRETERM["final"][k], k
    assert FETAL["dataset"][k] - FETAL["qa"][k] == FETAL["final"][k], k

# --------------------------------------------------------------------------------------------
# Layout spec: list of primitives (all in pt, y down)
# --------------------------------------------------------------------------------------------
prims = []   # ("rect", x, y, w, h, fill, stroke, lw) | ("line", x1,y1,x2,y2,color,lw) | ("arrow", x1,y1,x2,y2,color,lw) | ("text", x, y_baseline, s, size, bold)

def box_lines(title, body, width, size_title=FS_TITLE, size=FS):
    """Return [(text, size, bold)] with wrapping inside the box's inner width."""
    inner = width - 2 * PAD
    out = [(l, size_title, True) for l in wrap(title, size_title, inner, True)]
    for b in body:
        out += [(l, size, False) for l in wrap(b, size, inner, False)]
    return out

def box_height(lines):
    return 2 * PAD + sum(sz * LEAD for _, sz, _ in lines) - (LEAD - 1) * lines[-1][1] * 0.0

def add_box(x, y, w, lines, fill, stroke, lw=LW_BOX):
    h = box_height(lines)
    prims.append(("rect", x, y, w, h, fill, stroke, lw))
    cy = y + PAD
    for i, (t, sz, bold) in enumerate(lines):
        cy += sz * LEAD
        base = cy - sz * (LEAD - 1) - sz * 0.22   # baseline inside the line box
        prims.append(("text", x + PAD, base, t, sz, bold))
    return h

def n_lines(d):
    return [f"Participants: n{NB}={NB}{d['participants']}", f"Scans: n{NB}={NB}{d['scans']}"]

x_pt_excl = MARGIN
x_pt_main = x_pt_excl + EXCL_W + ARROW_GAP
x_f_main = x_pt_main + MAIN_W + CENTER_GAP
x_f_excl = x_f_main + MAIN_W + ARROW_GAP

top = MARGIN
# ---- main boxes (preterm) ----
L1 = box_lines("Preterm dataset", n_lines(PRETERM["dataset"]), MAIN_W)
L2 = box_lines("Usable preterm data", n_lines(PRETERM["usable"]), MAIN_W)
L3 = box_lines("Final preterm sample", n_lines(PRETERM["final"]), MAIN_W)
F1 = box_lines("Fetal dataset", n_lines(FETAL["dataset"]), MAIN_W)
F3 = box_lines("Final fetal sample", n_lines(FETAL["final"]), MAIN_W)

h1, h2, h3 = box_height(L1), box_height(L2), box_height(L3)
y1 = top
y2 = y1 + h1 + V_GAP
y3 = y2 + h2 + V_GAP
hf1 = box_height(F1)

add_box(x_pt_main, y1, MAIN_W, L1, BLUE_TINT, BLUE)
add_box(x_pt_main, y2, MAIN_W, L2, BLUE_TINT, BLUE)
add_box(x_pt_main, y3, MAIN_W, L3, BLUE_TINT, BLUE)
add_box(x_f_main, y1, MAIN_W, F1, RED_TINT, RED)
add_box(x_f_main, y3, MAIN_W, F3, RED_TINT, RED)

# ---- main arrows ----
cx_pt, cx_f = x_pt_main + MAIN_W / 2, x_f_main + MAIN_W / 2
prims.append(("arrow", cx_pt, y1 + h1, cx_pt, y2, BLUE, LW_ARROW))
prims.append(("arrow", cx_pt, y2 + h2, cx_pt, y3, BLUE, LW_ARROW))
prims.append(("arrow", cx_f, y1 + hf1, cx_f, y3, RED, LW_ARROW))

# ---- exclusion boxes: branch from the main arrow, into a box that states the reason and the counts ----
def excl(title, reason, d, who):
    return box_lines(title, [reason, f"Scans: n{NB}={NB}{d['scans']}", f"{who}: n{NB}={NB}{d['participants']}"], EXCL_W)

E1 = excl("Excluded: data quality", QA_REASON, PRETERM["qa"], "Participants excluded")
E2 = excl("Excluded: brain injury", INJURY_REASON, PRETERM["injury"], "Participants")
E3 = excl("Excluded: data quality", QA_REASON, FETAL["qa"], "Participants excluded")

j1 = y1 + h1 + V_GAP / 2              # junction heights (middle of each arrow)
j2 = y2 + h2 + V_GAP / 2
for (lines, yj, xbox, xmain_c, color, side) in [
        (E1, j1, x_pt_excl, cx_pt, BLUE, "left"),
        (E2, j2, x_pt_excl, cx_pt, BLUE, "left"),
        (E3, j1, x_f_excl, cx_f, RED, "right")]:
    hh = box_height(lines)
    add_box(xbox, yj - hh / 2, EXCL_W, lines, "#ffffff", color)
    xedge = xbox + EXCL_W if side == "left" else xbox
    prims.append(("arrow", xmain_c, yj, xedge, yj, color, LW_ARROW))

H = y3 + h3 + MARGIN

# --------------------------------------------------------------------------------------------
# Renderers
# --------------------------------------------------------------------------------------------
def arrow_head(x1, y1, x2, y2):
    import math
    dx, dy = x2 - x1, y2 - y1
    L = math.hypot(dx, dy); ux, uy = dx / L, dy / L
    bx, by = x2 - ux * HEAD_L, y2 - uy * HEAD_L
    px, py = -uy, ux
    return [(x2, y2), (bx + px * HEAD_W / 2, by + py * HEAD_W / 2), (bx - px * HEAD_W / 2, by - py * HEAD_W / 2)], (bx, by)

def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

def render_svg(path):
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}pt" height="{H:.1f}pt" viewBox="0 0 {W} {H:.2f}">',
         f'<rect width="{W}" height="{H:.2f}" fill="#ffffff"/>',
         '<g font-family="Arial, Helvetica, \'Liberation Sans\', sans-serif" fill="#000000">']
    for p in prims:
        if p[0] == "rect":
            _, x, y, w, h, fill, stroke, lw = p
            o.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" fill="{fill}" stroke="{stroke}" stroke-width="{lw}"/>')
        elif p[0] == "arrow":
            _, x1, y1, x2, y2, color, lw = p
            tri, (bx, by) = arrow_head(x1, y1, x2, y2)
            o.append(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{bx:.2f}" y2="{by:.2f}" stroke="{color}" stroke-width="{lw}"/>')
            o.append('<polygon points="' + " ".join(f"{a:.2f},{b:.2f}" for a, b in tri) + f'" fill="{color}"/>')
        elif p[0] == "text":
            _, x, y, s, sz, bold = p
            o.append(f'<text x="{x:.2f}" y="{y:.2f}" font-size="{sz}" font-weight="{"bold" if bold else "normal"}">{esc(s)}</text>')
    o.append("</g></svg>")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(o))

def render_raster(path_tif, path_png, dpi=600, ss=2):
    k = dpi * ss / 72.0
    img = Image.new("RGB", (round(W * k), round(H * k)), "#ffffff")
    d = ImageDraw.Draw(img)
    fonts = {}
    def font(sz, bold):
        key = (sz, bold)
        if key not in fonts:
            fonts[key] = ImageFont.truetype(_font_path(bold), round(sz * k))
        return fonts[key]
    for p in prims:
        if p[0] == "rect":
            _, x, y, w, h, fill, stroke, lw = p
            half = lw * k / 2
            d.rectangle([x * k - half, y * k - half, (x + w) * k + half, (y + h) * k + half], fill=stroke)
            d.rectangle([x * k + half, y * k + half, (x + w) * k - half, (y + h) * k - half], fill=fill)
        elif p[0] == "arrow":
            _, x1, y1, x2, y2, color, lw = p
            tri, (bx, by) = arrow_head(x1, y1, x2, y2)
            d.line([x1 * k, y1 * k, bx * k, by * k], fill=color, width=max(1, round(lw * k)))
            d.polygon([(a * k, b * k) for a, b in tri], fill=color)
        elif p[0] == "text":
            _, x, y, s, sz, bold = p
            d.text((x * k, y * k), s, font=font(sz, bold), fill=INK, anchor="ls")
    img = img.resize((round(W * dpi / 72), round(H * dpi / 72)), Image.LANCZOS)
    img.save(path_tif, compression="tiff_lzw", dpi=(dpi, dpi))
    img.save(path_png, dpi=(dpi, dpi))

if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    render_svg(os.path.join(here, "Flowchart_redrawn.svg"))
    render_raster(os.path.join(here, "Flowchart_redrawn.tif"), os.path.join(here, "Flowchart_redrawn.png"))
    print(f"Wrote Flowchart_redrawn.svg/.tif/.png  ({W/72:.2f} x {H/72:.2f} in)")
