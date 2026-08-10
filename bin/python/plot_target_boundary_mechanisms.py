#!/usr/bin/env python3
"""Figure 7 -- causal model for outcome changes under interval subsetting.

The figure is drawn in the same visual language as Figure 1 of the manuscript:
the same 180 mm Frontiers double-column width, the same feature palette, flat
fills with a thin dark outline, bold row labels in a left gutter, and a single
bottom legend that defines every colour used.  Where Figure 1 stacks a gene
model over three target tracks, this figure stacks a truth record and a called
record over the target track, so a reader who has seen Figure 1 can read the
rows without re-learning the layout.

All text is at least 8 pt at the printed size, and no element is distinguished
by colour alone: every record is also identified by its row label, and every
outcome is stated in words.

This module is self-contained on purpose.  It is copied into the analysis
Singularity image (see containers/Singularity.python-r-analysis) and run there
as ``plot-target-boundary-mechanisms``, where the sv-article figure package is
not available.  ``sv-article/figures_revision/src/figure7_boundary.py`` imports
``schematic_figure`` from here so that the submission package and the pipeline
render exactly the same artwork.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Patch, Rectangle

# --------------------------------------------------------------------------
# Geometry and type sizes, identical to the Frontiers figure package
# --------------------------------------------------------------------------
MM_PER_IN = 25.4
WIDTH_MM = 180.0          # Frontiers double column
HEIGHT_MM = 104.0
DPI = 400                 # -> 2835 px at 180 mm

FS_TICK = 8               # smallest text on the page
FS_AXIS = 9               # panel titles
FS_LEGEND = 8
FS_PANEL = 11             # (A), (B), (C), (D)

# --------------------------------------------------------------------------
# Palette, continued from Figure 1 so the two schematics read as one pair
# --------------------------------------------------------------------------
C_TARGET_FILL, C_TARGET_EDGE = "#3387d5", "#18459c"   # Figure 1 target blue
C_PAD_FILL = "#c3ddf6"                                 # padded flank
C_TRUTH_FILL, C_TRUTH_EDGE = "#49aac9", "#0f6e86"      # Figure 1 exon teal
C_CALL_FILL, C_CALL_EDGE = "#f4a63c", "#b9741d"        # Figure 1 5' UTR orange
C_PASS = "#2f6b2f"                                     # Figure 1 3' UTR edge
C_FAIL = "#b3261e"
C_INK = "#3f4a52"                                      # Figure 1 intergenic

# --------------------------------------------------------------------------
# Shared layout.  Only the target row changes between panels A, B and D; the
# truth and called records are drawn at identical coordinates throughout, so
# the reader can see that nothing about the caller output changes.
# --------------------------------------------------------------------------
X_MIN, X_MAX = -25.0, 103.0
Y_TITLE = 9.45
Y_TRUTH = 7.10
Y_CALL = 5.30
Y_TARGET = 3.25
Y_CAP1 = 1.35
Y_CAP2 = 0.45

H_REC = 0.60
H_TARGET = 0.54

TRUTH = (30.0, 58.0)
CALL = (44.0, 68.0)
FRAGMENT = (18.0, 40.0)
BROAD = (2.0, 100.0)
PAD = 8.0

ROW_LABELS = ((Y_TRUTH, "Truth"), (Y_CALL, "Call"), (Y_TARGET, "Target"))


def _bar(ax, x0, x1, y, height, fill, edge, linestyle="solid", zorder=4):
    patch = Rectangle(
        (x0, y - height / 2.0), x1 - x0, height,
        facecolor=fill, edgecolor=edge, linewidth=0.7,
        linestyle=linestyle, zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def _panel(fig, col, row, letter, title):
    """One panel axes plus its title, positioned in a 2 x 2 grid."""
    left_mm, gutter_mm, right_mm = 2.0, 5.0, 2.0
    panel_w_mm = (WIDTH_MM - left_mm - right_mm - gutter_mm) / 2.0
    panel_h_mm = 43.0
    top_mm = 1.0
    row_gap_mm = 3.0

    x_mm = left_mm + col * (panel_w_mm + gutter_mm)
    y_mm = HEIGHT_MM - top_mm - (row + 1) * panel_h_mm - row * row_gap_mm

    ax = fig.add_axes([
        x_mm / WIDTH_MM, y_mm / HEIGHT_MM,
        panel_w_mm / WIDTH_MM, panel_h_mm / HEIGHT_MM,
    ])
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(0.0, 10.2)
    ax.axis("off")

    ax.text(X_MIN, Y_TITLE, f"({letter})", ha="left", va="center",
            fontsize=FS_PANEL, fontweight="bold", color="black")
    ax.text(X_MIN + 13.0, Y_TITLE, title, ha="left", va="center",
            fontsize=FS_AXIS, fontweight="bold", color="black")

    for y, name in ROW_LABELS:
        if name == "Call" and letter == "C":
            continue          # panel C is about the truth denominator only
        ax.text(-2.0, y, name, ha="right", va="center",
                fontsize=FS_TICK, fontweight="bold", color="black")
    return ax


def _caption(ax, line1, line2):
    ax.text(X_MIN + 13.0, Y_CAP1, line1, ha="left", va="center",
            fontsize=FS_TICK, color="black")
    ax.text(X_MIN + 13.0, Y_CAP2, line2, ha="left", va="center",
            fontsize=FS_TICK, color="black")


def _edge_guides(ax, edges, top=Y_TRUTH + 0.75):
    for x in edges:
        ax.plot([x, x], [Y_TARGET, top], color=C_TARGET_EDGE, linestyle=":",
                linewidth=0.8, zorder=2)


def _truth_record(ax):
    _bar(ax, *TRUTH, Y_TRUTH, H_REC, C_TRUTH_FILL, C_TRUTH_EDGE)


def _called_record(ax, removed=False):
    if removed:
        _bar(ax, *CALL, Y_CALL, H_REC, "white", C_FAIL, linestyle=(0, (2.4, 1.4)))
        ax.plot([CALL[0], CALL[1]],
                [Y_CALL - H_REC / 2.0, Y_CALL + H_REC / 2.0],
                color=C_FAIL, linewidth=0.9, zorder=5)
    else:
        _bar(ax, *CALL, Y_CALL, H_REC, C_CALL_FILL, C_CALL_EDGE)


def _match_link(ax, broken=False):
    x = (CALL[0] + TRUTH[1]) / 2.0
    colour = "#9aa4ac" if broken else C_INK
    ax.add_patch(FancyArrowPatch(
        (x, Y_TRUTH - H_REC / 2.0 - 0.06), (x, Y_CALL + H_REC / 2.0 + 0.06),
        arrowstyle="<->", mutation_scale=7, linewidth=0.9, color=colour,
        linestyle=(0, (2.0, 1.2)) if broken else "solid", zorder=5,
    ))
    if broken:
        ax.text(x, (Y_TRUTH + Y_CALL) / 2.0, "✕", ha="center", va="center",
                fontsize=FS_TICK, fontweight="bold", color=C_FAIL,
                bbox=dict(boxstyle="circle,pad=0.12", facecolor="white",
                          edgecolor="none"), zorder=6)
        label = "no pair"
    else:
        label = "matched"
    ax.text(x - 2.5, (Y_TRUTH + Y_CALL) / 2.0, label, ha="right", va="center",
            fontsize=FS_TICK, color=colour if broken else C_INK, zorder=6,
            bbox=dict(boxstyle="square,pad=0.12", facecolor="white",
                      edgecolor="none"))


def _outcome(ax, headline, detail, good):
    """Uniform outcome badge in the right margin, clear of every record bar."""
    colour = C_PASS if good else C_FAIL
    x, y = 87.0, (Y_TRUTH + Y_CALL) / 2.0
    ax.add_patch(Rectangle(
        (x - 15.0, y - 1.02), 30.0, 2.04,
        facecolor="white", edgecolor=colour, linewidth=0.7, zorder=5,
    ))
    ax.text(x, y + 0.36, headline, ha="center", va="center",
            fontsize=FS_AXIS, fontweight="bold", color=colour, zorder=6)
    ax.text(x, y - 0.52, detail, ha="center", va="center",
            fontsize=FS_TICK, color=colour, zorder=6)


def _panel_a(fig):
    ax = _panel(fig, 0, 0, "A", "Matching in the broad region")
    _bar(ax, *BROAD, Y_TARGET, H_TARGET, C_TARGET_FILL, C_TARGET_EDGE)
    ax.text((BROAD[0] + BROAD[1]) / 2.0, Y_TARGET, "broad benchmark region",
            ha="center", va="center", fontsize=FS_TICK, fontweight="bold",
            color="white", zorder=6)
    _truth_record(ax)
    _called_record(ax)
    _match_link(ax)
    _outcome(ax, "TP", "true positive", True)
    _caption(ax,
             "Truth and call are displaced, but both stay inside",
             "the region and meet the ≤500 bp / ≥0.7 criteria.")


def _panel_b(fig):
    ax = _panel(fig, 1, 0, "B", "Restriction before matching")
    _bar(ax, *FRAGMENT, Y_TARGET, H_TARGET, C_TARGET_FILL, C_TARGET_EDGE)
    _edge_guides(ax, FRAGMENT)
    _truth_record(ax)
    _called_record(ax, removed=True)
    _match_link(ax, broken=True)
    _outcome(ax, "FN", "false negative", False)
    ax.text((CALL[0] + CALL[1]) / 2.0, Y_CALL - H_REC / 2.0 - 0.52,
            "excluded by the target", ha="center", va="top",
            fontsize=FS_TICK, color=C_FAIL, zorder=6,
            bbox=dict(boxstyle="square,pad=0.12", facecolor="white",
                      edgecolor="none"))
    _caption(ax,
             "The fragmented target keeps the truth and removes",
             "the call, so the same event now scores as a miss.")


def _panel_c(fig):
    ax = _panel(fig, 0, 1, "C", "Membership rule sets the denominator")
    _bar(ax, *FRAGMENT, Y_TARGET, H_TARGET, C_TARGET_FILL, C_TARGET_EDGE)
    _edge_guides(ax, FRAGMENT)
    _truth_record(ax)
    ax.text((TRUTH[0] + TRUTH[1]) / 2.0, Y_TRUTH + H_REC / 2.0 + 0.42,
            "deletion crossing the target edge", ha="center", va="bottom",
            fontsize=FS_TICK, color=C_TRUTH_EDGE, zorder=6,
            bbox=dict(boxstyle="square,pad=0.12", facecolor="white",
                      edgecolor="none"))

    for y, text, colour in ((Y_CALL + 0.35, "containment: 165 truth SVs", C_FAIL),
                            (Y_CALL - 1.15, "overlap: 261 truth SVs", C_PASS)):
        ax.text(70.0, y, text, ha="center", va="center", fontsize=FS_TICK,
                fontweight="bold", color=colour, zorder=6,
                bbox=dict(boxstyle="round,pad=0.30", facecolor="white",
                          edgecolor=colour, linewidth=0.7))
    _caption(ax,
             "Containment drops deletions crossing a target edge;",
             "overlap keeps them, raising the denominator by 58%.")


def _panel_d(fig):
    ax = _panel(fig, 1, 1, "D", "Padding restores the excluded call")
    padded = (FRAGMENT[0] - PAD, FRAGMENT[1] + PAD)
    _bar(ax, padded[0], padded[1], Y_TARGET, H_TARGET, C_PAD_FILL,
         C_TARGET_EDGE, linestyle=(0, (2.4, 1.4)), zorder=3)
    _bar(ax, *FRAGMENT, Y_TARGET, H_TARGET, C_TARGET_FILL, C_TARGET_EDGE)
    _edge_guides(ax, padded)
    _truth_record(ax)
    _called_record(ax)
    _match_link(ax)
    _outcome(ax, "TP", "pair restored", True)


    y_pad = Y_TARGET - H_TARGET / 2.0 - 0.60
    for x0, x1 in ((padded[0], FRAGMENT[0]), (FRAGMENT[1], padded[1])):
        ax.add_patch(FancyArrowPatch(
            (x0, y_pad), (x1, y_pad), arrowstyle="<->", mutation_scale=6,
            linewidth=0.8, color=C_TARGET_EDGE, zorder=5,
        ))
    ax.text(padded[1] + 2.0, y_pad, "+500 bp", ha="left", va="center",
            fontsize=FS_TICK, color=C_TARGET_EDGE)
    _caption(ax,
             "Padding by the 500 bp match distance restores the",
             "pair: 70 of 70 on GRCh37 and 111 of 111 on GRCh38.")


def schematic_figure():
    """Build and return the Figure 7 matplotlib figure."""
    fig = plt.figure(figsize=(WIDTH_MM / MM_PER_IN, HEIGHT_MM / MM_PER_IN))
    fig.patch.set_facecolor("white")

    _panel_a(fig)
    _panel_b(fig)
    _panel_c(fig)
    _panel_d(fig)

    handles = [
        Patch(facecolor=C_TRUTH_FILL, edgecolor=C_TRUTH_EDGE, linewidth=0.7,
              label="Truth SV"),
        Patch(facecolor=C_CALL_FILL, edgecolor=C_CALL_EDGE, linewidth=0.7,
              label="Called SV"),
        Patch(facecolor="white", edgecolor=C_FAIL, linewidth=0.7,
              linestyle=(0, (2.4, 1.4)), label="Call removed by target"),
        Patch(facecolor=C_TARGET_FILL, edgecolor=C_TARGET_EDGE, linewidth=0.7,
              label="Benchmarked target interval"),
        Patch(facecolor=C_PAD_FILL, edgecolor=C_TARGET_EDGE, linewidth=0.7,
              linestyle=(0, (2.4, 1.4)), label="Target padding"),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.6 / HEIGHT_MM),
        ncol=5,
        frameon=False,
        handlelength=1.6,
        handleheight=1.0,
        columnspacing=1.1,
        handletextpad=0.45,
        fontsize=FS_LEGEND,
    )
    return fig


def apply_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "font.size": FS_TICK,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--transitions", type=Path, required=False,
        help="record-level transition table; accepted for pipeline "
             "compatibility and not required to render the schematic",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    apply_style()
    fig = schematic_figure()
    fig.savefig(args.output_dir / "target_boundary_mechanisms.png", dpi=DPI)
    plt.close(fig)


if __name__ == "__main__":
    main()
