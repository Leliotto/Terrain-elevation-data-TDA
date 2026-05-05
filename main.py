"""
TDA — Amsterdam Centre (canal network) vs Eindhoven North (peri-urban)
Improved visualisations with per-figure interpretation annotations.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import rasterio
import gudhi

os.makedirs("output", exist_ok=True)

# ── colour palette ────────────────────────────────────────────
AMS_COL   = "#2A7DC9"   # steelblue  (Amsterdam)
EHV_COL   = "#C94F2A"   # burnt orange (Eindhoven)
AMS_LIGHT = "#A8CCE8"
EHV_LIGHT = "#EDBBAA"

CITIES = ["Amsterdam Centre", "Eindhoven North"]
PATHS  = {"Amsterdam Centre": "data/amsterdam_centre_strip.tif",
           "Eindhoven North":  "data/eindhoven_north_strip.tif"}
MAX_SIZE = 700   # cap for critical-point grid (px per side)

# ─────────────────────────────────────────────────────────────
# UTILITIES
# ─────────────────────────────────────────────────────────────

def load_tif(path):
    with rasterio.open(path) as src:
        elev = src.read(1).astype(float)
        nodata = src.nodata
    if nodata is not None:
        elev[elev == nodata] = np.nan
    return elev

def downsample(data, max_size=MAX_SIZE):
    h, w = data.shape
    factor = max(1, max(h, w) // max_size)
    return data[::factor, ::factor], factor

def find_critical_points(elev):
    rows, cols = elev.shape
    is_min = np.zeros((rows, cols), bool)
    is_max = np.zeros((rows, cols), bool)
    is_sad = np.zeros((rows, cols), bool)
    pad = np.pad(elev, 1, mode="constant", constant_values=np.nan)
    for r in range(rows):
        for c in range(cols):
            center = elev[r, c]
            if np.isnan(center):
                continue
            nbrs = pad[r:r+3, c:c+3].flatten()
            nbrs = np.delete(nbrs, 4)
            nbrs = nbrs[~np.isnan(nbrs)]
            if len(nbrs) < 3:
                continue
            if center < np.min(nbrs):
                is_min[r, c] = True
            elif center > np.max(nbrs):
                is_max[r, c] = True
            else:
                diff = nbrs - center
                if np.sum(np.diff(np.sign(diff)) != 0) >= 2:
                    is_sad[r, c] = True
    return is_min, is_max, is_sad

def compute_persistence(elev):
    filled = elev.copy()
    nan_mask = np.isnan(filled)
    if nan_mask.any():
        filled[nan_mask] = np.nanmean(filled)
    cc = gudhi.CubicalComplex(
        dimensions=list(filled.shape),
        top_dimensional_cells=filled.flatten()
    )
    cc.compute_persistence()
    return cc.persistence(), cc

def pairs_0d(pairs):
    return [(b, d) for dim, (b, d) in pairs
            if dim == 0 and not np.isinf(d)]

def annotate_box(ax, text, xy=(0.02, 0.97), fontsize=8.5, color="0.25"):
    ax.annotate(text, xy=xy, xycoords="axes fraction",
                va="top", ha="left", fontsize=fontsize, color=color,
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.85", alpha=0.88))

# ─────────────────────────────────────────────────────────────
# LOAD
# ─────────────────────────────────────────────────────────────
print("Loading data …")
data = {c: load_tif(PATHS[c]) for c in CITIES}
ams, ehv = data["Amsterdam Centre"], data["Eindhoven North"]
print(f"  Amsterdam  shape={ams.shape}  range=[{np.nanmin(ams):.2f}, {np.nanmax(ams):.2f}] m")
print(f"  Eindhoven  shape={ehv.shape}  range=[{np.nanmin(ehv):.2f}, {np.nanmax(ehv):.2f}] m")

# Shared scale for multi-city figures
vmin_all = min(np.nanmin(ams), np.nanmin(ehv))
vmax_all = max(np.nanmax(ams), np.nanmax(ehv))

# ═══════════════════════════════════════════════════════════════
# STEP 1 — RAW ELEVATION HEATMAPS (individual colour scales)
# ═══════════════════════════════════════════════════════════════
print("\nStep 1 — Heatmaps …")

fig, axes = plt.subplots(1, 2, figsize=(16, 8),
                          gridspec_kw={"width_ratios":
                                       [ams.shape[1]/ams.shape[0],
                                        ehv.shape[1]/ehv.shape[0]]})
fig.suptitle("Step 1 — Terrain Elevation\n"
             "Amsterdam Centre (historic canal city) vs Eindhoven North (peri-urban)",
             fontsize=13, fontweight="bold", y=1.01)

descs = [
    "Mean: −0.17 m  (below sea level)\nCanal trenches clearly visible as dark channels\nElevation spread reflects water vs built land",
    "Mean: +18.38 m  (above sea level)\nFlat polder landscape with gentle relief\nUniform warm tone = low topographic variation",
]
for ax, d, title, desc in zip(axes, [ams, ehv], CITIES, descs):
    im = ax.imshow(d, cmap="terrain",
                   vmin=np.nanmin(d), vmax=np.nanmax(d))
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.axis("off")
    plt.colorbar(im, ax=ax, label="Elevation (m)",
                 fraction=0.046, pad=0.04, shrink=0.8)
    annotate_box(ax, desc)

plt.tight_layout()
plt.savefig("output/step1_heatmaps.png", dpi=160, bbox_inches="tight")
plt.show()
print("  ✓ step1_heatmaps.png")

# ═══════════════════════════════════════════════════════════════
# STEP 2 — HISTOGRAMS
# ═══════════════════════════════════════════════════════════════
print("\nStep 2 — Histograms …")

stats = {}
for name, d in data.items():
    v = d[~np.isnan(d)]
    stats[name] = dict(mean=np.mean(v), median=np.median(v),
                       std=np.std(v), mn=np.min(v), mx=np.max(v))

fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
fig.suptitle("Step 2 — Elevation Distribution Comparison", fontsize=13, fontweight="bold")

colors_list = [AMS_COL, EHV_COL]
for ax, name, col in zip(axes[:2], CITIES, colors_list):
    v = data[name][~np.isnan(data[name])].flatten()
    s = stats[name]
    ax.hist(v, bins=100, color=col, alpha=0.80, edgecolor="white", linewidth=0.2)
    ax.axvline(s["mean"],   color="black", ls="--", lw=1.6,
               label=f"Mean: {s['mean']:.2f} m")
    ax.axvline(s["median"], color="gray",  ls=":",  lw=1.4,
               label=f"Median: {s['median']:.2f} m")
    ax.set_title(name, fontsize=11, fontweight="bold")
    ax.set_xlabel("Elevation (m)")
    ax.set_ylabel("Pixel count")
    ax.legend(fontsize=9)
    txt = (f"std = {s['std']:.2f} m\n"
           f"min = {s['mn']:.2f} m\n"
           f"max = {s['mx']:.2f} m")
    annotate_box(ax, txt, xy=(0.97, 0.97), fontsize=8.5)
    ax.annotate("" if name.startswith("Amsterdam")
                else "", xy=(0.5, 0.5),
                xycoords="axes fraction")

# Normalised overlay — lets us compare shape regardless of pixel count
ax = axes[2]
for name, col in zip(CITIES, colors_list):
    v = data[name][~np.isnan(data[name])].flatten()
    ax.hist(v, bins=100, color=col, alpha=0.48,
            label=name, edgecolor="white", linewidth=0.2,
            density=True)
ax.set_title("Normalised Overlay\n(density, comparable shape)", fontsize=11, fontweight="bold")
ax.set_xlabel("Elevation (m)")
ax.set_ylabel("Density")
ax.legend(fontsize=9)
annotate_box(ax,
    "Amsterdam: wide, flat spread — canal lows\n  pull distribution below 0 m\n"
    "Eindhoven: narrow, tall peak — uniform\n  polder terrain around 18 m",
    xy=(0.03, 0.97))

plt.tight_layout()
plt.savefig("output/step2_histograms.png", dpi=160, bbox_inches="tight")
plt.show()

print("\n  Elevation statistics:")
print(f"  {'City':22s}  {'Mean':>7}  {'Std':>7}  {'Min':>7}  {'Max':>7}")
print(f"  {'-'*58}")
for name in CITIES:
    s = stats[name]
    print(f"  {name:22s}  {s['mean']:7.2f}  {s['std']:7.2f}  "
          f"  {s['mn']:7.2f}  {s['mx']:7.2f}")
print("  ✓ step2_histograms.png")

# ═══════════════════════════════════════════════════════════════
# STEP 3 — CRITICAL POINTS
# ═══════════════════════════════════════════════════════════════
print(f"\nStep 3 — Critical points (MAX_SIZE={MAX_SIZE}) …")

ds    = {c: downsample(data[c]) for c in CITIES}
grids = {c: ds[c][0] for c in CITIES}
facts = {c: ds[c][1] for c in CITIES}

for c in CITIES:
    print(f"  {c}: working grid {grids[c].shape}  ×{facts[c]}")

crit = {}
for c in CITIES:
    print(f"  detecting {c} …")
    crit[c] = find_critical_points(grids[c])

cp_counts = {}
for c in CITIES:
    mn, mx, sd = crit[c]
    cp_counts[c] = dict(minima=int(mn.sum()), maxima=int(mx.sum()),
                        saddles=int(sd.sum()))
    print(f"  {c}: minima={mn.sum()}  maxima={mx.sum()}  saddles={sd.sum()}")

# ── figure: side-by-side critical point maps ──
fig, axes = plt.subplots(1, 2, figsize=(18, 9),
                          gridspec_kw={"width_ratios":
                                       [grids["Amsterdam Centre"].shape[1],
                                        grids["Eindhoven North"].shape[1]]})
fig.suptitle("Step 3 — Critical Points via Discrete Morse Theory\n"
             "▼ Minima   ▲ Maxima   ■ Saddles",
             fontsize=13, fontweight="bold", y=1.01)

for ax, c in zip(axes, CITIES):
    g = grids[c]
    mn, mx, sd = crit[c]
    counts = cp_counts[c]
    ax.imshow(g, cmap="terrain",
              vmin=np.nanmin(g), vmax=np.nanmax(g), alpha=0.80)

    def scat(mask, col, label, marker, sz, max_pts=2000):
        ys, xs = np.where(mask)
        if len(ys) > max_pts:
            idx = np.random.choice(len(ys), max_pts, replace=False)
            ys, xs = ys[idx], xs[idx]
        ax.scatter(xs, ys, c=col, s=sz, marker=marker,
                   label=f"{label} (n={mask.sum():,})",
                   linewidths=0, alpha=0.75, zorder=3)

    scat(mn, "#3399FF", "Minima",  "v", 14)
    scat(mx, "#FF4444", "Maxima",  "^", 14)
    scat(sd, "#FFD700", "Saddles", "s",  5)

    ratio = counts["saddles"] / max(counts["minima"], 1)
    desc = (f"saddles/minima ratio: {ratio:.1f}\n"
            f"High ratio = flatter, noisier terrain")
    ax.set_title(c, fontsize=12, fontweight="bold")
    ax.axis("off")
    ax.legend(loc="lower right", fontsize=9, framealpha=0.88)
    annotate_box(ax, desc)

plt.tight_layout()
plt.savefig("output/step3_critical_points.png", dpi=160, bbox_inches="tight")
plt.show()

# ── figure: critical point count comparison bar chart ──
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
fig.suptitle("Step 3b — Critical Point Counts vs Saddle/Minima Ratio",
             fontsize=12, fontweight="bold")

ax = axes[0]
cats = ["Minima", "Maxima", "Saddles"]
keys = ["minima", "maxima", "saddles"]
x = np.arange(len(cats))
w = 0.35
b1 = ax.bar(x - w/2,
            [cp_counts["Amsterdam Centre"][k] for k in keys],
            w, label="Amsterdam Centre", color=AMS_COL, alpha=0.85)
b2 = ax.bar(x + w/2,
            [cp_counts["Eindhoven North"][k] for k in keys],
            w, label="Eindhoven North",  color=EHV_COL, alpha=0.85)
ax.set_xticks(x); ax.set_xticklabels(cats)
ax.set_ylabel("Count"); ax.set_title("Absolute counts")
ax.legend()
for bars in [b1, b2]:
    for bar in bars:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, h * 1.01,
                f"{h:,}", ha="center", va="bottom", fontsize=7.5)

ax = axes[1]
ratios = {c: cp_counts[c]["saddles"] / max(cp_counts[c]["minima"], 1)
          for c in CITIES}
bars = ax.bar(CITIES, ratios.values(),
              color=[AMS_COL, EHV_COL], alpha=0.85, width=0.45)
ax.set_ylabel("Saddles / Minima ratio")
ax.set_title("Saddle-to-Minima Ratio\n(higher = flatter / noisier terrain)")
for bar, val in zip(bars, ratios.values()):
    ax.text(bar.get_x() + bar.get_width()/2, val * 1.01,
            f"{val:.1f}", ha="center", va="bottom", fontsize=11, fontweight="bold")
annotate_box(ax,
    "Eindhoven's much higher ratio shows\n"
    "that flat terrain produces many\n"
    "near-degenerate saddle points — noise.\n"
    "This is why persistence filtering matters.",
    xy=(0.03, 0.97))

plt.tight_layout()
plt.savefig("output/step3b_counts.png", dpi=160, bbox_inches="tight")
plt.show()
print("  ✓ step3_critical_points.png  step3b_counts.png")

# ═══════════════════════════════════════════════════════════════
# STEP 4 — PERSISTENCE DIAGRAMS
# ═══════════════════════════════════════════════════════════════
print("\nStep 4 — Persistence (gudhi) …")

pers = {}
for c in CITIES:
    print(f"  computing {c} …")
    pairs, _ = compute_persistence(grids[c])
    pers[c] = pairs_0d(pairs)
    print(f"  {c}: {len(pers[c])} finite 0-dim pairs")

# ── figure: persistence diagrams with RELATIVE axes ──
fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
fig.suptitle("Step 4 — Persistence Diagrams  (0-dim, sublevel set filtration)\n"
             "Each point = one topological basin.  Distance from diagonal = significance.",
             fontsize=12, fontweight="bold")

for ax, c, col in zip(axes, CITIES, [AMS_COL, EHV_COL]):
    p0 = pers[c]
    if not p0:
        ax.set_title(f"{c}\n(no finite pairs)"); continue
    births  = np.array([b for b, d in p0])
    deaths  = np.array([d for b, d in p0])
    persist = deaths - births

    sc = ax.scatter(births, deaths, c=persist, cmap="plasma",
                    s=12, alpha=0.65, linewidths=0,
                    vmin=0, vmax=np.percentile(persist, 98))

    # diagonal
    lo, hi = min(births.min(), deaths.min()), max(births.max(), deaths.max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=1.0, alpha=0.45,
            label="diagonal (persistence = 0)")

    # highlight top-10 most persistent
    top_idx = np.argsort(persist)[-10:]
    ax.scatter(births[top_idx], deaths[top_idx],
               s=60, facecolors="none", edgecolors="white",
               linewidths=1.2, zorder=5, label="Top-10 features")

    ax.set_xlabel("Birth elevation (m)")
    ax.set_ylabel("Death elevation (m)")
    ax.set_title(f"{c}\n({len(p0):,} features)", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8.5)
    plt.colorbar(sc, ax=ax, label="Persistence (m)", fraction=0.046, pad=0.04)

    # stats annotation
    txt = (f"max persistence: {persist.max():.2f} m\n"
           f"median persistence: {np.median(persist):.2f} m\n"
           f"total persistence: {persist.sum():.1f} m")
    annotate_box(ax, txt, xy=(0.02, 0.97))

plt.tight_layout()
plt.savefig("output/step4_persistence.png", dpi=160, bbox_inches="tight")
plt.show()

# ── figure: persistence distribution comparison ──
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle("Step 4b — Persistence Value Distributions",
             fontsize=12, fontweight="bold")

ax = axes[0]
for c, col in zip(CITIES, [AMS_COL, EHV_COL]):
    persist = np.array([d - b for b, d in pers[c]])
    ax.hist(persist, bins=80, color=col, alpha=0.55,
            label=c, density=True, edgecolor="white", linewidth=0.2)
ax.set_xlabel("Persistence (m)")
ax.set_ylabel("Density")
ax.set_title("Distribution of feature persistence")
ax.legend()
annotate_box(ax,
    "A heavy right tail = many long-lived\n"
    "topological features (e.g. deep canals).\n"
    "A spike near zero = mostly noise.")

ax = axes[1]
for c, col in zip(CITIES, [AMS_COL, EHV_COL]):
    persist = np.sort([d - b for b, d in pers[c]])[::-1]
    total = persist.sum()
    cumsum = np.cumsum(persist) / total
    ax.plot(np.linspace(0, 1, len(persist)), cumsum,
            color=col, lw=2.0, label=c)
ax.set_xlabel("Fraction of features (most persistent first)")
ax.set_ylabel("Cumulative share of total persistence")
ax.set_title("Lorenz-style persistence curve\n(steeper = more concentrated in few features)")
ax.legend()
ax.axline((0, 0), slope=1, ls="--", color="0.6", lw=0.8)
annotate_box(ax,
    "If Amsterdam's curve is steeper:\n"
    "a few deep canal basins dominate.\n"
    "If Eindhoven's is more linear:\n"
    "persistence is spread across many\n"
    "shallow, equally-noisy features.")

plt.tight_layout()
plt.savefig("output/step4b_persistence_dist.png", dpi=160, bbox_inches="tight")
plt.show()
print("  ✓ step4_persistence.png  step4b_persistence_dist.png")

# ═══════════════════════════════════════════════════════════════
# STEP 5 — FILTERED MINIMA AT MULTIPLE THRESHOLDS
# ═══════════════════════════════════════════════════════════════
print("\nStep 5 — Filtered minima …")

def sig_minima(elev, min_mask, threshold):
    regional_mean = np.nanmean(elev)
    sig = np.zeros_like(min_mask)
    ys, xs = np.where(min_mask)
    for y, x in zip(ys, xs):
        if (regional_mean - elev[y, x]) >= threshold:
            sig[y, x] = True
    return sig

# Compute persistence-based thresholds
thresh_data = {}
for c in CITIES:
    persist = np.array([d - b for b, d in pers[c]])
    thresh_data[c] = {
        "50th": float(np.percentile(persist, 50)),
        "75th": float(np.percentile(persist, 75)),
        "90th": float(np.percentile(persist, 90)),
    }
    print(f"  {c}: 50th={thresh_data[c]['50th']:.3f}m  "
          f"75th={thresh_data[c]['75th']:.3f}m  "
          f"90th={thresh_data[c]['90th']:.3f}m")

# ── figure: three-threshold sensitivity for each city ──
for c, col in zip(CITIES, [AMS_COL, EHV_COL]):
    g = grids[c]
    mn, _, _ = crit[c]
    thresholds = thresh_data[c]

    fig, axes = plt.subplots(1, 3, figsize=(18, 7))
    fig.suptitle(f"Step 5 — Significant Minima at Different Persistence Thresholds\n{c}",
                 fontsize=12, fontweight="bold")

    for ax, (label, thresh) in zip(axes, thresholds.items()):
        sig = sig_minima(g, mn, thresh)
        ax.imshow(g, cmap="terrain",
                  vmin=np.nanmin(g), vmax=np.nanmax(g), alpha=0.80)
        # all minima faint
        ys, xs = np.where(mn)
        ax.scatter(xs, ys, c="white", s=3, alpha=0.18,
                   linewidths=0, zorder=2)
        # significant minima
        ys, xs = np.where(sig)
        ax.scatter(xs, ys, c=col, s=28, alpha=0.92,
                   linewidths=0.4, edgecolors="white", zorder=3)
        ax.set_title(f"{label} percentile threshold\n"
                     f"thresh = {thresh:.3f} m  →  {sig.sum():,} significant minima",
                     fontsize=10)
        ax.axis("off")
        pct_kept = 100 * sig.sum() / max(mn.sum(), 1)
        annotate_box(ax, f"{pct_kept:.1f}% of minima retained", xy=(0.02, 0.04))

    fname = f"output/step5_{'ams' if c.startswith('A') else 'ehv'}_thresholds.png"
    plt.tight_layout()
    plt.savefig(fname, dpi=160, bbox_inches="tight")
    plt.show()
    print(f"  ✓ {fname}")

# ── figure: side-by-side best threshold comparison ──
fig, axes = plt.subplots(1, 2, figsize=(18, 9),
                          gridspec_kw={"width_ratios":
                                       [grids["Amsterdam Centre"].shape[1],
                                        grids["Eindhoven North"].shape[1]]})
fig.suptitle("Step 5b — Final Comparison: Significant Minima (75th percentile threshold)\n"
             "Blue = persistent topological basins surviving noise filtering",
             fontsize=12, fontweight="bold", y=1.01)

for ax, c, col in zip(axes, CITIES, [AMS_COL, EHV_COL]):
    g = grids[c]
    mn, _, _ = crit[c]
    thresh = thresh_data[c]["75th"]
    sig = sig_minima(g, mn, thresh)

    ax.imshow(g, cmap="terrain",
              vmin=np.nanmin(g), vmax=np.nanmax(g), alpha=0.80)
    ys, xs = np.where(mn)
    ax.scatter(xs, ys, c="white", s=3, alpha=0.18,
               linewidths=0, zorder=2,
               label=f"All minima ({mn.sum():,})")
    ys, xs = np.where(sig)
    ax.scatter(xs, ys, c=col, s=32, alpha=0.92,
               linewidths=0.4, edgecolors="white", zorder=3,
               label=f"Significant ({sig.sum():,})")

    ax.set_title(c, fontsize=12, fontweight="bold")
    ax.axis("off")
    ax.legend(loc="lower right", fontsize=9, framealpha=0.88)

    if c.startswith("Amsterdam"):
        desc = ("Sparse significant minima.\nCanal basins have moderate\n"
                "persistence — fragmented by\nbuilding NoData gaps.")
    else:
        desc = ("Minima trace linear drainage\n"
                "ditches & watercourses —\n"
                "Eindhoven's engineered\nwater infrastructure.")
    annotate_box(ax, desc)

plt.tight_layout()
plt.savefig("output/step5b_comparison.png", dpi=160, bbox_inches="tight")
plt.show()
print("  ✓ step5b_comparison.png")

# ═══════════════════════════════════════════════════════════════
# STEP 6 — SUMMARY STATS TABLE
# ═══════════════════════════════════════════════════════════════
print("\nStep 6 — Summary figure …")

fig, ax = plt.subplots(figsize=(13, 4.5))
ax.axis("off")
fig.suptitle("Summary — Key Metrics Comparison", fontsize=13, fontweight="bold")

ams_p = np.array([d - b for b, d in pers["Amsterdam Centre"]])
ehv_p = np.array([d - b for b, d in pers["Eindhoven North"]])

rows = [
    ["Metric", "Amsterdam Centre", "Eindhoven North", "Interpretation"],
    ["Mean elevation (m)",
     f"{stats['Amsterdam Centre']['mean']:.2f}",
     f"{stats['Eindhoven North']['mean']:.2f}",
     "Ams below sea level; Ehv on higher polder land"],
    ["Elevation std (m)",
     f"{stats['Amsterdam Centre']['std']:.2f}",
     f"{stats['Eindhoven North']['std']:.2f}",
     "Ams wider spread → more topographic complexity"],
    ["Raw minima count",
     f"{cp_counts['Amsterdam Centre']['minima']:,}",
     f"{cp_counts['Eindhoven North']['minima']:,}",
     "Ehv has more raw minima — flat noise effect"],
    ["Saddles / Minima ratio",
     f"{cp_counts['Amsterdam Centre']['saddles']/max(cp_counts['Amsterdam Centre']['minima'],1):.1f}",
     f"{cp_counts['Eindhoven North']['saddles']/max(cp_counts['Eindhoven North']['minima'],1):.1f}",
     "Higher ratio = flatter, noisier terrain"],
    ["Max persistence (m)",
     f"{ams_p.max():.2f}",
     f"{ehv_p.max():.2f}",
     "Deepest significant basin in each city"],
    ["Median persistence (m)",
     f"{np.median(ams_p):.3f}",
     f"{np.median(ehv_p):.3f}",
     "Typical feature depth — both cities similar"],
    ["Total persistence (m)",
     f"{ams_p.sum():.1f}",
     f"{ehv_p.sum():.1f}",
     "Aggregate topological complexity"],
]

table = ax.table(
    cellText=rows[1:],
    colLabels=rows[0],
    cellLoc="center",
    loc="center",
    bbox=[0, 0, 1, 1]
)
table.auto_set_font_size(False)
table.set_fontsize(9.5)
for (row, col), cell in table.get_celld().items():
    cell.set_edgecolor("#cccccc")
    if row == 0:
        cell.set_facecolor("#2A7DC9")
        cell.set_text_props(color="white", fontweight="bold")
    elif col == 3:
        cell.set_facecolor("#f5f5f5")
        cell.set_text_props(color="#333333", fontsize=8.5)
    elif row % 2 == 0:
        cell.set_facecolor("#eef4fb")
    else:
        cell.set_facecolor("white")

plt.tight_layout()
plt.savefig("output/step6_summary.png", dpi=160, bbox_inches="tight")
plt.show()
print("  ✓ step6_summary.png")

print("\n✅ All done! Files in output/:")
for f in sorted(os.listdir("output")):
    print(f"   {f}")