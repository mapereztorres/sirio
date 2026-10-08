import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
import re
import matplotlib as mpl

beta_values = np.logspace(-4, -2, 50)
#beta0 = 1e-3
beta0=10**((np.log10(BETA_EFF_MAX) + np.log10(BETA_EFF_MIN))/2)
#if 'LHS1610' in starname:

#Reference value for detection threshold
sigma3 = 0.1 # mJy

print('starname')
print(starname)
print(sigma3)

models = [
    ("alfven_wing", "_alfven_wing_model.csv",   "Alfvén Wing"),
    ("reconnection", "_reconnection_model.csv", "Reconnection"),
    ("sb", "_sb_model.csv",                     "Stretch and Break")
]

for tag, suffix, title_text in models:

    print(f"\n=== Processing {tag} model ===")

    csv_file = os.path.join(FOLDER, "CSV", outfile + suffix)
    print("Reading CSV:", csv_file)

    if not os.path.isfile(csv_file):
        print(" -> File not found, skipping.")
        continue

    df = pd.read_csv(csv_file)
    print("Loaded CSV with", len(df), "rows")

    df["M_DOT"] = pd.to_numeric(df["M_DOT"], errors="coerce")
    df["FLUX"]  = pd.to_numeric(df["FLUX"],  errors="coerce")

    mdot_vals = np.sort(df["M_DOT"].dropna().unique().astype(float))
    x_super = float(x_superalfv)

    # BUILD FLUX GRID
    flux_grid = np.zeros((len(beta_values), len(mdot_vals)))
    for j, mdot in enumerate(mdot_vals):
        flux0 = df.loc[df["M_DOT"] == mdot, "FLUX"].values[0]
        flux_grid[:, j] = flux0 * (beta_values / beta0)

    # APPLY MASK (grey out super-alfvenic region)
    mask = (mdot_vals > x_super)[np.newaxis, :].repeat(len(beta_values), axis=0)
    masked_flux = np.ma.array(flux_grid, mask=mask)

    # PLOT
    fig, ax = plt.subplots(figsize=(8, 6))

    c = ax.pcolormesh(
        mdot_vals, beta_values, masked_flux,
        shading="auto",
        norm=LogNorm(vmin=1e-2, vmax=1e2),
        linewidth=0, edgecolors="none", zorder=0
    )
    c.set_rasterized(True)

    # COLORBAR
    cb = plt.colorbar(c, ax=ax, label="Flux (mJy)")
    cb.ax.axhline(sigma3, color="white", linewidth=3)

    # --- Hatching on the colorbar for the excluded region ---
    cb_ymin = cb.norm.vmin   # 1e-2
    cb_ymax = cb.norm.vmax   # 1e2

    frac_sigma3 = (np.log10(sigma3) - np.log10(cb_ymin)) / (np.log10(cb_ymax) - np.log10(cb_ymin))
    frac_sigma3 = np.clip(frac_sigma3, 0, 1)

    mpl.rcParams['hatch.linewidth'] = 5.5
    hatch_rect = Rectangle(
        (0, 0), 1, frac_sigma3,
        transform=cb.ax.transAxes,
        hatch="/",
        facecolor="none",
        edgecolor="white",
        linewidth=0,
        zorder=5,
        clip_on=True
    )
    cb.ax.add_patch(hatch_rect)

    # --- "Undetectable" label to the right of the colorbar ---
    cb.ax.annotate(
        "Undetectable",
        xy=(1.0, frac_sigma3 / 2),
        xycoords="axes fraction",
        xytext=(1.8, frac_sigma3 / 2),
        textcoords="axes fraction",
        fontsize=13,
        color="black",
        ha="left", va="center",
        arrowprops=dict(arrowstyle="-", color="black", lw=1.2),
        annotation_clip=False,
    )

    # CONTOURS
    levels = [1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4]
    cs2 = ax.contour(
        mdot_vals, beta_values, masked_flux,
        levels=levels,
        colors="black", linestyles="dashed", linewidths=4, zorder=3
    )
    cs1 = ax.contour(
        mdot_vals, beta_values, masked_flux,
        levels=[sigma3],
        colors="white", linestyles="solid", linewidths=4, zorder=4
    )

    def label_contours_at_midpoint(ax, cs, fmt, fontsize=20, color="white", threshold=None):
        for i, level in enumerate(cs.levels):
            # matplotlib >= 3.8: .collections fue eliminado; usar allsegs
            try:
                segs = cs.allsegs[i]
            except AttributeError:
                segs = [path.vertices for path in cs.collections[i].get_paths()]
            for seg in segs:
                if len(seg) == 0:
                    continue
                mid = seg[len(seg) // 2]
                if threshold is not None and level < threshold:
                    bbox = dict(boxstyle="round,pad=0.2", fc="black", alpha=0.8, ec="none")
                else:
                    bbox = dict(boxstyle="round,pad=0.2", fc="none", alpha=0.0, ec="none")
                ax.text(mid[0], mid[1], fmt(level),
                        fontsize=fontsize, color=color,
                        ha="center", va="center", rotation=0, bbox=bbox)

    label_contours_at_midpoint(ax, cs2, fmt=lambda x: f"{x:g} mJy", color="white", threshold=sigma3)

    ax.set_xscale("log")
    ax.set_yscale("log")

    # =========================================================
    # HATCHED FILL — white stripes over regions where flux < sigma3
    # =========================================================
    mpl.rcParams['hatch.linewidth'] = 5.5

    cf = ax.contourf(
        mdot_vals, beta_values, masked_flux.filled(np.inf),
        levels=[0, sigma3],
        hatches=["/"],
        colors=["none"],
        zorder=0.1
    )
    for collection in cf.collections:
        collection.set_edgecolor("white")
        collection.set_linewidth(0)

    # --- "Undetectable" label inside the hatched region of the main plot ---

    valid_mdot = mdot_vals[mdot_vals <= x_super]

    if len(valid_mdot) == 0:
        print(f" -> No sub-Alfvénic mdot values for {tag} ({starname}); skipping 'Undetectable' label.")
    else:
        label_x = np.sqrt(valid_mdot[0] * valid_mdot[-1])   # geometric mean in x

        j_mid = np.argmin(np.abs(mdot_vals - label_x))
        flux_col = flux_grid[:, j_mid]
        below_sigma = np.where(flux_col < sigma3)[0]

        # Only draw if the hatched region is large enough and label_x is inside the plot
        if len(below_sigma) >= 3 and label_x > mdot_vals[0] * 2:
            bottom_half = below_sigma[:len(below_sigma) // 2]
            i_label = bottom_half[len(bottom_half) // 2]
            label_y = beta_values[i_label]

            ax.text(
                label_x, label_y,
                "Undetectable",
                fontsize=14,
                color="white",
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.3", fc="black", alpha=1.0, ec="none"),
                zorder=20
            )

    # LABELS / TITLE
    ax.set_xlabel(r"Mass Loss rate [$\dot{M}_\odot$]")
    ax.set_ylabel(r"$\beta$")

    geom2 = geometry.strip("-").replace("-", " ")
    geom2 = re.sub(r"Bstar[0-9.\[\]A-Za-z]*", "", geom2)
    geom2 = " ".join(w.capitalize() for w in geom2.split())
    if "Pfss" in geom2:
        geom2 = "PFSS"

    ax.set_title(f"{title_text}: {geom2}")

    # SAVE
    outdir = os.path.join(FOLDER, "contours")
    os.makedirs(outdir, exist_ok=True)

    pdf_path = os.path.join(outdir, outfile + suffix.replace("_model.csv", "_contours.pdf"))
    png_path = os.path.join(outdir, outfile + suffix.replace("_model.csv", "_contours.png"))

    fig.canvas.draw()
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Saved figures:")
    print(" ", pdf_path)
    print(" ", png_path)
