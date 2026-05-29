import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import rasterio
import gudhi
import matplotlib.patches as mpatches
import scipy.ndimage as ndimage
from scipy.spatial import cKDTree

MAX_SIZE = 2000

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
    

def plot_data_binary(data_dict, cities_list, output_path="output/combined_missing_data_mask.png"):
    valid_cities = [c for c in cities_list if c in data_dict]
    num_cities = len(valid_cities)

    fig, axes = plt.subplots(1, num_cities, figsize=(9 * num_cities, 11), facecolor='white')

    if num_cities == 1:
        axes = [axes]

    cmap_holes = mcolors.ListedColormap(['#111111', '#FF3333'])

    # Loop over the axes and matching city data simultaneously
    for ax, city in zip(axes, valid_cities):
        roi_data = data_dict[city]

        # Auto-downsample before building the mask.
        # matplotlib converts the boolean array to an RGBA float64 buffer
        # (4 channels × pixels × 8 bytes) when rendering: a full 25000×10000
        # array would require ~7.5 GiB. Capping the longest axis at 2000 px
        # keeps the render buffer under ~130 MB.
        # NaN % statistics are still computed on the full-resolution array.
        vis_factor = max(1, max(roi_data.shape) // 2000)
        missing_mask = np.isnan(roi_data[::vis_factor, ::vis_factor])

        # Calculate statistical attributes on the full array for accuracy
        total_pixels = roi_data.size
        missing_pixels = int(np.isnan(roi_data).sum())
        missing_percentage = (missing_pixels / total_pixels) * 100

        ax.imshow(missing_mask, cmap=cmap_holes)

        # Build legend configurations for each individual subplot channel
        dark_patch = mpatches.Patch(color='#111111', label='Available elevation data')
        red_patch = mpatches.Patch(color='#FF3333', label='Missing data (NoData / LiDAR Holes)')
        ax.legend(handles=[dark_patch, red_patch], loc='upper right', fontsize=10)

        # Title strings with embedded formatting
        ax.set_title(f"{city}\nData Gaps (NoData: {missing_percentage:.3f}%)", fontsize=13, pad=12)
        ax.axis("off")

    # Tighten up structural bounds to avoid overlapping margins or label truncations
    fig.tight_layout()

    # Save — dpi=150 is sufficient for full-raster overview plots
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=150, facecolor='white', bbox_inches='tight')
    print(f"\n[SUCCESS] Saved unified column visualization to: {output_path}\n")

    # Force Matplotlib to output and render the image inline within the Jupyter Notebook cell
    plt.show()


def plot_topological_grids(data_dict, step_title="Terrain Elevation", vis_factor=None, save_filename=None, cmap_name="gist_earth"):
    """
    Unified, memory-safe plotting framework for urban DTM/DSM topological datasets.
    Keeps legends, scale ranges, and description stickers standardized across all targets
    for precise cross-city comparisons.
    
    Parameters:
    ----------
    data_dict : dict
        Dictionary where keys are city names/identifiers and values are the 2D arrays.
    step_title : str
        The subtitle indicator used to track your pipeline progression stage.
    vis_factor : int or None
        Downsampling index step to prevent MemoryErrors.
    save_filename : str or None
        File path destination for the PNG asset.
    cmap_name : str
        The Matplotlib colormap string (e.g., 'gist_earth', 'terrain', 'Earth').
    """
    keys = list(data_dict.keys())
    num_plots = len(keys)
    
    # 1. Automatically calculate dynamic width aspect ratios
    width_ratios = []
    for k in keys:
        shape = data_dict[k].shape
        width_ratios.append(shape[1] / shape[0])
        
    # 2. Establish figure architecture
    if num_plots > 1:
        fig, axes = plt.subplots(1, num_plots, figsize=(10 * num_plots, 9), 
                                 gridspec_kw={"width_ratios": width_ratios}, facecolor='white')
        if num_plots == 2:
            main_title = f"{step_title}\n{keys[0]} vs {keys[1]}"
        else:
            main_title = f"{step_title}\nComparative Dataset Analysis"
    else:
        fig, axes = plt.subplots(1, 1, figsize=(10, 20), facecolor='white')
        axes = [axes]
        main_title = f"{keys[0]} DTM Layout\n{step_title}"
        
    fig.suptitle(main_title, fontsize=14, fontweight="bold", y=1.02)
    
    # 3. Iterate over each array configuration element symmetrically
    for ax, city_name in zip(axes, keys):
        grid_data = data_dict[city_name]
        
        if vis_factor is None:
            current_vis = 4 if max(grid_data.shape) > 10000 else 1
        else:
            current_vis = vis_factor
            
        render_buffer = grid_data[::current_vis, ::current_vis]
        
        # Calculate consistent metadata statistics for the summary sticker
        nan_clean = render_buffer[~np.isnan(render_buffer)]
        mean_elev = np.mean(nan_clean) if len(nan_clean) > 0 else 0.0
        vmin_use = np.nanmin(grid_data)
        vmax_use = np.nanmax(grid_data)
        
        desc_text = f"Mean: {mean_elev:+.2f}m NAP\nLocal Range: [{vmin_use:.1f}m, {vmax_use:.1f}m]"
            
        # Display the terrain matrix array using a unified scale
        im = ax.imshow(render_buffer, cmap=cmap_name, vmin=vmin_use, vmax=vmax_use)
        ax.set_title(city_name, fontsize=13, fontweight="bold", pad=8)
        ax.axis("off")
        
        # Standardized colorbar scaling for direct, 1-to-1 comparison tracking
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink=0.8 if num_plots > 1 else 1.0)
        cbar.set_label("Elevation (m)", fontsize=11, labelpad=8)
            
        # Inject matching information boxes across all target axis panels
        try:
            annotate_box(ax, desc_text)
        except NameError:
            ax.text(0.05, 0.05, desc_text, transform=ax.transAxes, color="black",
                    bbox=dict(facecolor="white", alpha=0.8, boxstyle="round,pad=0.5"), fontsize=10)

    plt.tight_layout()
    
    if save_filename:
        os.makedirs("output", exist_ok=True)
        plt.savefig(f"output/{save_filename}", dpi=300, facecolor='white', bbox_inches='tight')
        print(f"Publication-ready canvas successfully written out to output/{save_filename}")
        
    plt.show()