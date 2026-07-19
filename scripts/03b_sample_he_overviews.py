#!/usr/bin/env python3
# ============================================================================
# 03b_sample_he_overviews.py
# Generate full-sample H&E overview PNGs (warped into Giotto coordinate space)
# ============================================================================
#
# Called by 03_figure_1_main.R to produce H&E panels for Figure 1.
# Reads a JSON config with per-sample spatial extents, loads H&E images
# at a low pyramid level, and applies the affine alignment.
#
# Usage: python3 03b_sample_he_overviews.py config.json
#
# Python dependencies (install once):
#   pip install numpy tifffile scipy Pillow
# ============================================================================

import sys
import json
import os
import csv
import gc
import numpy as np
import tifffile
from PIL import Image
from scipy.ndimage import affine_transform

PIXEL_SIZE = 0.2125
HE_LEVEL = 2          # 1/4 resolution (good balance for overview)
OUTPUT_WIDTH = 1200    # pixels wide (height computed from aspect ratio)

SAMPLE_OFFSETS = {
    "18321_a": (0, 0),
    "24774_d": (12000, 0),
    "54023_b": (0, -8000),
    "56784_c": (12000, -8000)
}


def load_alignment_matrix(path):
    with open(path) as f:
        return np.array([list(map(float, row)) for row in csv.reader(f)])


def generate_overview(sample_id, working_dir, output_dir, extent):
    """Generate an H&E overview for one sample.

    extent: dict with xmin, xmax, ymin, ymax in Giotto µm coordinates.
    """
    x_off, y_off = SAMPLE_OFFSETS[sample_id]
    sample_path = os.path.join(working_dir, f"0027420_{sample_id}")
    he_file = os.path.join(sample_path, f"{sample_id}_he_image.ome.tif")
    matrix_file = os.path.join(sample_path, "he_alignment", "matrix.csv")

    if not os.path.exists(he_file) or not os.path.exists(matrix_file):
        print(f"  {sample_id}: H&E or alignment matrix not found, skipping")
        return

    raw_affine = load_alignment_matrix(matrix_file)
    scale_factor = 2 ** HE_LEVEL

    print(f"  {sample_id}: loading H&E level {HE_LEVEL}...", end=" ", flush=True)
    with tifffile.TiffFile(he_file) as tif:
        he_img = tif.series[0].levels[HE_LEVEL].asarray()
    print(f"shape={he_img.shape}")

    # Output dimensions from spatial extent
    xmin, xmax = extent["xmin"], extent["xmax"]
    ymin, ymax = extent["ymin"], extent["ymax"]
    width_um = xmax - xmin
    height_um = ymax - ymin

    out_w = OUTPUT_WIDTH
    out_h = int(out_w * height_um / width_um)

    # Build affine: output pixel (py, px) -> H&E level pixel (row, col)
    # Same math as 07a_panel_c_he_extract.py extract_he_crop()
    inv_aff = np.linalg.inv(raw_affine)
    ia = inv_aff[:2, :2]
    it = inv_aff[:2, 2]

    step_x = width_um / out_w
    step_y = height_um / out_h
    ps = PIXEL_SIZE
    sf = scale_factor

    # Output pixel (0,0) = top-left = Giotto (xmin, ymax)
    x0 = xmin - x_off
    y0 = -(ymax) + y_off

    # Matrix maps [row_out, col_out] -> [row_in, col_in]
    # row_out increments y downward (Giotto y decreasing)
    # col_out increments x rightward (Giotto x increasing)
    matrix = np.array([
        [ia[1, 1] * step_y / (ps * sf), ia[1, 0] * step_x / (ps * sf)],
        [ia[0, 1] * step_y / (ps * sf), ia[0, 0] * step_x / (ps * sf)]
    ])
    offset = np.array([
        (ia[1, 0] * x0 + ia[1, 1] * y0) / (ps * sf) + it[1] / sf,
        (ia[0, 0] * x0 + ia[0, 1] * y0) / (ps * sf) + it[0] / sf
    ])

    out_img = np.full((out_h, out_w, 3), 255, dtype=np.uint8)
    for c in range(3):
        out_img[:, :, c] = affine_transform(
            he_img[:, :, c], matrix, offset,
            output_shape=(out_h, out_w),
            order=1, mode='constant', cval=255
        )

    out_path = os.path.join(output_dir, f"he_overview_{sample_id}.png")
    Image.fromarray(out_img).save(out_path, optimize=True)
    print(f"  {sample_id}: saved {out_path} ({out_w}x{out_h})")

    del he_img
    gc.collect()


def main(config_file):
    with open(config_file) as f:
        config = json.load(f)

    working_dir = config["working_dir"]
    output_dir = config["output_dir"]
    os.makedirs(output_dir, exist_ok=True)

    print("Generating H&E sample overviews...")
    for sample in config["samples"]:
        generate_overview(
            sample["sample_id"], working_dir, output_dir, sample["extent"]
        )
    print("Done.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <config.json>")
        sys.exit(1)
    main(sys.argv[1])
