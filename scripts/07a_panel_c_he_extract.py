#!/usr/bin/env python3
# ============================================================================
# 07a_panel_c_he_extract.py
# Extract H&E and DAPI image crops for supplementary ROI figures
# ============================================================================
#
# PURPOSE:
#   Reads ROI specifications from a JSON config file (produced by
#   07_panel_c_he_extract.R) and extracts aligned H&E and DAPI image
#   crops for each ROI. Handles the affine alignment of H&E images
#   and coordinate transforms between Giotto and Xenium spaces.
#
# USAGE:
#   python3 07a_panel_c_he_extract.py roi_config.json
#
# DEPENDENCIES:
#   tifffile, numpy, scipy, PIL (Pillow)
#
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

# --- Constants ---
PIXEL_SIZE = 0.2125       # Xenium pixel size in microns
HE_PYRAMID_LEVEL = 1     # Half-resolution H&E (good quality/memory balance)
OUTPUT_PX = 600           # Output crop size in pixels

SAMPLE_OFFSETS = {
    "18321_a": (0, 0),
    "24774_d": (12000, 0),
    "54023_b": (0, -8000),
    "56784_c": (12000, -8000)
}


def load_alignment_matrix(matrix_file):
    """Load 3x3 affine matrix from Xenium he_alignment/matrix.csv."""
    with open(matrix_file) as f:
        reader = csv.reader(f)
        return np.array([list(map(float, row)) for row in reader])


def extract_dapi_crop(dapi_img, roi_x, roi_y, roi_size, x_off, y_off):
    """Extract and normalize a DAPI crop.

    DAPI pixels map directly to Xenium coordinates:
      x_pixel = (x_giotto - x_off) / PIXEL_SIZE
      y_pixel = -(y_giotto - y_off) / PIXEL_SIZE
    """
    half = roi_size / 2
    cx_g = roi_x + half
    cy_g = roi_y + half

    # Giotto -> Xenium local (um) -> pixels
    cx_local = cx_g - x_off
    cy_local = -(cy_g - y_off)
    cx_px = cx_local / PIXEL_SIZE
    cy_px = cy_local / PIXEL_SIZE
    half_px = half / PIXEL_SIZE

    # Crop bounds (row = y, col = x)
    r_min = max(0, int(round(cy_px - half_px)))
    r_max = min(dapi_img.shape[0], int(round(cy_px + half_px)))
    c_min = max(0, int(round(cx_px - half_px)))
    c_max = min(dapi_img.shape[1], int(round(cx_px + half_px)))

    crop = dapi_img[r_min:r_max, c_min:c_max]
    if crop.size == 0:
        return None

    # Normalize uint16 to uint8 with contrast stretching
    p1, p99 = np.percentile(crop, [1, 99.5])
    norm = np.clip((crop.astype(np.float32) - p1) / max(p99 - p1, 1) * 255,
                   0, 255).astype(np.uint8)

    img = Image.fromarray(norm)
    return img.resize((OUTPUT_PX, OUTPUT_PX), Image.LANCZOS)


def extract_he_crop(he_img, raw_affine, roi_x, roi_y, roi_size,
                    x_off, y_off, scale_factor):
    """Extract an H&E crop using affine warp into Giotto coordinate space.

    raw_affine maps H&E level-0 pixels -> Xenium pixels.
    scale_factor: pyramid downscale (e.g. 2 for level 1).

    The transform chain for each output pixel (py_out, px_out):
      Giotto um -> Xenium local um -> Xenium px -> H&E px -> H&E level px
    is computed as a single 2x2 affine + offset for scipy.
    """
    inv_aff = np.linalg.inv(raw_affine)
    ia = inv_aff[:2, :2]   # 2x2
    it = inv_aff[:2, 2]    # translation

    step = roi_size / OUTPUT_PX
    ps = PIXEL_SIZE
    sf = scale_factor
    k = step / (ps * sf)

    # x0, y0: Xenium-local um at the output origin
    # Output pixel (0,0) = top-left = Giotto (roi_x, roi_y + roi_size)
    x0 = roi_x - x_off
    y0 = -(roi_y + roi_size) + y_off

    # scipy affine_transform: input[row,col] = matrix @ output[row,col] + offset
    # row = H&E_y, col = H&E_x; output row = py_out, col = px_out
    matrix = np.array([
        [ia[1, 1] * k, ia[1, 0] * k],
        [ia[0, 1] * k, ia[0, 0] * k]
    ])
    offset = np.array([
        (ia[1, 0] * x0 + ia[1, 1] * y0) / (ps * sf) + it[1] / sf,
        (ia[0, 0] * x0 + ia[0, 1] * y0) / (ps * sf) + it[0] / sf
    ])

    out_img = np.zeros((OUTPUT_PX, OUTPUT_PX, 3), dtype=np.uint8)
    for c in range(3):
        out_img[:, :, c] = affine_transform(
            he_img[:, :, c],
            matrix,
            offset,
            output_shape=(OUTPUT_PX, OUTPUT_PX),
            order=1,        # bilinear interpolation
            mode='constant',
            cval=255        # white background for out-of-bounds
        )

    return Image.fromarray(out_img)


def main(config_file):
    with open(config_file) as f:
        config = json.load(f)

    working_dir = config["working_dir"]
    output_dir = config["output_dir"]
    rois = config["rois"]

    os.makedirs(output_dir, exist_ok=True)

    # Group ROIs by sample for efficient I/O
    by_sample = {}
    for roi in rois:
        by_sample.setdefault(roi["sample_id"], []).append(roi)

    for sample_id, sample_rois in by_sample.items():
        x_off, y_off = SAMPLE_OFFSETS[sample_id]
        sample_path = os.path.join(working_dir, f"0027420_{sample_id}")

        print(f"\n{'=' * 60}")
        print(f"  {sample_id}: {len(sample_rois)} ROIs")
        print(f"{'=' * 60}")

        # --- Load DAPI (plain TIFF, ~730 MB) ---
        dapi_file = os.path.join(sample_path, "morphology_focus",
                                 "tif_exports", "morphology_focus_0000.tif")
        dapi_img = None
        if os.path.exists(dapi_file):
            print("  Loading DAPI...", end=" ", flush=True)
            dapi_img = tifffile.imread(dapi_file)
            print(f"shape={dapi_img.shape}")
        else:
            print("  WARNING: DAPI not found")

        # --- Load H&E (pyramid level, ~1 GB for level 1) ---
        he_file = os.path.join(sample_path,
                               f"{sample_id}_he_image.ome.tif")
        he_img = None
        raw_affine = None
        he_scale = 2 ** HE_PYRAMID_LEVEL

        if os.path.exists(he_file):
            matrix_file = os.path.join(sample_path,
                                       "he_alignment", "matrix.csv")
            if os.path.exists(matrix_file):
                raw_affine = load_alignment_matrix(matrix_file)

            print(f"  Loading H&E level {HE_PYRAMID_LEVEL}...", end=" ",
                  flush=True)
            with tifffile.TiffFile(he_file) as tif:
                he_img = tif.series[0].levels[HE_PYRAMID_LEVEL].asarray()
            print(f"shape={he_img.shape}")
        else:
            print("  WARNING: H&E not found")

        # --- Extract crops ---
        for roi in sample_rois:
            cl = roi["cluster_id"]
            rx, ry = roi["roi_x"], roi["roi_y"]
            rs = roi["roi_size"]
            print(f"\n  Cluster {cl:2d}: ROI ({rx:.0f}, {ry:.0f}), "
                  f"size {rs} um")

            if dapi_img is not None:
                dapi_crop = extract_dapi_crop(dapi_img, rx, ry, rs,
                                              x_off, y_off)
                if dapi_crop is not None:
                    p = os.path.join(output_dir, f"cluster_{cl:02d}_dapi.png")
                    dapi_crop.save(p)
                    print(f"    DAPI -> {os.path.basename(p)}")

            if he_img is not None and raw_affine is not None:
                he_crop = extract_he_crop(he_img, raw_affine, rx, ry, rs,
                                          x_off, y_off, he_scale)
                p = os.path.join(output_dir, f"cluster_{cl:02d}_he.png")
                he_crop.save(p)
                print(f"    H&E  -> {os.path.basename(p)}")

        del dapi_img, he_img
        gc.collect()

    print(f"\n{'=' * 60}")
    print("  All crops extracted successfully.")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <roi_config.json>")
        sys.exit(1)
    main(sys.argv[1])
