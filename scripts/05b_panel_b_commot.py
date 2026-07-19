#!/usr/bin/env python
"""
Per-sample COMMOT (optimal-transport CCC) on Xenium data — spatial-aware
cross-validation of CellChat Panel B findings.

Input:  data/processed/commot/per_sample_input/{sample_id}/{counts.mtx,
        logcounts.mtx, genes.csv, cells.csv} (prepared from raw Xenium data;
        prep step not included in this pipeline)
Output: data/processed/commot/per_sample/{sample_id}/{commot_results.csv,
        summary.txt} — consumed by scripts/05c_panel_b_commot_aggregate.R

Parameters:
  - L-R database: CellChatDB (same as CellChat for direct comparison)
  - species: human
  - dis_thr (distance threshold): 250 µm (matches CellChat interaction.range)
  - heteromeric: True (handles receptor complexes)

"""

import os
import sys
import time
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.io as sio
import anndata as ad
import commot as ct

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import PROCESSED_DIR
BASE_DIR = str(PROCESSED_DIR)
INPUT_DIR  = os.path.join(BASE_DIR, "commot", "per_sample_input")
OUTPUT_DIR = os.path.join(BASE_DIR, "commot", "per_sample")
os.makedirs(OUTPUT_DIR, exist_ok=True)

SAMPLES = [
    ("18321_a", "irHepatitis"),
    ("24774_d", "irHepatitis"),
    ("54023_b", "AIH"),
    ("56784_c", "AIH"),
]

# COMMOT parameters
DIS_THR = 250.0  # µm — match CellChat interaction.range
SPECIES = "human"
DB_NAME = "CellChat"    # CellChat = CellChatDB in COMMOT's naming (for same L-R pairs as CellChat)


def load_sample(sample_id):
    """Load per-sample AnnData from exported R data."""
    sdir = os.path.join(INPUT_DIR, sample_id)
    counts = sio.mmread(os.path.join(sdir, "counts.mtx")).tocsr().T  # cells × genes
    logcounts = sio.mmread(os.path.join(sdir, "logcounts.mtx")).tocsr().T
    genes = pd.read_csv(os.path.join(sdir, "genes.csv"))["gene"].tolist()
    cells = pd.read_csv(os.path.join(sdir, "cells.csv"))

    adata = ad.AnnData(
        X=logcounts,
        obs=cells.set_index("cell_ID"),
        var=pd.DataFrame(index=genes),
    )
    adata.layers["counts"] = counts

    # Spatial coordinates required by COMMOT
    adata.obsm["spatial"] = cells[["sdimx", "sdimy"]].values.astype(float)

    return adata


def run_one_sample(sample_id, condition):
    print(f"\n========== {sample_id} ({condition}) ==========", flush=True)
    t0 = time.time()

    adata = load_sample(sample_id)
    print(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes", flush=True)

    # Load CellChatDB L-R pairs (filtered to genes in panel)
    df_lr = ct.pp.ligand_receptor_database(
        species=SPECIES, signaling_type="Secreted Signaling", database=DB_NAME
    )
    df_lr_cc = ct.pp.ligand_receptor_database(
        species=SPECIES, signaling_type="Cell-Cell Contact", database=DB_NAME
    )
    df_lr_ecm = ct.pp.ligand_receptor_database(
        species=SPECIES, signaling_type="ECM-Receptor", database=DB_NAME
    )
    # Combine all categories
    df_lr_all = pd.concat([df_lr, df_lr_cc, df_lr_ecm], ignore_index=True)
    print(f"  CellChatDB L-R pairs (all categories): {len(df_lr_all)}", flush=True)

    # Filter to pairs where both ligand and receptor are in panel
    df_lr_filt = ct.pp.filter_lr_database(df_lr_all, adata, min_cell_pct=0.05)
    print(f"  After filtering to panel genes with min 5% expression: {len(df_lr_filt)}", flush=True)

    # Run COMMOT
    print(f"  Running ct.tl.spatial_communication (dis_thr={DIS_THR})...", flush=True)
    ct.tl.spatial_communication(
        adata,
        database_name=DB_NAME,
        df_ligrec=df_lr_filt,
        dis_thr=DIS_THR,
        heteromeric=True,
        pathway_sum=True,
    )
    print(f"  Done in {(time.time() - t0) / 60:.1f} min", flush=True)

    # Extract cell-pair × L-R results
    # COMMOT stores in adata.obsp with keys like "commot-CellChatDB-pathway" etc.
    # We want per-cell L-R communication scores and cell-type-aggregate summaries
    out_sdir = os.path.join(OUTPUT_DIR, sample_id)
    os.makedirs(out_sdir, exist_ok=True)

    # Save full adata
    adata.write_h5ad(os.path.join(out_sdir, "adata_with_commot.h5ad"))

    # Cell-type aggregation: sum signals by source → target cell-type pairs
    cell_labels = adata.obs["cluster_name"].astype(str)
    unique_labels = sorted(cell_labels.unique())

    # Extract per-L-R per-cell-type-pair aggregate
    records = []
    for key in list(adata.obsp.keys()):
        if not key.startswith(f"commot-{DB_NAME}-"):
            continue
        lr_id = key.replace(f"commot-{DB_NAME}-", "")
        # adata.obsp[key] is cell × cell sparse matrix of communication flows
        # Aggregate by cell-type pair
        mat = adata.obsp[key]
        for src in unique_labels:
            src_idx = np.where(cell_labels == src)[0]
            if len(src_idx) < 10:
                continue
            for tgt in unique_labels:
                tgt_idx = np.where(cell_labels == tgt)[0]
                if len(tgt_idx) < 10:
                    continue
                sub = mat[src_idx, :][:, tgt_idx]
                total_flow = float(sub.sum())
                if total_flow > 0:
                    records.append({
                        "lr_id": lr_id,
                        "source": src,
                        "target": tgt,
                        "total_flow": total_flow,
                        "mean_flow": total_flow / (len(src_idx) * len(tgt_idx)),
                        "n_source": len(src_idx),
                        "n_target": len(tgt_idx),
                    })

    res = pd.DataFrame(records)
    res["sample_id"] = sample_id
    res["condition"] = condition
    res.to_csv(os.path.join(out_sdir, "commot_results.csv"), index=False)
    print(f"  Saved commot_results.csv ({len(res)} rows)", flush=True)

    # Summary
    with open(os.path.join(out_sdir, "summary.txt"), "w") as f:
        f.write(f"Sample: {sample_id} ({condition})\n")
        f.write(f"Cells: {adata.n_obs}\n")
        f.write(f"Genes (panel): {adata.n_vars}\n")
        f.write(f"L-R pairs tested: {len(df_lr_filt)}\n")
        f.write(f"Distance threshold: {DIS_THR} µm\n")
        f.write(f"Runtime (min): {(time.time() - t0) / 60:.1f}\n")
        f.write(f"Cell-type-pair × L-R records with non-zero flow: {len(res)}\n")

    return res


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", default=None,
                        help="Run only this sample_id; otherwise iterate over all.")
    args = parser.parse_args()

    if args.sample is not None:
        matches = [(s, c) for s, c in SAMPLES if s == args.sample]
        if not matches:
            raise SystemExit(f"Unknown sample '{args.sample}'. Known: {[s for s,_ in SAMPLES]}")
        sid, cond = matches[0]
        try:
            run_one_sample(sid, cond)
        except Exception as e:
            print(f"  ERROR on {sid}: {e}", flush=True)
            import traceback; traceback.print_exc()
            raise
        print("\nDone.")
    else:
        all_results = []
        for sid, cond in SAMPLES:
            try:
                res = run_one_sample(sid, cond)
                all_results.append(res)
            except Exception as e:
                print(f"  ERROR on {sid}: {e}", flush=True)
                import traceback; traceback.print_exc()

        if all_results:
            master = pd.concat(all_results, ignore_index=True)
            master.to_csv(os.path.join(OUTPUT_DIR, "commot_all_samples.csv"), index=False)
            print(f"\nWrote master: {len(master)} rows across {len(all_results)} samples")

        print("\nDone.")
