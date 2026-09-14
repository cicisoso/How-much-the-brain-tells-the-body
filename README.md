# How much the brain tells the body: the structural bandwidth of the *Drosophila* neck connective

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22758902.svg)](https://doi.org/10.5281/zenodo.22758902)

Code, per-axon tables and figures for the manuscript *"How much the brain tells the body: the structural bandwidth of the Drosophila neck connective"* (Liu, Zong, Chen and Xiong; School of Information Technology, Zhejiang Financial College).

We measured the cross-sectional area of every brain–body axon in the neck connective of the male *Drosophila* central nervous system connectome (MaleCNS v1.0) directly from the 16 nm segmentation, validated the measurement (8 nm re-measurement, bilateral homologues, mesh sections, giant fiber), and repeated it on the female BANC connectome by sectioning neuron meshes at the curated neck plane. The tables in `data_tables/` are Table S1 of the paper.

## Contents

| Folder | What it holds |
|---|---|
| `code/` | Analysis and figure scripts (Python 3.11). See the pipeline below. |
| `data_tables/` | Table S1: per-axon calibers, classes, targets and model delays (`Table_S1.xlsx` with README/MaleCNS/BANC sheets; CSV copies). |
| `results/` | Intermediate results used by the figure scripts (JSON summaries and parquet tables). |
| `figures/` | Main Figures 1–4 and Figure S1 (PNG and PDF). |

## Data sources (not included; public)

- MaleCNS v1.0: https://male-cns.janelia.org (CC-BY 4.0). Files used: `body-annotations-male-cns-v1.0-minconf-0.5.feather`, `body-stats-male-cns-v1.0-minconf-0.5.feather`, `syn-partners-male-cns-v1.0-minconf-0.5.feather`, SWC skeletons (`gs://flyem-male-cns/v1.0/segmentation/skeletons-malecns/skeletons-swc/`), and the segmentation (`gs://flyem-male-cns/v1.0/segmentation`, read with cloud-volume).
- BANC v888: public bucket `gs://lee-lab_brain-and-nerve-cord-fly-connectome` (files `compiled_data/banc_888/banc_888_meta.feather`, `neuron_annotations/v888/neck_connective_y92500.parquet`, and the `neuron_meshes` precomputed mesh source).

## Setup

```bash
pip install -r requirements.txt
export FLYBRAIN_DATA=/path/to/data   # folder holding the MaleCNS feather files, swc/ skeletons and banc/ files
```

Scripts resolve the data folder from `FLYBRAIN_DATA` (default `../data`) and write to `../results` and `figures/` relative to `code/`.

## Pipeline (run from `code/`)

1. `neck_waist.py` – locate the connective along the body axis from skeleton spread.
2. `male_neck_slabs.py` – download three 16 nm slabs of the segmentation and count voxels per segment and plane.
3. `caliber_analysis.py` – per-axon areas, tilt correction, per-body caliber table.
4. `male_neck_slab_mip0.py` + `mip_validation.py` – 8 nm validation slab.
5. `male_mesh_validation.py` – mesh-section versus voxel validation on 175 axons.
6. `delay_budget.py` – skeleton path lengths from the neck to brain and VNC arbors; conduction-delay model.
7. `budget_male.py` – information-budget tables and class summaries (writes `results/budget_male.json`).
8. `roi_shares.py` – synapse-level target neuropils of DNs and ANs from the partner table.
9. `banc_neck_sections.py` – female mesh sections at the BANC neck plane (multiprocessing; ~30 min).
10. `compare_sex.py`, `compare_sex2.py` – male–female comparison, relative calibers, Table S1 female sheet.
11. `fig1_census.py` … `fig4_sex.py`, `figS1_validation.py` – figures (`figlib.py` holds shared styling; the optional panel-alignment gate is skipped if the `audit_panel_alignment` module is absent).

## Citation

Archived version: Zenodo, https://doi.org/10.5281/zenodo.22758902. If you use these tables or code, please cite the Zenodo record and the paper (reference to be added on publication), together with the MaleCNS and BANC connectome papers listed in the manuscript.

## License

Code: MIT License. Data tables and figures: CC BY 4.0.
