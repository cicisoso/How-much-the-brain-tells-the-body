# The structural bandwidth of the neck connective between the *Drosophila* brain and body

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22794997.svg)](https://doi.org/10.5281/zenodo.22794997)

Code, per-axon tables, source data and figures for the manuscript *"The structural bandwidth of the neck connective between the Drosophila brain and body"* (Liu, Zong, Chen and Xiong; School of Information Technology, Zhejiang Financial College; submitted to *Communications Biology*). Version 1.0 of this repository accompanied an earlier, shorter version of the manuscript titled "How much the brain tells the body".

We measured the cross-sectional area of every brain–body axon in the neck connective of the male *Drosophila* central nervous system connectome (MaleCNS v1.0) directly from the 16 nm segmentation, validated the measurement four ways (8 nm re-measurement, bilateral homologues, mesh sections, giant fiber), related each axon's caliber to its synaptic inputs and outputs, its cell-type population size and its predicted neurotransmitter, and repeated the caliber census on the female BANC connectome by sectioning neuron meshes at the curated neck plane.

## Contents

| Folder | What it holds |
|---|---|
| `code/` | Analysis and figure scripts (Python 3.11). See the pipeline below. |
| `data_tables/` | Supplementary Data 1 (`Supplementary_Data_1.xlsx`: README, MaleCNS and BANC sheets; CSV copies of the two data sheets) and Source Data (`Source_Data.xlsx`: one sheet per figure panel and for Table 1). |
| `results/` | Intermediate results used by the figure scripts (JSON summaries, CSV tables and parquet tables, including `connectivity_scaling.json`, `robustness.json`, `budget_table.csv` and `channel_axons_male_connectivity.parquet`). |
| `figures/` | Figures 1–5 and Supplementary Figures 1–2 (PNG and PDF). |

## Data sources (not included; public)

- MaleCNS v1.0: https://male-cns.janelia.org (CC-BY 4.0). Files used: `body-annotations-male-cns-v1.0-minconf-0.5.feather`, `body-stats-male-cns-v1.0-minconf-0.5.feather`, `body-neurotransmitters-male-cns-v1.0.feather`, `syn-partners-male-cns-v1.0-minconf-0.5.feather`, SWC skeletons (`gs://flyem-male-cns/v1.0/segmentation/skeletons-malecns/skeletons-swc/`), and the segmentation (`gs://flyem-male-cns/v1.0/segmentation`, read with cloud-volume).
- BANC v888: public bucket `gs://lee-lab_brain-and-nerve-cord-fly-connectome` (files `compiled_data/banc_888/banc_888_meta.feather`, `neuron_annotations/v888/neck_connective_y92500.parquet`, and the `neuron_meshes` precomputed mesh source).

## Setup

```bash
pip install -r requirements.txt
export FLYBRAIN_DATA=/path/to/data   # folder holding the MaleCNS feather files, swc/ skeletons, neck_slabs/ and banc/ files
```

Scripts resolve the data folder from `FLYBRAIN_DATA` (default `../data`) and read and write `results/` and `figures/` relative to the working directory; run them from a folder that contains (or links to) the `results/` and `figures/` folders of this repository.

## Pipeline

1. `neck_waist.py` – locate the connective along the body axis from skeleton spread.
2. `male_neck_slabs.py` – download three 16 nm slabs of the segmentation and count voxels per segment and plane.
3. `caliber_analysis.py` – per-axon areas, tilt correction, per-body caliber table.
4. `male_neck_slab_mip0.py` + `mip_validation.py` – 8 nm validation slab.
5. `male_mesh_validation.py` – mesh-section versus voxel validation on 175 axons.
6. `delay_budget.py` – skeleton path lengths from the neck to brain and VNC arbors; conduction-delay model.
7. `budget_male.py` – information-budget tables and class summaries (`results/budget_male.json`).
8. `roi_shares.py` – synapses made and received by neck neurons, extracted from the partner table.
9. `connectivity_scaling.py` – caliber versus synaptic outputs and inputs by region, strong partners, cell-type population size and predicted neurotransmitter (`results/connectivity_scaling.json`, `results/channel_axons_male_connectivity.parquet`, `results/type_population_*.csv`).
10. `banc_neck_sections.py` – female mesh sections at the BANC neck plane (multiprocessing; ~30 min).
11. `compare_sex.py`, `compare_sex2.py` – male–female comparison, relative calibers.
12. `fig1_census.py`, `fig2_budget.py`, `fig3_allocation.py`, `fig4_connectivity.py`, `fig5_sex.py`, `figS1_validation.py`, `figS2_robustness.py` – figures (`figlib.py` holds shared styling; the optional panel-alignment gate is a no-op if the `audit_panel_alignment` module of the nature-figure toolkit is absent).

## Citation

Archived version 1.1.0: Zenodo, https://doi.org/10.5281/zenodo.22794997 (version 1.0, Report manuscript: https://doi.org/10.5281/zenodo.22758902). If you use these tables or code, please cite the Zenodo record and the paper (reference to be added on publication), together with the MaleCNS and BANC connectome papers listed in the manuscript.

## License

Code: MIT License. Data tables and figures: CC BY 4.0.
