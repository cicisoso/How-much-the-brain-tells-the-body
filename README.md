# Scaling laws of axon calibre and the capacity of the brain–body channel in *Drosophila*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22758901.svg)](https://doi.org/10.5281/zenodo.22758901)

Code, per-axon tables, figure source data and figures for the manuscript *"Scaling laws of axon calibre and the capacity of the brain–body channel in Drosophila"* (Liu, Zong, Chen and Xiong; School of Information Technology, Zhejiang Financial College; under submission to *Journal of the Royal Society Interface*). Earlier releases accompanied earlier versions of the manuscript (see Versions).

We measured the cross-sectional area of every brain–body axon in the neck connective of the male *Drosophila* central nervous system connectome (MaleCNS v1.0) directly from the 16 nm segmentation. The measurement was validated four ways: 8 nm re-measurement, bilateral homologues, mesh sections and the giant fibre. We then related each axon's calibre to its synaptic outputs and inputs, its arbor and release sites, its cell-type population size and its predicted neurotransmitter, and repeated the calibre census on the female BANC connectome by sectioning neuron meshes at the curated neck plane. Version 1.2.0 adds the arbor analysis: synaptic output scales as *d*^1.4 and decomposes exactly into target-region cable (∝ *d*^0.7–0.8), release sites per µm (∝ *d*^0.6) and contacts per release site (constant). It also adds the descending/ascending capacity ratio as a function of the rate–calibre exponent in both sexes.

## Contents

| Folder | What it holds |
|---|---|
| `code/` | Analysis scripts (Python 3.11) and the figure scripts of version 1.1.0. See the pipeline below. |
| `code/rsif/` | Figure scripts for the current manuscript (figures 1–6 and S1–S2, Royal Society style; `figlib_rsif.py` holds the shared style). |
| `data_tables/` | Current manuscript: Data S1 (`Data_S1_per_axon_tables.xlsx`: README, MaleCNS and BANC sheets, with CSV copies `DataS1_*.csv`) and Data S2 (`Data_S2_figure_source_data.xlsx`: one sheet per figure panel and table 1). Version 1.1.0: `Supplementary_Data_1.xlsx`, `Source_Data.xlsx` and their CSV copies. |
| `results/` | Intermediate results used by the figure scripts (JSON summaries, CSV and parquet tables), including `arbor_scaling.json`, `channel_axons_male_arbor.parquet`, `capacity_sexes.json` and `connectivity_scaling.json`. |
| `figures/rsif/` | Figures 1–6 and S1–S2 of the current manuscript (PNG and PDF; the scripts also write 600 dpi TIFF files, which are not tracked). |
| `figures/` | Figures 1–5 and S1–S2 of version 1.1.0 (PNG and PDF). |

## Data sources (not included; public)

- MaleCNS v1.0: https://male-cns.janelia.org (CC-BY 4.0). Files used: `body-annotations-male-cns-v1.0-minconf-0.5.feather`, `body-stats-male-cns-v1.0-minconf-0.5.feather`, `body-neurotransmitters-male-cns-v1.0.feather`, `syn-partners-male-cns-v1.0-minconf-0.5.feather`, SWC skeletons (`gs://flyem-male-cns/v1.0/segmentation/skeletons-malecns/skeletons-swc/`), and the segmentation (`gs://flyem-male-cns/v1.0/segmentation`, read with cloud-volume).
- BANC v888: public bucket `gs://lee-lab_brain-and-nerve-cord-fly-connectome` (files `compiled_data/banc_888/banc_888_meta.feather`, `neuron_annotations/v888/neck_connective_y92500.parquet`, and the `neuron_meshes` precomputed mesh source).

## Setup

```bash
pip install -r requirements.txt
export FLYBRAIN_DATA=/path/to/data   # folder holding the MaleCNS feather files, swc/ skeletons, neck_slabs/ and banc/ files
```

Scripts resolve the data folder from `FLYBRAIN_DATA` (default `../data`) and read and write `results/`, `figures/` and `data_tables/` relative to the working directory. Run them from the repository root, for example `python3 code/arbor_scaling.py`.

## Pipeline

1. `neck_waist.py`: locate the connective along the body axis from skeleton spread.
2. `male_neck_slabs.py`: download three 16 nm slabs of the segmentation and count voxels per segment and plane.
3. `caliber_analysis.py`: per-axon areas, tilt correction and the per-body calibre table.
4. `male_neck_slab_mip0.py` and `mip_validation.py`: 8 nm validation slab.
5. `male_mesh_validation.py`: mesh-section versus voxel validation on 175 axons.
6. `delay_budget.py`: skeleton path lengths from the neck to the brain and VNC arbors; conduction-delay model.
7. `budget_male.py`: information-budget tables and class summaries (`results/budget_male.json`).
8. `roi_shares.py`: synapses made and received by neck neurons, extracted from the partner table (writes `neck_syn_pre.parquet` and `neck_syn_post.parquet` to the data folder).
9. `connectivity_scaling.py`: calibre versus synaptic outputs and inputs by region, strong partners, cell-type population size and predicted neurotransmitter (`results/connectivity_scaling.json`, `results/channel_axons_male_connectivity.parquet`, `results/type_population_*.csv`).
10. `banc_neck_sections.py`: female mesh sections at the BANC neck plane (multiprocessing; about 30 min).
11. `compare_sex.py` and `compare_sex2.py`: male–female comparison and relative calibres.
12. `arbor_scaling.py` (new in 1.2.0): for every channel axon, measures the cable, terminal tips, branch points and mean skeleton calibre in the brain and the VNC, and the release sites (distinct T-bars) in each region. It then computes:
    - the exact decomposition of the output exponent into cable, release-site density and contacts per site;
    - the partial exponent of diameter at fixed target cable, whole-neuron cable, tip count, or target class plus cable;
    - the sensory-afferent controls and the tip-count exponents.

    Outputs: `results/arbor_scaling.json`, `results/arbor_scaling.parquet` and `results/channel_axons_male_arbor.parquet`. The skeleton step is cached; set `REBUILD=1` to recompute it.
13. `capacity_sexes.py` (new in 1.2.0): descending/ascending capacity ratio Σ*d*^β as a function of β, and the parity exponent, in both connectomes, with bootstrap 95% CIs from 2000 resamples (`results/capacity_sexes.json`).
14. `code/rsif/fig*.py` (new in 1.2.0): figures 1–6 and S1–S2 of the current manuscript, written to `figures/rsif/`.
15. `make_data_tables_rsif.py` (new in 1.2.0): Data S1 and Data S2, written to `data_tables/`.

The version 1.1.0 figure scripts (`fig1_census.py`, `fig2_budget.py`, `fig3_allocation.py`, `fig4_connectivity.py`, `fig5_sex.py`, `figS1_validation.py`, `figS2_robustness.py`, with `figlib.py`) are kept for provenance and write the earlier figure set to `figures/`. The optional panel-alignment gate in both figure libraries is a no-op if the `audit_panel_alignment` module of the nature-figure toolkit is absent; for `code/rsif/`, set `NATURE_FIGURE_SCRIPTS` to the toolkit's scripts folder to enable it. Version 1.2.0 also fixes an indentation error in `code/figlib.py` that stopped the version 1.1.0 figure scripts from running.

## Versions

| Version | Date | Manuscript | Zenodo DOI |
|---|---|---|---|
| 1.2.0 | 2026-10-02 | *Scaling laws of axon calibre and the capacity of the brain–body channel in Drosophila* (J. R. Soc. Interface) | 10.5281/zenodo.23097576 |
| 1.1.0 | 2026-09-16 | *The structural bandwidth of the neck connective between the Drosophila brain and body* (Communications Biology) | 10.5281/zenodo.22794997 |
| 1.0.0 | 2026-09-14 | *How much the brain tells the body* (Report) | 10.5281/zenodo.22758902 |

## Citation

Please cite the Zenodo record of the version you used. Version 1.2.0, which accompanies the current manuscript, is https://doi.org/10.5281/zenodo.23097576; the concept DOI https://doi.org/10.5281/zenodo.22758901 always resolves to the latest version. Please also cite the paper (reference to be added on publication), together with the MaleCNS and BANC connectome papers listed in the manuscript.

## License

Code: MIT License. Data tables and figures: CC BY 4.0.
