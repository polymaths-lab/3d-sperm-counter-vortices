# Repository for Manuscript at Nature Communications

This repository contains all codes used in the manuscript currently under consideration at *Nature Communications*.  
An earlier version of the manuscript is available on bioRxiv:

https://www.biorxiv.org/content/10.1101/2024.05.29.596379v1.full

> Note: The bioRxiv version is not the latest version of the manuscript; the most recent version corresponds to the submission to *Nature Communications*.

---

## Overview of the Repository Structure

The repository is organised into four main types of folders:

- `Data*`
- `MT_Fig*`
- `SI_Fig*`
- `SI_Movie*`

These correspond to different stages of data processing and figure/movie generation.

---

## 1. `Data*` — Experimental Data Processing, Numerical Reconstruction, and Idealised Model Generation

Example folders:
- `Data0_Experimental Sperm`
- `Data1_Reconstruction`

These folders contain codes for:

- Pre-processing of experimental sperm cell data (e.g. smoothing spline).
- Numerical reconstruction of experimentally observed sperm motility.
- Computation of the induced flow fields.
- Generation of idealised sperm models, including:
  - waveform,
  - motility,
  - flow fields.

### Important notes on data availability

- The **raw experimental source data are provided** in this repository.
- However, **calculated data are not provided**, including:
  - reconstructed motility data,
  - flow field data,
  - idealised model waveform data.

This is because:
1. These data can be reproduced by readers using the provided source data and codes.
2. The total size of these calculated datasets is very large.
3. These datasets are intermediate results and are intended to be generated locally by users.

⚠️ **Crucial**:  
Although not provided, these calculated datasets are **essential inputs** for running the codes in the other three types of folders (`MT_Fig*`, `SI_Fig*`, `SI_Movie*`).

---

## 2. `MT_Fig*` — Main Text Figures

Example folders:
- `MT_Fig1_Aim Fig`
- `MT_Fig2_Recon Sperm Motility`

These folders contain all codes used to generate the **main text figures**, corresponding to:

> Figures 1–7 in the main manuscript.

They rely on the processed and reconstructed data generated in the `Data*` folders.

---

## 3. `SI_Fig*` — Supplementary Figures

Example folders:
- `SI_Fig1_Recon Sperm Motility`
- `SI_Fig2_ReconFlow AverageNear CrossSection`

These folders contain all codes used to generate the **Supplementary Figures**, corresponding to:

> Supplementary Figures 1–8.

They also require the data generated from the `Data*` folders.

---

## 4. `SI_Movie*` — Supplementary Videos

Example folders:
- `SI_Movie1234_Traj WF`
- `SI_Movie5_ReconNearFlow`

These folders contain all codes used to generate the **Supplementary Movies**, corresponding to:

> Supplementary Movies 1–8.

Again, these rely on the motility and flow data produced in the `Data*` folders.

---

## Reproducibility

In summary, the full computational pipeline is:

1. Start from experimental source data in `Data*`.
2. Run reconstruction and modelling codes in `Data*` to generate intermediate datasets.
3. Use these datasets as inputs for:
   - main text figures (`MT_Fig*`),
   - supplementary figures (`SI_Fig*`),
   - supplementary movies (`SI_Movie*`).

This design allows full reproducibility of all results in the manuscript while keeping the repository size manageable.
