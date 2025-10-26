# ES_Project
**Module-Selection Balance in the Evolution of Modular Organisms**  
Minkyu Kim, Sarah M. Ardell, Sergey Kryazhimskiy (2025)

---

<!-- =============================
SECTION: OVERVIEW
You can keep this as-is except to update the final publication details.
================================= -->
### Overview
This repository contains all MATLAB R2025a code used in  
Kim, Ardell & Kryazhimskiy (2025), *Module-Selection Balance in the Evolution of Modular Organisms*.  
It includes implementations of both **Pleiotropic** and **Modular** Genotype–Phenotype–Fitness Models (GPFMs) and all figure-generation scripts used in the paper.

All analyses, figures, and simulations were performed using MATLAB R2025a.  
Repository maintained by Minkyu Kim (Cornell University).

---

<!-- =============================
SECTION: DIRECTORY STRUCTURE
Update this if you reorganize or rename folders.
================================= -->
### Directory Structure

| Folder | Description |
| ------- | ----------- |
| `code/ModularGPFM/` | Modular GPFM simulations and figure scripts (Figs 3–6) |
| `code/PleiotropicGPFM/` | Pleiotropic GPFM simulations (Fig 2) |
| `code/utils/` | Shared plotting, seeding, and validation utilities |
| `code/run_examples/` | Minimal working demos to test environment |
| `data/` | Example `.mat` files and initial phenotype data |
| `Figures/` | Generated EPS/PDF figures |
| `Supplementary/` | Supplementary information document |
| `make_all_figures.m` | Master script to reproduce all figures |

---

<!-- =============================
SECTION: REQUIREMENTS
Modify if toolboxes change or you add dependencies.
================================= -->
### Requirements

- **MATLAB R2025a**
- Required Toolboxes:
  - Statistics and Machine Learning Toolbox  
  - Parallel Computing Toolbox (optional, for batch simulations)
- Tested on:
  - macOS 14.5 (Apple Silicon)
  - Windows 11
  - Linux (RHEL 9)

---

<!-- =============================
SECTION: QUICK START
You can leave the steps mostly as-is.
================================= -->
### Quick Start

1. **Clone or download the repository**
   ```bash
   git clone https://github.com/mkkim1894/ES_Project.git
   cd ES_Project
   ``` 
2. **Open MATLAB and add paths**
   ```matlab
   addpath(genpath('code'))
   ```
3. **Run demo**
   ```matlab
   run_examples/demo_modular
   ```
4. **Reproduce all figures**
   ```matlab
   make_all_figures
   ```

---

### Citation

If you use this code, please cite:

> Kim, M., Ardell, S. M., & Kryazhimskiy, S. (2025).  
> *Module-Selection Balance in the Evolution of Modular Organisms.*  
> [Journal Name, Volume, Pages].  
> DOI: [replace_with_final_DOI]

Machine-readable metadata is included in the file `CITATION.cff`.

---

### License

This project is released under the **MIT License**.  
See the [LICENSE](./LICENSE) file for full text.

---

### Contact

**Minkyu Kim**  
Department of Computational Biology, Cornell University  
Email: `mk2687@cornell.edu`  
GitHub: [mkkim1894](https://github.com/mkkim1894)

---

### Repository DOI

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.xxxxxxx.svg)](https://doi.org/10.5281/zenodo.xxxxxxx)
