# VectorFit-related Analysis

This repository contains adapted versions of the VectorFit scripts used for the characterization of optical aberrations across different image-splitter configurations.

## Attribution & Licensing

The core framework in this repository is adapted from **VectorFit**, developed by the **ImPhys group at TU Delft**. 
- **Original Source:** [TU Delft GitLab (VectorFit)](https://gitlab.tudelft.nl/imphys/ci/vectorfit)
- **License:** Distributed under the [Apache-2.0 license](https://www.apache.org/licenses/LICENSE-2.0).

### Adaptations in this Repository:
This version includes modified implementations and results based on:
- `fit_beads_aberrations.m`
- `get_NAT_inversion.m` (Inversion of **Nodal Aberration Theory** models)
- `segmentation_beads_difference_uniform_cropped.m`
- Several `set_parameters.m` variants optimized for hardware-specific configurations.

---

## Academic Context

The adaptations and results hosted in this repository are part of the following Master's Thesis research:

**Title:** *AUTOMATED COMPUTATIONAL CHARACTERIZATION OF ABERRATIONS IN FLUORESCENCE MICROSCOPY: A COMPARATIVE ANALYSIS*  
**Author:** Aleksandra Will  
**Institution:** IMC Hochschule für Angewandte Wissenschaften Krems (University of Applied Sciences)  
**Submitted:** 15.03.2026

### Research Focus
The scripts in this repository facilitate a comparative analysis of field-dependent optical aberrations. The research focused on characterizing how different image-splitter configurations impact PSF quality across the Field of View (FOV):

- **NAT Modeling**: **Nodal Aberration Theory (NAT)** was used as a predictive model to define the theoretical behavior and shape of aberrations across the field.
- **Field Mapping**: Experimental bead measurements were processed using **cubic interpolation** to generate continuous maps of aberration magnitudes across the FOV.
- **Characterization**: By comparing empirical PSF measurements with NAT predictions, the study provides a detailed characterization of the optical performance and alignment of various image-splitter hardware setups.
