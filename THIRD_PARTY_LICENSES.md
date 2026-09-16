# Third-Party Software Acknowledgement


## CageCavityCalc (C3)


DSP-Cage incorporates and modifies several computational components
originally developed in CageCavityCalc.


Original repository:

https://github.com/VicenteMartiCentelles/CageCavityCalc


Original author:

Vicente Marti-Centelles


Copyright:

Copyright (c) 2021 VicenteMartiCentelles


License:

MIT License


The original license is preserved in:

LICENSE_C3


---

## Adapted Components

DSP-Cage does not reuse the complete implementation of CageCavityCalc.
Instead, selected low-level modules, data structures, and scientific
parameters from CageCavityCalc were incorporated to support molecular
structure preprocessing, atomic property representation, and molecular
property calculations.

The core methodology of DSP-Cage, including the Dynamic Splitting of Probes
(DSP) strategy, adaptive probe sampling mechanism, cavity exploration
algorithm, and GUI-oriented visualization workflow, was independently
developed in this work and is fundamentally different from the original
CageCavityCalc approach.

The relationship between DSP-Cage and CageCavityCalc is summarized below:


| Component | DSP-Cage File | Description |
|---|---|---|
| Molecular structure input and atom data organization | Cavity_Calculation.py, Input_method.py | Selected file reading functions and molecular atom naming/data organization structures were reused for initial molecular structure input and preprocessing. |
| Atomic parameter database | data.py | Atomic physical parameters, including van der Waals radii, atomic masses, and third-party molecular property parameters (hydrophValuesGhose1998 and hydrophValuesCrippen1999), were retained without modification. These parameters are used for distance calculation and molecular property evaluation. |
| Molecular hydrophobicity calculation | hydrophobicity.py | The fundamental principles and mathematical formulation of molecular hydrophobicity potential (MHP) calculation were referenced. The sampling strategy and spatial mapping procedure were redesigned to fit the DSP-Cage framework and GUI visualization requirements. |
| Electrostatic property calculation | electrostatics.py | The theoretical basis and calculation principles of electrostatic potential (ESP) estimation were referenced. The implementation was modified and integrated with the DSP-Cage dynamic probe sampling strategy. |



---

## Modification and Extension


The adapted components were modified and integrated into DSP-Cage.


DSP-Cage additionally introduces:


- Dynamic Splitting of Probes strategy;
- adaptive cavity exploration;
- guest-aware cavity threshold determination;
- dynamic cavity tracking.


These components were newly developed for DSP-Cage.


---
## Differences from CageCavityCalc

Although DSP-Cage incorporates several low-level components from
CageCavityCalc, the two frameworks use fundamentally different cavity
characterization strategies.


CageCavityCalc mainly performs cavity identification through a static
grid-based sampling strategy combined with geometric angle-based criteria.
The cavity region is determined by predefined grid points and spatial
geometric evaluation.


In contrast, DSP-Cage introduces a Dynamic Splitting of Probes strategy,
where cavity characterization is performed through adaptive probe generation,
dynamic spatial exploration, and probe-based cavity reconstruction.

The major methodological differences include:


1. Sampling strategy

CageCavityCalc:
- Uses static grid points;
- Determines cavity points through grid occupancy and geometric angle
  evaluation.

DSP-Cage:
- Uses dynamically generated probes;
- Performs adaptive probe splitting and spatial exploration;
- Does not rely on the original grid-angle sampling strategy.


2. Cavity reconstruction mechanism

CageCavityCalc:
- Reconstructs cavity regions based on predefined spatial grids.

DSP-Cage:
- Reconstructs cavity regions through dynamic probe evolution and adaptive
  sampling.


3. Visualization and interaction

DSP-Cage additionally introduces GUI-oriented spatial mapping and
visualization mechanisms, which were not included in the original
CageCavityCalc implementation.


Therefore, the reused components mainly serve as basic computational
infrastructure, while the central algorithmic contributions of DSP-Cage were
newly developed.
---

## Citation


Users should cite both:


DSP-Cage:

" DSP-Cage: Supramolecular Cage Cavity Characterization and Visualization via Dynamic Splitting of Probes "


and


CageCavityCalc:

DOI:
10.1021/acs.jcim.4c00355


---

## License Notice


The original CageCavityCalc components remain under the MIT License.

All other newly developed DSP-Cage components are released for academic and non-commercial research use.

Commercial use requires permission from the authors.
