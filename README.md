# GS-DeePC Project  
Gain-Scheduling Data-Enabled Predictive Control for Nonlinear Systems

---

## Publication Reference  

**Gain-Scheduling Data-Enabled Predictive Control for Nonlinear Systems with Linearized Operating Regions**  
Sebastian Zieglmeier, Mathias Hudoba de Badyn, Narada D. Warakagoda, Thomas R. Krogstad, Paal Engelstad  

This work presents a **Gain-Scheduled Data-Enabled Predictive Control (GS-DeePC)** framework for nonlinear systems based on multiple locally linear data representations. Instead of relying on a single global Hankel matrix, the operating range of a measurable scheduling variable is partitioned into regions, and regional Hankel matrices are constructed from persistently exciting data. To ensure smooth transitions between linearization regions and suppress region-induced chattering, **composite (overlapping) regions** are introduced by merging neighboring data sets, enabling a robust switching mechanism. The proposed approach preserves the original DeePC optimization structure and can reduce computational complexity by requiring only short, locally informative data sequences. Experimental results on a nonlinear DC motor with an unbalanced disc demonstrate significantly improved performance compared to standard DeePC.

---

## Overview  

This repository contains the MATLAB implementation used in the paper for **Gain-Scheduled DeePC (GS-DeePC)**. The method extends classical DeePC to nonlinear systems by exploiting locally valid linear input–output behavior across multiple operating regions, while remaining fully data-driven.

Key features:
- Gain-scheduled DeePC using multiple regional Hankel matrices  
- Region-based partitioning of a measurable scheduling variable  
- Composite (overlapping) regions for smooth scheduling and robust switching  
- Preservation of the standard DeePC optimization problem structure  
- Experimental validation on a nonlinear DC motor with an unbalanced disc  

---

## Requirements & Dependencies  

| Tool | Version | Notes |
|-----|--------|-------|
| MATLAB | R2024b (tested) | Core development environment |
| MOSEK | 10.x or newer | Convex optimization solver (academic license supported) |

---

## Citation  

If you use this code in your research, please cite:

    @article{Zieglmeier2025_GSDeePC,
      author  = {Sebastian Zieglmeier and Mathias Hudoba de Badyn and Narada D. Warakagoda and Thomas R. Krogstad and Paal Engelstad},
      title   = {Gain-Scheduling Data-Enabled Predictive Control for Nonlinear Systems with Linearized Operating Regions},
      journal = {arXiv},
      year    = {2025}
    }

---

## Contact  

**Sebastian Zieglmeier**  
Email: seb.zieglmeier@gmail.com  
Email: s.g.zieglmeier@its.uio.no  
GitHub: SebsDevLab  

**Mathias Hudoba de Badyn**  
Email: mathias.hudoba.de.badyn@its.uio.no  
