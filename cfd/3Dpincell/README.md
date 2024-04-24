# 3D Pincell
This directory contains all of the cases necessary for modelling a full 3D pincell. Given the time constraints, this was the only case that was modelled, but it was initially planned to model an entire assembly.

## Methodology
This is described in more detail in [the final report](../../finalDeliverables/Final_Report.pdf), but to avoid needing a multiregion solver to accurately model the coolant salt flowing over the tube as well as the fuel salt flowing _inside_ the tube, the problem was split into _two_ calculations: one for the coolant salt with a fixed heat flux on the boundary, and one for the fuel salt with a fixed temperature on the boundary and a given volumetric heat generation term.