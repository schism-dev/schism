# HydrologicExchange development rules

- Follow the established SCHISM Fortran style, especially the hydrodynamic core in `src/Hydro/schism_step.F90` and `src/Hydro/schism_init.F90`.
- Match nearby naming, declarations, indentation, comments, module organization, and build conventions. Avoid unrelated reformatting.
- Keep code and interfaces lean, direct, and human readable. Prefer the smallest structure that supports the agreed behavior.
- Do not add excessive validation, defensive checks, generic abstractions, or adapter infrastructure before the component structure and data contracts are decided. Add checks where they are needed for correctness or an agreed interface.
- Keep process calculations separate from code that applies fluxes to SCHISM hydrodynamics.
