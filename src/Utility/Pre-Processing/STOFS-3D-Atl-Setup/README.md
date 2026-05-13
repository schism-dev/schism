# stofs3d-setup

Scripts for setting up **STOFS-3D Atlantic** and related compound flooding modeling applications based on SCHISM.

These scripts automate common preprocessing tasks such as mesh preparation, forcing setup, and workflow configuration.

This package replaces the legacy STOFS-3D-Atl-shadow-VIMS package. Previous versions of the original implementation can be retrieved from commit 71dd94b9 and earlier Git history.

---

## Dependencies

The scripts are tested on Python 3.10. 

---

## Installation

### 1. Install `pyschism` (development version)

Install the latest development version from GitHub:

```bash
git clone https://github.com/schism-dev/pyschism.git
cd pyschism
pip install -e . --no-build-isolation
```

The `--no-build-isolation` flag is currently required due to packaging issues in the upstream repository.

---


### 2. Install experimental pylibs utilities

```bash
pip install git+https://github.com/wzhengui/pylibs.git#subdirectory=pylib_experimental
```

---

### 3. Temporary dependency

This dependency will be removed in a future update:

```bash
pip install git+https://github.com/feiye-vims/schism_py_pre_post.git
```

---

### 4. This package:
pip install "git+https://github.com/schism-dev/schism.git@master#subdirectory=src/Utility/Pre-Processing/STOFS-3D-Atl-Setup"

---

## Notes

- This repository assumes familiarity with the **SCHISM/STOFS modeling workflow**.
- Several dependencies are currently transitional and will be simplified as upstream packages stabilize.


## Recommended Usage and Configuration Workflow

To ensure consistency, reproducibility, and maintainability across STOFS-3D setups, **users are strongly required to follow the standardized workflow provided in this repository**.  

While it may be tempting to reuse isolated scripts or extract portions of the code to generate individual inputs, **this practice is discouraged**. Doing so breaks traceability between configurations, versions, and generated artifacts, making it difficult to debug issues, compare experiments, or reproduce results.

> **Key Principle:**  
> Even if this workflow is not always the most flexible or fastest for one-off tasks, it establishes a consistent structure that allows us to track model evolution, configuration changes, and experiment history in a reliable way.

---

### 1. Define a Version via Factory Method

Each new STOFS-3D version should be defined explicitly in  
`stofs3d_setup/config/stofs3d_atl_config.py` using a factory method.

This ensures that **all version-specific assumptions are centralized and documented**.

Example:

```python
@classmethod
def v7p4(cls):
    '''Factory method to create a configuration for STOFS3D-v7.4 3D setup'''
    return cls(
        ocean_bnd_ids=[0, 1],
        elev2d_uniform_shift=-0.42,
        nudging_zone_width=7.3,
        shapiro_zone_width=11.5,
        shapiro_tilt=3.5,
        feeder_info_file=None,
        hgrid_without_feeders=None,
        relocate_source=True,
        mandatory_sources_coor=rsf.v19p2_for_sms_v27_mandatory_sources_coor,
        nwm_cache_folder=None,
        source_ele_replace_dict={
            # 53: 3552194,
        },
        bc_flags=[
            [5, 5, 4, 4],
            [5, 5, 4, 4],
            [0, 1, 1, 2],
        ],
        bc_const=[
            [None, None, None, None],
            [None, None, None, None],
            [None, None, None, 0.0],
        ],
        bc_relax=[
            [None, None, 0.5, 0.5],
            [None, None, 0.5, 0.5],
            [None, None, 0.01, 1.0],
        ],
        tvd_regions=None
    )
```

#### Rationale
- Keeps **version logic centralized**
- Prevents silent divergence between runs
- Enables clear comparison across versions (e.g., v7.3 vs v7.4)

---

### 2. Define Runs via YAML Configuration

All run-specific settings must be defined in a YAML file rather than hard-coded or scripted ad hoc.

Example: `v7p4_sample.yml`

```yaml
profile: "v7.4_sample"  # give any name

run:
  project_dir: "/sciclone/schism10/feiye/STOFS3D-v7.4/"
  runid: "100x"
  hgrid_path: "/sciclone/schism10/hjyoo/task/task10_Atlantic/RUN100d/src/hgrid/hgrid.gr3"

model:
  version: "v7p4"  # defined in src/stofs3d_setup/config/stofs3d_atl_config.py
  rnday: 3
  startdate: "2016-12-01"  # *
  replace_nwm_with_usgs: True
  nwm_cache_folder: None  # Use None or delete this line if you don't have saved NWM product

runtime:
  scr_dir: "/sciclone/scr10/feiye/STOFS3D-v7.4/"  # Use None if you don't want to save the outputs on a scratch disk.

inputs:
  bctides: false
  vgrid: false
  gr3: true
  nudge_gr3: false
  shapiro: true
  diffmin: true
  drag: true
  elev_ic: false
  flux_th: false
  source_sink: false
  soil: false
  "hotstart.nc": false
  "3D.th.nc": false
  "*nu.nc": false
  "elev2D.th.nc": false
  "*.prop": false
```

#### Rationale
- Provides a **single source of truth** for each run
- Makes experiments **portable and reproducible**
- Enables systematic comparison across scenarios

---

### 3. Build Using the Standard Recipe Interface

All setups must be generated through the provided recipe system.

Example:

```python
from stofs3d_setup.config.schema import Settings
from stofs3d_setup.recipes.generic import build 

cfg = Settings.from_yaml("/sciclone/data10/feiye/stofs3d-setup/configs/v7p4_sample.yml")
build(cfg)
```

This will generate three folders under the project directory, for example:

```bash
/sciclone/schism10/feiye/STOFS3D-v7.4/I01/
```

The `I` folder (`Inputs`) contains:

- a copy of all scripts
- the generated model input files
- a metadata.yml file containing version information (e.g., git commit hash and tag), which can be used to retrieve the exact script version from git

This folder serves as a reproducible snapshot of the input-generation workflow.

```bash
/sciclone/schism10/feiye/STOFS3D-v7.4/R01/
```

The `R` folder (`Run`) is where the model simulation is executed.  
It mainly contains symbolic links to the input files stored in the corresponding `I` folder.

```bash
/sciclone/schism10/feiye/STOFS3D-v7.4/O01/
```

The `O` folder (`Outputs`) is intended for post-processing products and analysis outputs.  
You may ignore this folder if you prefer to store post-processed results elsewhere.

#### Rationale
- Ensures **all preprocessing steps are executed consistently**
- Avoids missing or inconsistent inputs
- Keeps logic aligned with repository updates

---

### Summary

```
Factory Method (version definition)
        ↓
YAML Config (run definition)
        ↓
Recipe Build (execution)
```

### Notes

- For STOFS version `v7.4` and later, use the latest `master` branch.

- For STOFS version `v7.3` and earlier, check out the corresponding historical setup scripts:
```bash
git checkout d769b3c66d6
```