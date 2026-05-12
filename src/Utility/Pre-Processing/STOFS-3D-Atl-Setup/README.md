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

#### Why this matters
- Keeps **version logic centralized**
- Prevents silent divergence between runs
- Enables clear comparison across versions (e.g., v7.3 vs v7.4)

---

### 2. Define Runs via YAML Configuration

All run-specific settings must be defined in a YAML file rather than hard-coded or scripted ad hoc.

Example: `v7p4_2017.yml`

```yaml
profile: "v7.4_2017"

run:
  project_dir: "/sciclone/schism10/feiye/STOFS3D-v7.4/"
  runid: "R01"
  hgrid_path: "/sciclone/schism10/Hgrid_projects/STOFS3D-v7.4/I01/hgrid.gr3"

model:
  version: "v7p4"
  rnday: 396
  startdate: "2016-12-01"
  replace_nwm_with_usgs: True

runtime:
  scr_dir: "/sciclone/scr10/feiye/STOFS3D-v7.4/"

inputs:
  bctides: true
  vgrid: false
  gr3: false
  nudge_gr3: false
  shapiro: false
  drag: false
  elev_ic: false
  flux_th: true
  source_sink: true
  soil: false
  "hotstart.nc": true
  "3D.th.nc": true
  "*nu.nc": true
  "elev2D.th.nc": true
  "*.prop": false
```

#### Why this matters
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

cfg = Settings.from_yaml("/sciclone/data10/feiye/stofs3d-setup/configs/v7p4_2017.yml")
build(cfg)
```

#### Why this matters
- Ensures **all preprocessing steps are executed consistently**
- Avoids missing or inconsistent inputs
- Keeps logic aligned with repository updates

---

### ❗ Important Requirement

**Do NOT:**
- Copy and modify individual scripts to generate inputs independently  
- Mix configurations from different versions manually  
- Bypass the YAML + factory method + recipe workflow  

---

### Summary

```
Factory Method (version definition)
        ↓
YAML Config (run definition)
        ↓
Recipe Build (execution)
```
