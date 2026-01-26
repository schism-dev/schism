# STOFS-3D-Setup Configurations

This folder holds **case-specific YAML configuration files** that describe how to set up and run STOFS-3D model experiments.

Each YAML file corresponds to a single setup recipe (see `src/stofs3d_setup/recipes/`), defining:
- which **recipe** to call (`profile:` key),
- model parameters (e.g., `rnday`, `startdate`, `version`),
- paths to **mesh**, **project**, and **scratch** directories,
- switches for which **input files** to generate.

### How to use
Run a case from the project root:
```bash
stofs3d-setup run --config configs/la_v8_2024_reforecast.yml
```
or inside Python (see test_recipe.py under the root folder):

```python
Copy code
from stofs3d_setup.config.schema import Settings
from stofs3d_setup.recipes.la_subdomain import build

cfg = Settings.from_yaml("configs/la_v8_2024_reforecast.yml")
build(cfg)
```

