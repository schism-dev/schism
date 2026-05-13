# test_recipe_direct.py
from stofs3d_setup.config.schema import Settings
from stofs3d_setup.recipes.generic import build

cfg = Settings.from_yaml("/sciclone/home/feiye/schism/src/Utility/Pre-Processing/STOFS-3D-Atl-Setup/configs/v7p4_sample.yml")
build(cfg)

pass
