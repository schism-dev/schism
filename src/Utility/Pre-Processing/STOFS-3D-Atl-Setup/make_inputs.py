# test_recipe_direct.py
from stofs3d_setup.config.schema import Settings
from stofs3d_setup.recipes.generic import build

cfg = Settings.from_yaml("/sciclone/data10/feiye/stofs3d-setup/configs/la_v8_2023_hindcast_R33.yml")
build(cfg)

pass
