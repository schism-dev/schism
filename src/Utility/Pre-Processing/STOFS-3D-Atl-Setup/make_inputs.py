# test_recipe_direct.py
from stofs3d_setup.config.schema import Settings
from stofs3d_setup.recipes.generic import build

cfg = Settings.from_yaml("/sciclone/data10/feiye/stofs3d-setup/configs/v7p4_2017.yml")
build(cfg)

pass
