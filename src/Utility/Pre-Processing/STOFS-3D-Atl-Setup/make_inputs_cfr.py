# test_recipe_direct.py
from stofs3d_setup.config.schema import Settings
from stofs3d_setup.recipes.generic import build 

cfg = Settings.from_yaml("/home/liuquncw/stofs3d-setup/configs/cfr_v1_hindcast.yml")
build(cfg)

pass
