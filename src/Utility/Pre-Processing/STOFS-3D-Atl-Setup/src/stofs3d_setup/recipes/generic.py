"""
STOFS3D Atlantic model recipe for Atlas hindcasts
"""
from ..config.schema import Settings
from ..core.stofs3d_atl_driver import stofs3d_atl_driver
from ..config.stofs3d_atl_config import ConfigStofs3dAtlantic

def get_model_cfg(version: str):
    return getattr(ConfigStofs3dAtlantic, version)()

def build(cfg: Settings):
    base_cfg = getattr(ConfigStofs3dAtlantic, cfg.model.version)()
    model_cfg = base_cfg.model_copy(update=cfg.model.model_dump(exclude_unset=True))

    return stofs3d_atl_driver(
        hgrid_path=cfg.run.hgrid_path,
        vgrid_path=cfg.run.vgrid_path,
        config=model_cfg,
        project_dir=cfg.run.project_dir,
        runid=cfg.run.runid,
        scr_dir=(cfg.runtime.scr_dir if cfg.runtime else None),
        input_files=cfg.inputs.model_dump(by_alias=True),
    )
