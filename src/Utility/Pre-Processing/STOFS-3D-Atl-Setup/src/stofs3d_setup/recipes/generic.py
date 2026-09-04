"""
STOFS3D Atlantic model recipe for Atlas hindcasts
"""
import shutil
from pathlib import Path

from ..config.schema import Settings
from ..core.stofs3d_atl_driver import stofs3d_atl_driver
from ..config.stofs3d_atl_config import ConfigStofs3dAtlantic

def get_model_cfg(version: str):
    return getattr(ConfigStofs3dAtlantic, version)()

def build(cfg: Settings):
    base_cfg = getattr(ConfigStofs3dAtlantic, cfg.model.version)()
    model_cfg = base_cfg.model_copy(update=cfg.model.model_dump(exclude_unset=True))
    if cfg.source_yaml_path is not None:
        input_dir = Path(cfg.run.project_dir) / f"I{cfg.run.runid}"
        input_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(cfg.source_yaml_path, input_dir / cfg.source_yaml_path.name)

    return stofs3d_atl_driver(
        hgrid_path=cfg.run.hgrid_path,
        vgrid_path=cfg.run.vgrid_path,
        config=model_cfg,
        project_dir=cfg.run.project_dir,
        runid=cfg.run.runid,
        scr_dir=(cfg.runtime.scr_dir if cfg.runtime else None),
        input_files=cfg.inputs.model_dump(by_alias=True),
    )
