# src/stofs3d_setup/config/schema.py
from pydantic import BaseModel, FilePath, DirectoryPath, Field, ConfigDict
from datetime import datetime
import yaml
from typing import Optional

class Run(BaseModel):
    project_dir: DirectoryPath
    runid: str
    hgrid_path: FilePath
    vgrid_path: Optional[FilePath] = None

class ModelCfg(BaseModel):
    model_config = ConfigDict(extra="allow", arbitrary_types_allowed=True)

    version: str
    rnday: int = Field(gt=0)
    startdate: datetime

class Runtime(BaseModel):
    scr_dir: Optional[DirectoryPath] = None

class Inputs(BaseModel):
    bctides: bool = False
    vgrid: bool = False
    gr3: bool = False
    nudge_gr3: bool = False
    shapiro: bool = False
    drag: bool = False
    flux_th: bool = False
    elev_ic: bool = False
    soil: bool = False
    source_sink: bool = False

    hotstart_nc: bool | None = Field(None, alias="hotstart.nc")
    th3d_nc: bool | None = Field(None, alias="3D.th.nc")
    elev2d_nc: bool | None = Field(None, alias="elev2D.th.nc")
    nu_glob: bool | None = Field(None, alias="*nu.nc")
    prop_glob: bool | None = Field(None, alias="*.prop")

class Settings(BaseModel):
    profile: str
    run: Run
    model: ModelCfg
    runtime: Runtime | None = None
    inputs: Inputs

    @classmethod
    def from_yaml(cls, p):
        with open(p, "r") as f:
            data = yaml.safe_load(f)
        return cls.model_validate(data)
