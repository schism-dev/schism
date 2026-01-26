"""
Configuration for STOFS-3D-ATL model,
used by stofs3d_atl_driver.py
"""

from __future__ import annotations
from pydantic import BaseModel, Field, ConfigDict
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Dict, Any
import numpy as np

from ..ops.Source_sink.Relocate import relocate_source_feeder as rsf


# ---------------------------------------------------------------------
#                               Classes
# ---------------------------------------------------------------------

class ConfigStofs3dAtlantic(BaseModel):
    """
    Pydantic-based configuration for STOFS-3D-Atlantic model.
    """

    model_config = ConfigDict(extra="ignore", arbitrary_types_allowed=True)  # ignore unknown keys gracefully

    # === Core simulation parameters ===
    startdate: datetime = Field(default_factory=lambda: datetime(2017, 12, 1))
    rnday: int = 60

    # === Boundary and forcing parameters ===
    ocean_bnd_ids: List[int] = Field(default_factory=lambda: [0])
    elev2d_uniform_shift: float = 0.0
    nudging_zone_width: float = 1.5
    nudging_day: float = 1.0
    shapiro_zone_width: float = 2.5
    shapiro_tilt: float = 2.0

    # === Boundary condition lists ===
    bc_flags: List[List[Optional[int]]] = Field(default_factory=list)
    '''
    Example:
    bc_flags = [
        [3, 3, 0, 0],  # tides for elev and vel
        [3, 3, 0, 0],  # tides for elev and vel
        [3, 3, 0, 0],  # tides for elev and vel
        [3, 3, 0, 0],  # tides for elev and vel
        [3, 3, 0, 0],  # tides for elev and vel
    ]
    '''
    bc_const: List[List[Optional[float]]] = Field(default_factory=list)
    '''
    Example:
    bc_const = [
        [None, None, None, None],
        [None, None, None, None],
        [None, None, None, None],
        [None, None, None, None],
        [None, None, None, None],
    ]
    '''
    bc_relax: List[List[Optional[float]]] = Field(default_factory=list)
    '''
    Example:
    self.bc_relax = [
        [None, None, None, None],
        [None, None, None, None],
        [None, None, None, None],
        [None, None, None, None],
        [None, None, None, None],
    ]
    '''

    # === Source/sink and feeder settings ===
    nwm_cache_folder: Optional[Path] = None
    relocate_source: bool = True
    feeder_info_file: Optional[Path] = None
    hgrid_without_feeders: Optional[Path] = None
    mandatory_sources_coor: Optional[np.ndarray] = None
    existing_source_json_path: Optional[Path] = None
    reuse_source_json: bool = False
    replace_nwm_with_usgs: bool = False
    source_ele_replace_dict: Dict[int, int] = None  # temporary fix for isolated feeder channels

    # === Miscellaneous values ===
    gr3_values: Dict[str, float] = Field(
        default_factory=lambda: {
            "albedo": 0.1,
            "diffmax": 1.0,
            "diffmin": 1e-6,
            "watertype": 1.0,
            "windrot_geo2proj": 0.0,
        }
    )
    tvd_regions: List[str] = Field(default_factory=list)

    # === Factory methods for different configurations ===
    @classmethod
    def v6(cls):
        '''Factory method to create a configuration for STOFS3D-v6'''
        return cls(
            nudging_zone_width=.3,  # very wide nudging zone
            shapiro_zone_width=11.5,  # very wide shapiro zone
            shapiro_tilt=3.5,  # very abrupt transition in the shapiro zone
            feeder_info_file=(
                '/sciclone/schism10/feiye/STOFS3D-v5/Inputs/v14/Parallel/'
                'SMS_proj/feeder/feeder.pkl'),
            mandatory_sources_coor=rsf.v19p2_mandatory_sources_coor,
            ocean_bnd_ids=[0, 1],
            bc_flags=[[5, 5, 4, 4]],
            tvd_regions=[
                'tvd0_1.reg', 'tvd0_2.reg', 'tvd0_3.reg', 'tvd0_4.reg',
                'tvd0_5.reg', 'tvd0_6.reg', 'tvd0_7.reg',
                'upwind_east_Carribbean.rgn', 'upwind_west_Carribbean.rgn'
            ]
        )

    @classmethod
    def v7(cls):
        '''Factory method to create a configuration for STOFS3D-v7'''
        return cls(
            ocean_bnd_ids=[0, 1],
            elev2d_uniform_shift=-0.42,  # add a uniform shift to elev2D
            nudging_zone_width=7.3,  # default nudging zone
            shapiro_zone_width=11.5,  # default shapiro zone
            shapiro_tilt=3.5,  # default abrupt transition in the shapiro zone
            feeder_info_file=(
                '/sciclone/schism10/Hgrid_projects/STOFS3D-v7/v20.0/Feeder/'
                'feeder_heads_bases.xy'),
            mandatory_sources_coor=rsf.v19p2_mandatory_sources_coor,
            nwm_cache_folder=Path('/sciclone/schism10/whuang07/schism20/NWM_v2.1/'),
            bc_flags=[
                [5, 5, 4, 4],  # Atlantic Ocean
                [5, 5, 4, 4],  # Gulf of St. Lawrence
                [0, 1, 2, 2],  # St. Lawrence River
            ],
            bc_const=[
                [None, None, None, None],  # Atlantic Ocean
                [None, None, None, None],  # Gulf of St. Lawrence
                [None, None, 10.0, 0.0],  # St. Lawrence River
            ],
            bc_relax=[  # relaxation timescale for each boundary variable
                [None, None, 0.5, 0.5],  # Atlantic Ocean
                [None, None, 0.5, 0.5],  # Gulf of St. Lawrence
                [None, None, 0.01, 1.0],  # St. Lawrence River
            ],
            tvd_regions=[
                'tvd0_1.reg', 'tvd0_2.reg', 'tvd0_3.reg', 'tvd0_4.reg',
                'tvd0_5.reg', 'tvd0_6.reg', 'tvd0_7.reg',
                'upwind_east_Carribbean.rgn', 'upwind_west_Carribbean.rgn'
            ]
        )

    @classmethod
    def v7_hercules_test(cls):
        '''Factory method to create a configuration for STOFS3D-v7 2D setup'''
        return cls(
            ocean_bnd_ids=[0, 1],
            elev2d_uniform_shift=-0.42,  # add a uniform shift to elev2D
            nudging_zone_width=7.3,  # default nudging zone
            shapiro_zone_width=11.5,  # default shapiro zone
            shapiro_tilt=3.5,  # default abrupt transition in the shapiro zone
            relocate_source=True,  # need the feeder info file
            feeder_info_file=(
                '/work/noaa/nosofs/feiye/STOFS-3D-Atl-Example-Setup/DATA/'
                'Feeder_channels/feeder_heads_bases.xy'),
            mandatory_sources_coor=rsf.v19p2_mandatory_sources_coor,
            nwm_cache_folder=None,
            bc_flags=[
                [3, 3, 0, 0],  # Atlantic Ocean
                [3, 3, 0, 0],  # Gulf of St. Lawrence
                [0, 1, 0, 0],  # St. Lawrence River
            ],
            tvd_regions=[
                'tvd0_1.reg', 'tvd0_2.reg', 'tvd0_3.reg', 'tvd0_4.reg',
                'tvd0_5.reg', 'tvd0_6.reg', 'tvd0_7.reg',
                'upwind_east_Carribbean.rgn', 'upwind_west_Carribbean.rgn'
            ]
        )

    @classmethod
    def v8_louisianna(cls):
        '''Factory method to create a configuration for STOFS3D-v8's local test in Louisianna'''
        return cls(
            ocean_bnd_ids=[0],
            nudging_zone_width=0,  # default nudging zone
            shapiro_zone_width=0,  # default shapiro zone
            shapiro_tilt=0,  # default abrupt transition in the shapiro zone
            feeder_info_file=(
                '/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v43s2_RiverMapper/'
                'v44/Feeder/feeder_heads_bases.xy'
            ),
            hgrid_without_feeders=None,
            mandatory_sources_coor=rsf.v45_s2_mandatory_sources_coor,
            relocate_source=False,
            nwm_cache_folder=None,
            bc_flags=[[5, 3, 0, 0]],
            bc_relax=[[None, None, None, None]],
            bc_const=[[None, None, None, None]],
        )

    @classmethod
    def v7p2(cls):
        '''Factory method to create a configuration for STOFS3D-v7.2 3D setup'''
        return cls(
            ocean_bnd_ids=[0, 1],
            elev2d_uniform_shift=-0.42,  # add a uniform shift to elev2D
            nudging_zone_width=7.3,  # default nudging zone
            shapiro_zone_width=11.5,  # default shapiro zone
            shapiro_tilt=3.5,  # default abrupt transition in the shapiro zone
            feeder_info_file=(
                # '/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v23.3/Feeder/'
                '/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v31/Feeder/'
                'feeder_heads_bases.xy'
            ),
            hgrid_without_feeders='/sciclone/schism10/feiye/STOFS3D-v8/R13p_v7/hgrid.gr3',
            relocate_source=True,
            mandatory_sources_coor=rsf.v19p2_for_sms_v27_mandatory_sources_coor,
            nwm_cache_folder=None,
            source_ele_replace_dict={
                53: 3552194,
                203837: 219533,
                253: 205745,
                277: 236142,
            },
            bc_flags=[
                [5, 5, 4, 4],  # Atlantic Ocean
                [5, 5, 4, 4],  # Gulf of St. Lawrence
                [0, 1, 1, 2],  # St. Lawrence River
            ],
            bc_const=[
                [None, None, None, None],  # Atlantic Ocean
                [None, None, None, None],  # Gulf of St. Lawrence
                [None, None, None, 0.0],  # St. Lawrence River
            ],
            bc_relax=[  # relaxation timescale for each boundary variable
                [None, None, 0.5, 0.5],  # Atlantic Ocean
                [None, None, 0.5, 0.5],  # Gulf of St. Lawrence
                [None, None, 0.01, 1.0],  # St. Lawrence River
            ],
            tvd_regions=[
                'tvd0_1.reg', 'tvd0_2.reg', 'tvd0_3.reg', 'tvd0_4.reg',
                'tvd0_5.reg', 'tvd0_6.reg', 'tvd0_7.reg',
                'upwind_east_Caribbean.rgn', 'upwind_west_Caribbean.rgn',
                # 'upwind_Honduras.reg'
            ]
        )

    @classmethod
    def v7p2_subset(cls):
        '''
        Factory method to create a configuration for STOFS3D-v7.2 3D setup
        for a composite mesh: coarse region + fine region (subset of the v7.2 mesh)
        '''
        return cls(
            ocean_bnd_ids=[0, 1],
            elev2d_uniform_shift=-0.42,  # add a uniform shift to elev2D
            nudging_zone_width=7.3,  # default nudging zone
            shapiro_zone_width=11.5,  # default shapiro zone
            shapiro_tilt=3.5,  # default abrupt transition in the shapiro zone
            feeder_info_file=(
                # '/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v23.3/Feeder/'
                '/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v31/Feeder/'
                'feeder_heads_bases.xy'
            ),
            hgrid_without_feeders=None,
            relocate_source=True,
            mandatory_sources_coor=rsf.v19p2_for_sms_v27_mandatory_sources_coor,
            nwm_cache_folder=None,
            bc_flags=[
                [5, 5, 4, 4],  # Atlantic Ocean
                [5, 5, 4, 4],  # Gulf of St. Lawrence
                [0, 1, 2, 2],  # St. Lawrence River
            ],
            bc_const=[
                [None, None, None, None],  # Atlantic Ocean
                [None, None, None, None],  # Gulf of St. Lawrence
                [None, None, 10.0, 0.0],  # St. Lawrence River
            ],
            bc_relax=[  # relaxation timescale for each boundary variable
                [None, None, 0.5, 0.5],  # Atlantic Ocean
                [None, None, 0.5, 0.5],  # Gulf of St. Lawrence
                [None, None, 0.01, 1.0],  # St. Lawrence River
            ],
            tvd_regions=[
                'tvd0_1.reg', 'tvd0_2.reg', 'tvd0_3.reg', 'tvd0_4.reg',
                'tvd0_5.reg', 'tvd0_6.reg', 'tvd0_7.reg',
                'upwind_east_Caribbean.rgn', 'upwind_west_Caribbean.rgn',
                # 'upwind_Honduras.reg'
            ]
        )

    @classmethod
    def v8(cls):
        '''Factory method to create a configuration for STOFS3D-v8 3D setup'''
        return cls(
            ocean_bnd_ids=[0, 1],
            elev2d_uniform_shift=-0.42,  # add a uniform shift to elev2D
            nudging_zone_width=7.3,  # default nudging zone
            shapiro_zone_width=11.5,  # default shapiro zone
            shapiro_tilt=3.5,  # default abrupt transition in the shapiro zone
            feeder_info_file=(
                # '/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v23.3/Feeder/'
                '/sciclone/schism10/Hgrid_projects/STOFS3D-v8/v31/Feeder/'
                'feeder_heads_bases.xy'
            ),
            relocate_source=True,
            mandatory_sources_coor=rsf.v24p4_mandatory_sources_coor,  # v23p3_mandatory_sources_coor,
            nwm_cache_folder=None,
            bc_flags=[
                [5, 5, 4, 4],  # Atlantic Ocean
                [5, 5, 4, 4],  # Gulf of St. Lawrence
                [0, 1, 2, 2],  # St. Lawrence River
            ],
            bc_const=[
                [None, None, None, None],  # Atlantic Ocean
                [None, None, None, None],  # Gulf of St. Lawrence
                [None, None, 10.0, 0.0],  # St. Lawrence River
            ],
            bc_relax=[  # relaxation timescale for each boundary variable
                [None, None, 0.5, 0.5],  # Atlantic Ocean
                [None, None, 0.5, 0.5],  # Gulf of St. Lawrence
                [None, None, 0.01, 1.0],  # St. Lawrence River
            ],
            tvd_regions=[
                'tvd0_1.reg', 'tvd0_2.reg', 'tvd0_3.reg', 'tvd0_4.reg',
                'tvd0_5.reg', 'tvd0_6.reg', 'tvd0_7.reg',
                'upwind_east_Carribbean.rgn', 'upwind_west_Carribbean.rgn'
            ]
        )

    @classmethod
    def v8_Hercules(cls):
        '''Factory method to create a configuration for STOFS3D-v8 3D setup'''
        return cls(
            elev2d_uniform_shift=-0.42,  # add a uniform shift to elev2D
            nudging_zone_width=7.3,  # default nudging zone
            shapiro_zone_width=11.5,  # default shapiro zone
            shapiro_tilt=3.5,  # default abrupt transition in the shapiro zone
            feeder_info_file=(
                '/work/noaa/nosofs/feiye/STOFS-3D-Atl-Example-Setup/DATA/'
                'Feeder_channels/feeder_heads_bases.xy'),
            mandatory_sources_coor=rsf.v19p2_mandatory_sources_coor,
            nwm_cache_folder=None,
            bc_flags=[
                [5, 5, 4, 4],  # Atlantic Ocean
                [5, 5, 4, 4],  # Gulf of St. Lawrence
                [0, 1, 2, 2],  # St. Lawrence River
            ],
            bc_const=[
                [None, None, None, None],  # Atlantic Ocean
                [None, None, None, None],  # Gulf of St. Lawrence
                [None, None, 10.0, 0.0],  # St. Lawrence River
            ],
            bc_relax=[  # relaxation timescale for each boundary variable
                [None, None, 0.5, 0.5],  # Atlantic Ocean
                [None, None, 0.5, 0.5],  # Gulf of St. Lawrence
                [None, None, 0.01, 1.0],  # St. Lawrence River
            ],
            tvd_regions=[
                'tvd0_1.reg', 'tvd0_2.reg', 'tvd0_3.reg', 'tvd0_4.reg',
                'tvd0_5.reg', 'tvd0_6.reg', 'tvd0_7.reg',
                'upwind_east_Carribbean.rgn', 'upwind_west_Carribbean.rgn'
            ]
        )
