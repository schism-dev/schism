"""
Configuration for STOFS-3D-ATL model,
used by stofs3d_atl_driver.py
"""


from pathlib import Path
from datetime import datetime

import Source_sink.Relocate.relocate_source_feeder as rsf


# ---------------------------------------------------------------------
#                               Classes
# ---------------------------------------------------------------------
class ConfigStofs3dAtlantic():
    '''A class to handle the configuration of STOFS-3D-ATL model,
    i.e., processing the parameters and storing the factory settings.
    '''
    def __init__(
        self,
        startdate=datetime(2017, 12, 1),  # start date of the model
        rnday=60,  # number of days to run the model
        ocean_bnd_ids=[0],  # list of open boundary ids for *D.th.nc
        elev2d_uniform_shift=0.0,  # add a uniform shift to elev2D
        nudging_zone_width=1.5,  # in degrees
        nudging_day=1.0,  # in days
        shapiro_zone_width=2.5,  # in degrees
        shapiro_tilt=2.0,  # more abrupt transition in the shapiro zone
        bc_flags=None,  # a list of lists, each sublist is a set of flags for an open boundary
        bc_const=None,  # a list of lists, each sublist is a set of constants for Eta, Vel, T, S at teach open boundary
        bc_relax=None,  # a list of lists, each sublist is a set of relaxations for Eta, Vel, T, S at each open boundary
        nwm_cache_folder=None,
        relocate_source=True,
        feeder_info_file=None,  # file containing feeder info,
                                # made by make_feeder_channel.py in RiverMapper
        no_feeder_hgrid=None,
        mandatory_sources_coor=None,  # a dictionary of mandatory sources' coordinates
        gr3_values=None,
        tvd_regions=None
    ):

        self.startdate = startdate
        self.rnday = rnday
        self.ocean_bnd_ids = ocean_bnd_ids
        self.elev2d_uniform_shift = elev2d_uniform_shift
        self.nudging_zone_width = nudging_zone_width
        self.nudging_day = nudging_day
        self.shapiro_zone_width = shapiro_zone_width
        self.shapiro_tilt = shapiro_tilt
        self.relocate_source = relocate_source
        self.nwm_cache_folder = nwm_cache_folder
        self.feeder_info_file = feeder_info_file
        self.no_feeder_hgrid = no_feeder_hgrid
        self.mandatory_sources_coor = mandatory_sources_coor

        if bc_flags is None:
            self.bc_flags = [
                [3, 3, 0, 0],  # tides for elev and vel
                [3, 3, 0, 0],  # tides for elev and vel
                [3, 3, 0, 0],  # tides for elev and vel
                [3, 3, 0, 0],  # tides for elev and vel
                [3, 3, 0, 0],  # tides for elev and vel
            ]
        else:
            self.bc_flags = bc_flags

        if bc_const is None:
            self.bc_const = [
                [None, None, None, None],
                [None, None, None, None],
                [None, None, None, None],
                [None, None, None, None],
                [None, None, None, None],
            ]
        else:
            self.bc_const = bc_const
            
        if bc_relax is None:
            self.bc_relax = [
                [None, None, None, None],
                [None, None, None, None],
                [None, None, None, None],
                [None, None, None, None],
                [None, None, None, None],
            ]
        else:
            self.bc_relax = bc_relax

        if gr3_values is None:
            self.gr3_values = {  # uniform gr3 values
                'albedo': 0.1,
                'diffmax': 1.0,
                'diffmin': 1e-6,
                'watertype': 1.0,
                'windrot_geo2proj': 0.0
            }
        else:
            self.gr3_values = gr3_values

        if tvd_regions is None:
            self.tvd_regions = []
        else:
            self.tvd_regions = tvd_regions

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
            feeder_info_file='',
            relocate_source=False,
            nwm_cache_folder=Path('/sciclone/schism10/whuang07/schism20/NWM_v2.1/'),
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
            no_feeder_hgrid='/sciclone/schism10/feiye/STOFS3D-v8/R13p_v7/hgrid.gr3',
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
                'upwind_east_Carribbean.rgn', 'upwind_west_Carribbean.rgn'
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
