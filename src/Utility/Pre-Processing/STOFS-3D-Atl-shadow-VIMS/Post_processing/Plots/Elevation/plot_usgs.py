from schism_py_pre_post.Plot.plot_usgs import plot_usgs


if __name__ == "__main__":
    plot_usgs(
        station_bp_file='/sciclone/schism10/feiye/STOFS3D-v5/Inputs/Stations/USGS_station.bp',
        model_start_day_str='2021-05-01 00:00:00',
        sec_per_time_unit=86400,
        elev_out_files={'RUN02b': '/sciclone/schism10/feiye/STOFS3D-v5/Outputs/O02b_JZ/fort.18',
                        'RUN01b': '/sciclone/schism10/feiye/STOFS3D-v5/Outputs/O01b_JZ/fort.18'},
        plot_start_day_str='2021-05-17 00:00:00',
        plot_end_day_str='2021-05-22 00:00:00',
        output_dir='$cdir/srv/www/htdocs/yinglong/feiye/STOFS3D-v5/RUN02b/USGS/',
    )
