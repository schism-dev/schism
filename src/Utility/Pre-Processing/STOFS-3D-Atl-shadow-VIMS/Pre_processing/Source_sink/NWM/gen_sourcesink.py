from datetime import datetime
from time import time
import numpy as np
import pathlib
import logging
import json

from pyschism.forcing.source_sink.nwm import NationalWaterModel, NWMElementPairings
from pyschism.mesh import Hgrid

logging.basicConfig(
    format="[%(asctime)s] %(name)s %(levelname)s: %(message)s",
    force=True,
)
logging.captureWarnings(True)

log_level = logging.DEBUG
logging.getLogger('pyschism').setLevel(log_level)


class SourceSinkIn:
    """ class for *.prop or other similar formats"""
    def __init__(self, filename=None, number_of_groups=2, ele_groups=[]):
        self.n_group = number_of_groups  # 0: source; 1: sink
        if filename is not None:
            """Read a usgs data file and initialize from its content """
            self.source_file = filename
            self.np_group = []
            self.ip_group = []

            with open(self.source_file) as fin:
                for k in range(0, self.n_group):
                    self.np_group.append(int(fin.readline().split()[0]))
                    print("Points in Group " + str(k+1) + ": " + str(self.np_group[k]))
                    self.ip_group.append(np.empty((self.np_group[k]), dtype=int))
                    for i in range(0, self.np_group[k]):
                        self.ip_group[k][i] = int(fin.readline())
                    fin.readline()  # empty line between groups
                    if self.np_group[k] > 0:
                        print("p first: " + str(self.ip_group[k][0]))
                        print("p last: " + str(self.ip_group[k][-1]))
        else:
            self.np_group = [len(x) for x in ele_groups]
            self.ip_group = [np.array(x) for x in ele_groups]

    def writer(self, filename=None):
        if filename is None:
            filename = self.source_file

        with open(filename, 'w') as fout:
            for k in range(0, self.n_group):
                print("Points in Group " + str(k+1) + ": " + str(self.np_group[k]))
                fout.write(f"{self.np_group[k]}\n")
                for i in range(0, self.np_group[k]):
                    fout.write(f"{self.ip_group[k][i]}\n")
                fout.write("\n")  # empty line

def reorder_parings(wdir):
    source_sink_in = SourceSinkIn(f'{wdir}/source_sink.in')

    source_parings = json.load(open(f'{wdir}/sources.json'))
    sinks_pairings = json.load(open(f'{wdir}/sinks.json'))

    assert np.array_equal( source_sink_in.ip_group[0], sorted(np.array(list(source_parings.keys())).astype(int)) ), 'Source parings do not match'
    assert np.array_equal( source_sink_in.ip_group[1], sorted(np.array(list(sinks_pairings.keys())).astype(int)) ), 'Sink parings do not match'

    # reorder source/sink parings
    ordered_source_parings = {k: source_parings[k] for k in source_sink_in.ip_group[0].astype(str)}
    ordered_sinks_pairings = {k: sinks_pairings[k] for k in source_sink_in.ip_group[1].astype(str)}
    
    # write to json
    json.dump(ordered_source_parings, open(f'{wdir}/sources_conus.json', 'w'), indent=4)
    json.dump(ordered_sinks_pairings, open(f'{wdir}/sinks_conus.json', 'w'), indent=4)

    pass
    

def gen_sourcesink(startdate:datetime, rnday:float, cache_folder:pathlib.Path=None):
    hgrid = Hgrid.open("./hgrid.gr3", crs="epsg:4326")

    t0 = time()

    #source/sink json files, if not exist, it will call NWMElementPairings to generate.
    sources_pairings = pathlib.Path.cwd() / 'sources.json'
    sinks_pairings = pathlib.Path.cwd() / 'sinks.json'
    output_directory = pathlib.Path.cwd()

    #input directory which saves nc files
    if cache_folder is not None and cache_folder.exists():  # if cache folder exists, use it
        # cache = pathlib.Path('./NWM_v2.1')
        cache = cache_folder
    else:  # if cache folder does not exist, create it
        cache = pathlib.Path(f'./{startdate.strftime("%Y%m%d")}')
        cache.mkdir(exist_ok=True, parents=True)

    # check if source/sink json file exists
    if all([sources_pairings.is_file(), sinks_pairings.is_file()]) is False:
        pairings = NWMElementPairings(hgrid)
        sources_pairings.parent.mkdir(exist_ok=True, parents=True)
        pairings.save_json(sources=sources_pairings, sinks=sinks_pairings)
    else:
        pairings = NWMElementPairings.load_json(
            hgrid,
            sources_pairings,
            sinks_pairings)

    #check nc files, if not exist will download
    nwm=NationalWaterModel(pairings=pairings, cache=cache)

    nwm.write(output_directory, hgrid, startdate, rnday, overwrite=True)
    print(f'It took {time()-t0} seconds to generate source/sink')

if __name__ == '__main__':
    # gen_sourcesink(datetime(2017, 12, 1), 10)

    wdir = '/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I12x/Source_sink/original_source_sink/'
    reorder_parings(wdir)