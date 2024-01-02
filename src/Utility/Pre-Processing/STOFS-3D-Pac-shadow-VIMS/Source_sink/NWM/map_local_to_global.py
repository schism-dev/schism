import json

import numpy as np

from pylib import read_schism_hgrid, read_schism_bpfile, inside_polygon

if __name__ == '__main__':

    layers = ['conus', 'alaska', 'hawaii']
    regs = ['columbia.reg', 'cowlitz.reg', 'lewis.reg', 'willamette.reg']

    gd = read_schism_hgrid('hgrid.gr3')
    gd.compute_ctr()

    for layer in layers:
        fname_source = f'sources_{layer}.json'
        fname_sink = f'sinks_{layer}.json'
        sources_fid = json.load(open(fname_source))
        sinks_fid = json.load(open(fname_sink))

        sources = {}
        sinks = {}
        
        local_to_global = json.load(open(f'local_to_global_{layer}.json'))

        #repalce keys with global element ids
        keys_source = list(sources_fid.keys())
        for key in keys_source: 
            sources[local_to_global.get(key)] = sources_fid.get(key)

        keys_sink = list(sinks_fid.keys())
        for key in keys_sink: 
            sinks[local_to_global.get(key)] = sinks_fid.get(key)


        #remove elements near willamette/columbia river
        if layer == 'conus':

            elems = []
            for reg in regs:
               bp = read_schism_bpfile(reg, fmt=1)
               eidxs = np.nonzero(inside_polygon(np.c_[gd.xctr, gd.yctr], bp.x, bp.y) == 1)[0] + 1
               elems.extend(eidxs.tolist())

            for elem in elems:
                if sources.get(str(elem)) is not None: 
                    del sources[str(elem)]
                if sinks.get(str(elem)) is not None:
                    del sinks[str(elem)]

        #write sources/sinks to new json
        with open(f'sources_{layer}_global.json', 'w') as f: 
            json.dump(sources, f)

        with open(f'sinks_{layer}_global.json', 'w') as f: 
            json.dump(sinks, f)
