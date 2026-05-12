"""
This is for temporary fix only.
Manually relocate the source element in sources.json and
other affected source/sink files, including source_sink.in
and source.nc
"""

import json
from copy import deepcopy
from pylib_experimental.schism_file import source_sink


def manual_relocation_json(source_json=None, relocate_dict=None, outdir=None):
    """
    Manually relocate the source element in sources.json and source_sink.in
    """

    # read sources.json
    with open(source_json, 'r', encoding='utf-8') as f:
        ele2fids_dict = json.load(f)
    
    # check if the keys in relocate_dict are in ele2fids_dict
    for ele in relocate_dict.keys():
        if str(ele) not in ele2fids_dict.keys():
            raise ValueError(f"Element {ele} not found in sources.json")

    # build a new dict to ensure order of the keys
    relocated_ele2fids_dict = {}
    for ele, fids in ele2fids_dict.items():
        if int(ele) in relocate_dict:
            new_ele = relocate_dict[int(ele)]
            relocated_ele2fids_dict[new_ele] = fids
            print(f"Relocating {ele} to {new_ele}")
        else:
            relocated_ele2fids_dict[int(ele)] = fids  # keep the original element
    
    # write the new sources.json
    with open(f'{outdir}/sources.json', 'w', encoding='utf-8') as f:
        json.dump(relocated_ele2fids_dict, f, indent=4)
    

def manual_relocation_source_sink(old_source_sink: source_sink, relocate_dict=None, outdir=None):
    """
    Manually relocate the source element in source_sink.in and write new source/sink
    """

    # check if the keys in relocate_dict are in source_sink.in
    for ele in relocate_dict.keys():
        if int(ele) not in old_source_sink.source_sink_in.ip_group[0]:
            raise ValueError(f"Element {ele} not found in source_sink.in")

    source_eles = old_source_sink.source_sink_in.ip_group[0].copy()
    for i, ele in enumerate(source_eles):
        if int(ele) in relocate_dict:
            new_ele = relocate_dict[int(ele)]
            print(f"Relocating {ele} to {new_ele}")
            source_eles[i] = new_ele

    # build a new source_sink.in to ensure order of the keys
    vsource = deepcopy(old_source_sink.vsource)
    vsource.df.columns = [str(ele) for ele in source_eles]

    relocated_ss = source_sink(
        vsource=vsource,
        vsink=old_source_sink.vsink,
        msource=old_source_sink.msource,
    )
    relocated_ss.writer(output_dir=outdir)


if __name__ == "__main__":
    old_ss = source_sink.from_files('/sciclone/schism10/feiye/STOFS3D-v8/I15c_v7/Source_sink/')
    source_json = '/sciclone/schism10/feiye/STOFS3D-v8/I15c_v7/Source_sink/sources.json'

    # relocation dict made by manual inspection in SMS based on vsource.xyz
    # (diagnostic output for source_sink during STOFS3D setup)
    relocate_dict = {
        53: 3557641,
        205: 125969,
        229: 135896,
        253: 201343,
        277: 236142,
        2: 3727442,
        182: 111909,
        152: 337777,
        137: 315422,
        178: 735784,
    }
    manual_relocation_json(
        source_json=source_json,
        relocate_dict=relocate_dict,
        outdir='/sciclone/schism10/feiye/STOFS3D-v8/I15c_v7/Source_sink/Manual_tweaks_isolated_feeders',
    )
    manual_relocation_source_sink(
        old_source_sink=old_ss,
        relocate_dict=relocate_dict,
        outdir='/sciclone/schism10/feiye/STOFS3D-v8/I15c_v7/Source_sink/Manual_tweaks_isolated_feeders/',
    )
    print('Done!')
