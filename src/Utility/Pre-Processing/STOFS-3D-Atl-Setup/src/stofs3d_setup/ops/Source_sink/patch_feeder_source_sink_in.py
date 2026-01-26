"""
Patch the source_sink.in to fix a few mislocated sources due to isolated feeder channels.
"""

import os
from pylib_experimental.schism_file import SourceSinkIn


def replace_list_item(lst, old_item, new_item):
    if old_item not in lst:
        raise ValueError(f'Item {old_item} not found in list.')
    for i, item in enumerate(lst):
        if item == old_item:
            lst[i] = new_item


def replace_ele_in_source_sink(ss_dir, source_ele_replace_dict):
    """
    Replace element ID in source_sink.in file.
    """

    source_sink_in = SourceSinkIn.from_file(f'{ss_dir}/source_sink.in')
    for old_ele, new_ele in source_ele_replace_dict.items():
        replace_list_item(source_sink_in.ele_groups[0], old_ele, new_ele)
        pass

    source_sink_in.writer(f'{ss_dir}/source_sink_fixed_isolated_feeders.in')
    os.system(f"mv {ss_dir}/source_sink.in {ss_dir}/source_sink.in.0")
    os.system(f"cp {ss_dir}/source_sink_fixed_isolated_feeders.in {ss_dir}/source_sink.in")
    

if __name__ == '__main__':
    ss_dir = '/sciclone/schism10/feiye/STOFS3D-v7.3/I21g/Source_sink/' 
    source_ele_replace_dict = {
        53: 3552194,
        203837: 219533,
        253: 205745,
        277: 236142,
    }

    replace_ele_in_source_sink(ss_dir, source_ele_replace_dict)
