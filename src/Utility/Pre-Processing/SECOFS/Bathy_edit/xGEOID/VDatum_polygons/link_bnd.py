import os
import glob

if __name__ == '__main__':
    path = '/sciclone/schism10/lcui01/schism20/ICOGS/SECOFS/Scripts/Hgrid/vdatum'
    dst = '/sciclone/schism10/hjyoo/task/task6_SECOFS/simulation/Whole_Domain/Geoid/VDatum_polygons'

    statename = ['NJ', 'DE', 'MD', 'VA', 'NC', 'SC', 'GA', 'FL', 'AL']
   
    dirs = glob.glob(f'{path}/*')

    subdirs = []
    for fname in dirs:
        st = fname.split('/')[-1][:2]
        if st in statename:
            src = f'{fname}/*.bnd'
            print(src)
            files = glob.glob(src)
            os.symlink(files[0], f"{dst}/{files[0].split('/')[-1]}")
