from pylib_essentials.schism_file import schism_grid

def gen_nudge(hgrid:schism_grid, outdir, rlmax = 1.5, rnu_day=0.25):
        """
        schism_grid (from pylib) must be in lon/lat coordinates
        set up nudge zone within rlmax distance from the ocean boundary;
        modify the nudging zone width rlmax.
        rlmax can be a uniform value, e.g., rl_max = 1.5,
        or a 2D array of the same size as the hgrid's number of nodes,
        e.g., rl_max = np.zeros(NP) with further tuning of the nudging zone width.
        """

        rnu_max = 1.0 / rnu_day / 86400.0

        #get nudge zone
        lon = hgrid.x
        lat = hgrid.y
        nudge_coeff = np.zeros(hgrid.x.shape, dtype=float)

        global_idxs = {}

        t0 = time()

        nudge_coeff = np.array([0.0]*NP, dtype=float)

        for i in self.ocean_bnd_ids:
            print(f'boundary {i}')
            bnd_idxs = gdf.iloc[i].indexes

            dis = abs((lon + 1j*lat)[:, None] - (lon[bnd_idxs] + 1j*lat[bnd_idxs])[None, :]).min(axis=1)

            out = (1-dis/rlmax)*rnu_max
            out[out<0] = 0
            out[out>rnu_max] = rnu_max
            fp = out>0
            nudge_coeff[fp] = np.maximum(out[fp], nudge_coeff[fp])

            idxs_nudge=np.zeros(NP, dtype=int)
            idxs=np.where(out > 0)[0]
            idxs_nudge[idxs]=1
  
            #expand nudging marker to neighbor nodes
            i34 = self.hgrid.elements.i34
            fp = i34==3
            idxs=np.where(np.max(out[elnode[fp, 0:3]], axis=1) > 0)[0]
            idxs_nudge[elnode[fp,0:3][idxs,:]]=1
            idxs=np.where(np.max(out[elnode[~fp, :]], axis=1) > 0)[0]
            idxs_nudge[elnode[~fp,:][idxs,:]]=1

            idxs=np.where(idxs_nudge == 1)[0]
            global_idxs[i] = idxs


        #logger.info(f'len of nudge idxs is {len(idxs)}')
        logger.info(f'It took {time() -t0} sencods to calcuate nudge coefficient')

        nudge = [f"rlmax={rlmax}, rnu_day={rnu_day}"]
        nudge.extend("\n")
        nudge.append(f"{NE} {NP}")
        nudge.extend("\n")
        hgrid = self.hgrid.to_dict()
        nodes = hgrid['nodes']
        elements = hgrid['elements']
        for idn, (coords, values) in nodes.items():
            line = [f"{idn}"]
            line.extend([f"{x:<.7e}" for x in coords])
            line.extend([f"{nudge_coeff[int(idn)-1]:<.7e}"])
            line.extend("\n")
            nudge.append(" ".join(line))

        for id, element in elements.items():
            line = [f"{id}"]
            line.append(f"{len(element)}")
            line.extend([f"{e}" for e in element])
            line.extend("\n")
            nudge.append(" ".join(line))

        with open(outdir / 'TEM_nudge.gr3','w+') as fid:
            fid.writelines(nudge)

        shutil.copy2(outdir / 'TEM_nudge.gr3', outdir / 'SAL_nudge.gr3')

        return global_idxs