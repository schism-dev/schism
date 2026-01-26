from pyschism.mesh.hgrid import Gr3

if __name__ == '__main__':
    hgrid = Gr3.open('./hgrid.gr3', crs='epsg:4326')
    grd = hgrid.copy()

    gr3_names = ['albedo', 'diffmax', 'diffmin', 'watertype', 'windrot_geo2proj']
    values = [0.1, 1.0, 1e-6, 1.0, 0.0]

    for name, value in zip(gr3_names, values):
       grd.description = name
       grd.nodes.values[:] = value
       grd.write(f'{name}.gr3', overwrite=True)

