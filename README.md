# pyvista-notes
some technical notes on some pyvista operations

Grid thesholding - cell operation
See https://github.com/pyvista/pyvista/discussions/1815

Note may want to increase Uniform Grid dimensions by one to put data on the cells

# xarray
```python
data = result['density'].reshape(result.dimensions[2],result.dimensions[1],result.dimensions[0])
XC = [result.bounds[0]+200*i for i in range(result.dimensions[0])]
YC = [result.bounds[2]+200*i for i in range(result.dimensions[1])]
ZC = [result.bounds[4]+200*i for i in range(result.dimensions[2])]
da = xr.DataArray(data=data, dims=["z", "y","x"], coords={"z": ZC, "y": YC, "x": XC})
```

https://github.com/pyvista/pyvista/discussions/2375
see the extremely useful https://github.com/pyvista/pyvista-xarray

# shapely polygons
```pythong
points_3d = np.pad(points_2d, [(0, 0), (0, 1)])  # shape (N, 3)
face = [N + 1] + list(range(N)) + [0]  # cell connectivity for a single cell
polygon = pv.PolyData(points_3d, faces=face)
# extrude along z and plot
body = polygon.extrude((0, 0, 0.5), capping=True)
body.plot(color='white', specular=1)
# extrude along z and plot
```
https://github.com/pyvista/pyvista/discussions/2325

# terrain following mesh
```python
def read_raster_depth(filename):
    """
    Helpful: http://xarray.pydata.org/en/stable/auto_gallery/plot_rasterio.html
    """
    # Read in the data
    data = xr.open_rasterio(filename)
    values = np.asarray(data)
    nans = values == data.nodatavals
    if np.any(nans):
        values = np.ma.masked_where(nans, values)
        values[nans] = np.nan
    # Make a mesh
    xx, yy = np.meshgrid(data['x'], data['y'])
    zz = values.reshape(xx.shape) # will make z-comp the values in the file
    # zz = np.zeros_like(xx) # or this will make it flat
    mesh = pv.StructuredGrid(xx, yy, zz)
    mesh['data'] = values.ravel(order='F')
    #mesh['data'] = mesh['data'] * -1
    return mesh


def extend_structured_ex(topgrid, bottomgrid, exaggerate):
    top = topgrid.points.copy()
    bottom = topgrid.points.copy()
    bottom[:,-1] = bottomgrid.points.copy()[:,-1]*exaggerate

    vol = pv.StructuredGrid()
    vol.points = np.vstack((top, bottom))
    vol.dimensions = [*topgrid.dimensions[0:2], 2]
    vol['base'] = vol.z.ravel(order='F')
    vol.hide_points(np.isnan(vol["base"]))
    return vol

volmes = extend_structured_ex(mesozoic, paleozoic, 50)        
```
