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

# depth slices
```python
from rasterio.crs import CRS
crs_4326 = CRS.from_string('EPSG:4326')    
mtisacoords = [138.3479499999999973,-23.5494699999999995, 141.5960300000000132,-18.0271800000000013]
xmin, ymin, xmax, ymax = mtisacoords

from pyproj import CRS

import rasterio
from rasterio.transform import Affine
from rasterio.crs import CRS
from rasterio.transform import from_origin

crs_4326 = CRS.from_string('EPSG:4326')    
mtisacoords = [bounds[0],bounds[1],bounds[2],bounds[3]]
xmin, ymin, xmax, ymax = mtisacoords

def depth_slices(mesh, array_name, save_dir):
    for z in mesh.z:
        print(z)
        slices = mesh.slice(normal=[0,0,1],origin=(mesh.center[0],mesh.center[1],z))
        #print (slices)

        try:
            slices02 = slices.cell_data_to_point_data()
            reshaped_arr_3 = slices02[array_name].reshape(mesh.dimensions[1], mesh.dimensions[0])

            reshaped_arr_3 = np.expand_dims(reshaped_arr_3, axis=0)
            da = xr.DataArray(data=reshaped_arr_3, dims=["band","y","x"], coords={"y":mesh.y,"x":mesh.x})
            #print(da)

            #transform = Affine.translation(gridX[0][0]-rRes/2, gridY[0][0]-rRes/2)*Affine.scale(rRes,rRes)
            transform = from_origin(da.x.min().item(),da.y.max().item(), 1000,1000)
            #print(transform)
            rasterCrs = CRS.from_epsg(3577)
            #print(rasterCrs.data)

            new_dataset = rasterio.open(save_dir + 'Depth_Slice_' + str(int(z)) + '.tif', 'w', driver='GTiff',
                height = mesh.dimensions[1], width = mesh.dimensions[0],
                count=1, dtype=str(reshaped_arr_3.dtype),
                crs=rasterCrs,
                nodata = -100,
                transform=transform)

            new_dataset.write(np.flipud(reshaped_arr_3.squeeze()), 1)
            new_dataset.close()

            raster = rioxarray.open_rasterio(save_dir + 'Depth_Slice_' + str(int(z)) + '.tif')
            raster.write_crs
            xdsc = raster.rio.reproject('EPSG:4326')
            xdsc = xdsc.rio.clip_box(
                minx=xmin,
                miny=ymin,
                maxx=xmax,
                maxy=ymax,
            )

            xdsc.rio.to_raster(save_dir + 'MtIsa_Depth_Slice_' + str(int(z)) + '.tif')

        except Exception as sliceE:
            print(z, sliceE)

        #break
        
        
depth_slices(mesh, 'NAC_gzinv3d.den', 'J:/MtIsa/145901_04_0/04_recombined_models/')
```


## Fence diagrams

https://github.com/pyvista/pyvista-support/issues/272

https://github.com/RichardScottOZ/Transform-2021/blob/main/Seismic-2D-Fence.ipynb


## Thresholding
- https://github.com/pyvista/pyvista/discussions/1815
- salient point is that it is Cellwise, NOT Pointwise

## Point Normals
- need polydata