import os
from osgeo import ogr

def multipoly2poly(in_lyr, out_lyr):
    for in_feat in in_lyr:
        geom = in_feat.GetGeometryRef()
        if geom.GetGeometryName() == 'MULTIPOLYGON':
            for geom_part in geom:
                addPolygon(geom_part.ExportToWkb(), out_lyr)
        else:
            addPolygon(geom.ExportToWkb(), out_lyr)

def addPolygon(simplePolygon, out_lyr):
    featureDefn = out_lyr.GetLayerDefn()
    polygon = ogr.CreateGeometryFromWkb(simplePolygon)
    out_feat = ogr.Feature(featureDefn)
    out_feat.SetGeometry(polygon)
    out_lyr.CreateFeature(out_feat)
    print 'Polygon added.'

from osgeo import gdal
gdal.UseExceptions()
driver = ogr.GetDriverByName('ESRI Shapefile')
in_ds = driver.Open('data/multipoly.shp', 0)
in_lyr = in_ds.GetLayer()
outputshp = 'data/poly.shp'
if os.path.exists(outputshp):
    driver.DeleteDataSource(outputshp)
out_ds = driver.CreateDataSource(outputshp)
out_lyr = out_ds.CreateLayer('poly', geom_type=ogr.wkbPolygon)
multipoly2poly(in_lyr, out_lyr)