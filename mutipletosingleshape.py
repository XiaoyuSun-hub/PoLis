import fiona
from shapely.geometry import shape, mapping
import csv
import glob
import os
from os import listdir
from os.path import isfile, join
txtfiles = []
filenames = []
averagescores = []
dictscore = {}
# in_cmpdir=r"C:/code/results/wholestudybrt/4bbagurban/unet16bcemax200best197/test"
# cmpfiles = [f for f in listdir(in_cmpdir) if isfile(join(in_cmpdir, f))]
# in_refdir = r"C:\code\processing\EnschedeImageDataset\urbanarea\urbanareabag\urbanareabag4B\raw\test\gt_polygonizedsingle"
in_refdir = r"C:\code\processing\EnschedeImageDataset\urbanarea\urbanareabrtbag\raw4b\testurban\gt_polygonizedsingle"
for referfile in glob.glob(join(in_refdir, "*.geojson")):
# in_refdir="polis/data/bag"
    # for referfile in glob.glob(join(in_refdir, "*.geojson")):
    # _252299_469504.poly_shapefile.acm.tol_0.125.shp format of reference ENS_252299_468992.geojson
         # ".shp"  #".geojson"
    # referfile = "ENS"+cmpfile.split("\\")[-1].split(".")[0]+".geojson"
    # in_refpath = os.path.join(in_refdir, referfile)
    # out_path = os.path.dirname(os.path.abspath(in_refpath))
    # out_path = out_path.replace(in_refpath.split(
            #  "\\")[-2], in_refpath.split("\\")[-2]+"_single")
    # out_path = "urbanareabag4B"
    out_path = "urbanareabag4B"
    try:
        os.stat(out_path)
    except:
        os.mkdir(out_path)
    out_file = join(out_path, referfile.split("\\")[-1])
    print(out_file)

# open the original MultiPolygon file
# mutifile = "polis/data/bag/"+"ENS_257163_471552.geojson"
# import fiona
# from fiona.crs import from_epsg
# source= fiona.open('shp/second_shp.shp', 'r', encoding = 'utf-8')

# with fiona.open('tool_shp_geojson/geojson_fiona.json','w',  driver ="GeoJSON", schema=source.schema, encoding = 'utf-8', crs=fiona.crs.from_epsg(4326)) as geojson:
#      geojson.write(feat)
    with fiona.open(referfile) as source:
        # create the new file: the driver and crs are the same
        # for the schema the geometry type is "Polygon" instead
        output_schema = dict(source.schema)  # make an independant copy
        output_schema['geometry'] = "Polygon"
        # for multi in source:
        #     print("xx")
        with fiona.open(out_file, 'w',
                    driver="GeoJSON",  # source.driver,
                    crs=source.crs,
                    schema=output_schema) as output:
        # print(len(list(source)))
         try:
            for multi in source:
                if multi['geometry']["type"] == 'MultiPolygon':
                  #    # extract polygons out of multipolygon
                  #         list = []
                  #         for polygon in boundary:
                  #                 list.append(polygon)
                    # extract each Polygon feature
                     for poly in shape(multi['geometry']):

                 # write the Polygon feature
                         output.write({
                        'properties': multi['properties'],
                        'geometry': mapping(poly)
                        })
                else:
                         output.write({
                        'properties': multi['properties'],
                        'geometry': mapping(shape(multi['geometry']))
                        })
         except:
           print("empty file")
           print(referfile)
