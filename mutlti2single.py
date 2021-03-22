
import fiona
from shapely.geometry import shape, mapping

# open the original MultiPolygon file
with fiona.open('multipolygons.shp') as source:
    # create the new file: the driver and crs are the same
    # for the schema the geometry type is "Polygon" instead
    output_schema = dict(source.schema)  # make an independant copy
    output_schema['geometry'] = "Polygon"

    with fiona.open('output.shp', 'w', 
                    driver=source.driver,
                    crs=source.crs,
                    schema=output_schema) as output:

        # read the input file
        for multi in source:

           # extract each Polygon feature
           for poly in shape(multi['geometry']):

              # write the Polygon feature
              output.write({
                  'properties': multi['properties'],
                  'geometry': mapping(poly)
              })