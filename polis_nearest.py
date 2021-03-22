import sys
from copy import copy
from collections import Counter

import click
import fiona
from shapely import geometry
from rtree import index


def compare_polys(poly_a, poly_b):
    """Compares two polygons via the "polis" distance metric.

    See "A Metric for Polygon Comparison and Building Extraction
    Evaluation" by J. Avbelj, et al.

    Input:
        poly_a: A Shapely polygon.
        poly_b: Another Shapely polygon.

    Returns:
        The "polis" distance between these two polygons.
    """
    bndry_a, bndry_b = poly_a.exterior, poly_b.exterior
    dist = polis(bndry_a.coords, bndry_b)
    dist += polis(bndry_b.coords, bndry_a)
    return dist


def polis(coords, bndry):
    """Computes one side of the "polis" metric.

    Input:
        coords: A Shapley coordinate sequence (presumably the vertices
                of a polygon).
        bndry: A Shapely linestring (presumably the boundary of
        another polygon).

    Returns:
        The "polis" metric for this pair.  You usually compute this in
        both directions to preserve symmetry.
    """
    sum = 0.0
    # Skip the last point (same as first)
    for pt in (geometry.Point(c) for c in coords[:-1]):
        # compare distance between point and polygon
        distoline = bndry.distance(pt)
        # besides point to line, still need to calculate min distance of point to point by xiaoyu
        mindistance = distoline
        for pt2 in (geometry.Point(p2) for p2 in bndry.coords[:-1]):
            dispp = pt.distance(pt2)
            if dispp < mindistance:
                mindistance = dispp
        if mindistance < distoline:
            distoline = mindistance

        sum += distoline
    return sum/float(2*len(coords))


def shp_to_list(shpfile):
    """Dumps all the geometries from the shapefile into a list.

    This makes it quick to build the spatial index!
    """
    with fiona.open(shpfile) as src:
        return [geometry.shape(rec['geometry']) for rec in src]


def path_leaf(path):
    import ntpath
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

# @click.command()
# @click.argument('in_ref', type=click.Path(exists=True, dir_okay=False, resolve_path=True))
# @click.argument('in_cmp', type=click.Path(exists=True, dir_okay=False, resolve_path=True))
# @click.argument('out', type=click.Path(dir_okay=False, writable=True, resolve_path=True))


def score(in_ref, in_cmp, out):
    """Given two polygon vector files, calculate the polis score
    between them. The third argument specifies the output file, which
    contains the same geometries as in_cmp, but with the polis score
    assigned to the geometry.

    We consider the first vector file to be the reference.

    Example (from the root of the git repo):

    $polis data/FullSubset.shp data/user_data.shp out.shp
    """
    # Read in all the geometries in the reference shapefile.
    ref_polys = shp_to_list(in_ref)
    # Build a spatial index of the reference shapes.
    # https://toblerity.org/rtree/class.html#rtree.index.Index
    #it is possible that there is no polygon in reference file, so need add a condition here
    if len(ref_polys) > 0:
        idx = index.Index((i, geom.bounds, None)
                          for i, geom in enumerate(ref_polys))
        # Input data to measure.
        hits = []
        # sum of all the min polis score
        sumscore = 0
        sumverticsdiff = 0
        with fiona.open(in_cmp) as src:
            meta = copy(src.meta)
            meta['schema']['properties'] = {
                'polis': 'float:15.2', 'vertices': 'int', 'vertices_ref': 'int', 'diff_ver': 'int'}
            with fiona.open(out, 'w', **src.meta) as sink:  # create a new file
                for rec in src:
                    cmp_poly = geometry.shape(rec['geometry'])
                    # give an index to compare polygon boundbox
                    # Returns the k-nearest objects to the given coordinates.default is one polygon
                    # return id of reference polygon which nearest to the bounds of compare polygon
                    #nearest(coordinates, num_results=1, objects=False) If two index entries are equidistant, both are returned.
                    ref_pindices = [i for i in idx.nearest(cmp_poly.bounds)]

                    # Limit how many we check, if given an excessive
                    # number of ties (aka if someone put in a huge
                    # bounding box that covered a ton of geometries.)
                    if len(ref_pindices) > 5:
                        ref_pindices = ref_pindices[:5]
                        # filter the polygons through the qquary window( cmp_poly.bounds)
                    scores = [compare_polys(cmp_poly, ref_polys[i])
                                  for i in ref_pindices]

                    polis_score = min(scores)
                        # record the nearest polygon id in reference polygon file
                    minscoreindex = ref_pindices[scores.index(polis_score)]
                    hits.append(minscoreindex)
                        # problems
                    vertices_ref = len(
                            ref_polys[minscoreindex].exterior.coords)
                    vertices = len(cmp_poly.exterior.coords)
                    vertices_dif = vertices-vertices_ref

                    sink.write({'geometry': rec['geometry'],
                            'properties': {'polis': polis_score,
                                           'vertices': vertices,
                                           'vertices_ref': vertices_ref,
                                           'diff_ver': vertices_dif
                                           },
                            })
                    sumscore = sumscore + polis_score
                    sumverticsdiff = sumverticsdiff + vertices_dif
                        # when hits ==0 ???
                # if len(hits) == 0:
                #     averagescore = 999
                #     aveverticdiff = 999
                # else:
                #     averagescore = sumscore/len(hits)
                #     aveverticdiff = sumverticsdiff/len(hits)
            averagescore = sumscore/len(hits)
            aveverticdiff = sumverticsdiff/len(hits)
        # Summarize results.
        # how many polygons in compare polygons find their correspond polygon in reference file,it could be one reference to multiple cmp polygons
        print("Number of matches: {}".format(len(hits)))
        print("Number of misses: {}".format(len(ref_polys) - len(hits)))
        print("Number of misses: {}".format(
            len(shp_to_list(in_cmp)) - len(hits)))
        # if the polygon in reference polygon file is match more than once, it is called Duplicate matches by xiaoyu
        print("Duplicate matches: {}".format(
            sum([1 for i in Counter(hits).values() if i > 1])))
    else:
        averagescore = 999
        aveverticdiff = 999

    # or use dictionay, key is filename, value is score
    filename = in_cmp
    return filename, averagescore, aveverticdiff


@click.command()
@click.argument('in_refdir', type=click.Path(exists=True, dir_okay=True, resolve_path=True))
@click.argument('in_cmpdir', type=click.Path(exists=True, dir_okay=True, resolve_path=True))
def ComputeAverPolis(in_refdir, in_cmpdir):

    import csv
    import glob
    import os
    from os import listdir
    from os.path import isfile, join
    txtfiles = []
    filenames = []
    averagescores = []
    dictscore = {}
    # cmpfiles = [f for f in listdir(in_cmpdir) if isfile(join(in_cmpdir, f))]
    # for  cmpfile in cmpfiles:
    # _252299_470784.poly_shapefile.acm.tol_0.125.shp
    # for cmpfile in glob.glob(join(in_cmpdir,"*poly_shapefile.acm.tol*.shp")):
    # for cmpfile in glob.glob(join(in_cmpdir,"*poly_shapefile.acm.tol*.shp")):
    # _257419_469248.poly_shapefile.acm.tol_0.125.shp
    # "*poly_shapefile.acm.tol*.shp")):
    for cmpfile in glob.glob(join(in_cmpdir, "*poly_shapefile.acm.tol*.shp")):
        # _252299_469504.poly_shapefile.acm.tol_0.125.shp format of reference ENS_252299_468992.geojson
        # ".shp"  #".geojson"
        referfile = "ENS"+cmpfile.split("\\")[-1].split(".")[0]+".geojson"
        in_refpath = os.path.join(in_refdir, referfile)
        in_cmppath = join(in_cmpdir, cmpfile)
        # "test" need to replate with a meaning ful name and path

        out_path = os.path.dirname(os.path.abspath(in_cmppath))
        out_path = out_path.replace(in_cmppath.split(
            "\\")[-2], in_cmppath.split("\\")[-2]+"_PoLis")
        out_file = join(out_path, cmpfile.split("\\")[-1])
        print(out_file)
        try:
            os.stat(out_path)
        except:
            os.mkdir(out_path)

        filename, averagescore, aveverticdiff = score(
            in_refpath, in_cmppath, out_file)
    # filenames.append(filename)
    # averagescores.append(averagescore)
        onlyfile = path_leaf(filename)
        dictscore[onlyfile] = [averagescore, aveverticdiff]

    split_name = "PoLis"
    stats_filepath = join(out_path, "{}.csv".format(split_name))
    # https://stackoverflow.com/questions/3348460/csv-file-written-with-python-has-blank-lines-between-each-row,add newline
    # otherwise there are blank lines between rows
    stats_file = open(stats_filepath, "w", newline='')
    fnames = ["name", "PoLis", "diffvertices"]
    writer = csv.DictWriter(stats_file, fieldnames=fnames)
    writer.writeheader()
    # for filename, averagescore in filenames,averagescores:
    for filename in dictscore.keys():
        writer.writerow({
            "name": filename,
            "PoLis": dictscore[filename][0],
            "diffvertices": dictscore[filename][1]
        })
    stats_file.close()


if __name__ == "__main__":
    # score()
    ComputeAverPolis()
