import sys
from copy import copy
from collections import Counter

import click
import fiona
from shapely import geometry
from rtree import index
import numpy as np


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


def calc(totalref, totalpre, tp_ref, tp_pre):
    totalpre_len = len(totalref)
    totalref_len = len(totalpre)
    precisionsum = 0
    recallsum = 0
    if totalpre_len == 0 or totalref_len == 0:
        return 0, 0
    # calculate the precision recall F1 for different threshold(which is iou)
    for key in tp_ref.keys():
        false_positive = totalpre_len-len(tp_pre[key])
        true_positive = len(tp_ref[key])
        if totalpre_len > 0 and totalref_len > 0:
            precisionsum += true_positive / totalpre_len
            recallsum += true_positive/totalref_len
        else:
            print("true postive:" + str(true_positive))
            print("false_positive:" + str(false_positive))
            print("false_negative:" + str(totalref_len-len(tp_ref[key])))
    # divided by 10 because from 0.5 to 0.95 there are 10 threshold
    AP = precisionsum/10
    AR = recallsum/10
    # if AP+AR > 0:
    #     F1 = 2*AP*AR/(AP+AR)
    # else:
    #     F1 = 0
    return AP, AR


def score(in_ref, in_cmp):
    """Given two polygon vector files, calculate the coco metrics
    between them.We consider the first vector file to be the reference.
    """
    # Threshold = 0.5
    idx = index.Index()
    tp_pre = {}
    tp_ref = {}
    tp_pre_small = {}
    tp_ref_small = {}
    tp_pre_medium = {}
    tp_ref_medium = {}
    tp_pre_large = {}
    tp_ref_large = {}
    pre_small = []
    ref_small = []
    pre_medium = []
    ref_medium = []
    pre_large = []
    ref_large = []
    tp_predicted = []
    tp_reference = []
    precision = {}
    recall = {}
    # F1all = {}
    # https://github.com/cocodataset/cocoapi/blob/master/PythonAPI/pycocotools/cocoeval.py
    # https://cocodataset.org/#detection-eval in coco,Area is measured as the number of pixels in the segmentation mask
    # in our data, the polygon is in RD_NEw coor sys, the unit is meter, which is extracted from the segmentation result of
    # VHR image with the resolution 0.25 m. so we need to convert calculated the area based on pixels to meters in vector data, area for 32*32pixel =32*32*0.25*0.25
    areaperpixel = 0.25**2
    smallarea = 32**2*areaperpixel
    mediumarea = 96**2*areaperpixel
    # iouthrs = np.arange(start=0.5, stop=0.95, step=0.05)
    iouthrs = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
    # precision = np.zero(len(iouthr))
    # recall = np.zero(len(iouthr))
    for threshold in iouthrs:
        tp_ref[str(threshold)] = []
        tp_pre[str(threshold)] = []
        tp_ref_small[str(threshold)] = []
        tp_pre_small[str(threshold)] = []
        tp_ref_medium[str(threshold)] = []
        tp_pre_medium[str(threshold)] = []
        tp_ref_large[str(threshold)] = []
        tp_pre_large[str(threshold)] = []
        precision[str(threshold)] = 0
        recall[str(threshold)] = 0
        # F1all[str(threshold)] = 0

    with fiona.open(in_ref) as refsrc:
        with fiona.open(in_cmp) as cmpsrc:
            for fid, feature in refsrc.items():
                poly = geometry.shape(feature['geometry'])
                idx.insert(fid, poly.bounds)
                # store the polygons in reference file with different size

                if poly.area < smallarea:
                    ref_small = ref_small+[fid]
                if poly.area >= smallarea and poly.area <= mediumarea:
                    ref_medium = ref_medium+[fid]
                if poly.area > mediumarea:
                    ref_large = ref_large+[fid]

            for fid_cmp, rec in cmpsrc.items():
                cmp_poly = geometry.shape(rec['geometry'])
                # store the polygons in predict file with different size
                if cmp_poly.area < smallarea:
                    pre_small = pre_small+[fid_cmp]
                if cmp_poly.area >= smallarea and cmp_poly.area <= mediumarea:
                    pre_medium = pre_medium+[fid_cmp]
                if cmp_poly.area > 96**2*areaperpixel:
                    pre_large = pre_large+[fid_cmp]

                # find the reference polygons intersect with predict polygons by their boundbox
                fids = [int(i) for i in idx.intersection(
                        cmp_poly.bounds)]
                # access the features that those fids reference
                if len(fids) > 0:
                    for fid in fids:
                        feature_ref = refsrc[fid]
                        ref_poly = geometry.shape(feature_ref['geometry'])
                        # check the geometries intersect, not just their bboxs
                        ref_polys = []
                        intersect_fids = []
                        if cmp_poly.intersects(ref_poly):
                            # calculate intersection
                            iou = cmp_poly.intersection(
                                ref_poly).area / cmp_poly.union(ref_poly).area

                            for threshold in iouthrs:
                                if iou >= threshold:
                                    tp_ref[str(threshold)] = tp_ref[str(
                                        threshold)] + [fid]
                                    tp_pre[str(threshold)] = tp_pre[str(
                                        threshold)]+[fid_cmp]
                                    # https://github.com/cocodataset/cocoapi/issues/36,both gt and dt should fufill the area restriction
                                    if cmp_poly.area < smallarea and ref_poly.area <smallarea:
                                        tp_ref_small[str(threshold)] = tp_ref_small[str(
                                            threshold)]+[fid]
                                        tp_pre_small[str(threshold)] = tp_pre_small[str(
                                            threshold)]+[fid]
                                    if cmp_poly.area >= smallarea and cmp_poly.area <= mediumarea and ref_poly.area >= smallarea and ref_poly.area <=mediumarea:
                                        tp_ref_medium[str(threshold)] = tp_ref_medium[str(
                                            threshold)]+[fid]
                                        tp_pre_medium[str(threshold)] = tp_pre_medium[str(
                                            threshold)]+[fid]
                                    if cmp_poly.area > 96 and ref_poly.area > 96:
                                        tp_ref_large[str(threshold)] = tp_ref_large[str(
                                            threshold)]+[fid]
                                        tp_pre_large[str(threshold)] = tp_pre_large[str(
                                            threshold)]+[fid]

            totalpre = len(shp_to_list(in_cmp))
            totalref = len(shp_to_list(in_ref))
            precisionsum = 0
            recallsum = 0
            # calculate the precision recall F1 for different threshold(which is iou)
            for key in tp_ref.keys():
                false_positive = totalpre-len(tp_pre[key])
                true_positive = len(tp_ref[key])
                if totalpre > 0 and totalref > 0:
                    precision[key] = true_positive / totalpre
                    recall[key] = true_positive/totalref
                    precisionsum += precision[key]
                    recallsum += recall[key]
                    # if precision[key] + recall[key] > 0:
                    #     F1all[key] = 2*precision[key] * recall[key] / \
                    #         (precision[key] + recall[key])
                    # else:
                    #     F1all[key] = 0
                else:
                    print("true postive:" + str(true_positive))
                    print("false_positive:" + str(false_positive))
                    print("false_negative:" + str(totalref-len(tp_ref[key])))
                    print(in_cmp)

            AP = precisionsum/len(iouthrs)
            AR = recallsum/len(iouthrs)
            if AP+AR > 0:
                F1 = 2*AP*AR/(AP+AR)
            else:
                F1 = 0
            AP_S, AR_S = calc(ref_small, pre_small, tp_ref_small, tp_pre_small)
            AP_M, AR_M = calc(ref_medium, pre_medium,
                              tp_ref_medium, tp_pre_medium)
            AP_L, AR_L = calc(ref_large, pre_large, tp_ref_large, tp_pre_large)
    filename = in_cmp
    return filename, AP, AR, F1, precision, recall, AP_S, AR_S,  AP_M, AR_M,  AP_L, AR_L


@ click.command()
@ click.argument('in_refdir', type=click.Path(exists=True, dir_okay=True, resolve_path=True))
@ click.argument('in_cmpdir', type=click.Path(exists=True, dir_okay=True, resolve_path=True))
def ComputeMetrics(in_refdir, in_cmpdir):

    import csv
    import glob
    import os
    from os import listdir
    from os.path import isfile, join
    txtfiles = []
    filenames = []
    averagescores = []
    dictscore = {}
    # _257675_469248.poly_shapefile.acm.tol_1.shp
    # for cmpfile in glob.glob(join(in_cmpdir, "_257675_469248.poly_shapefile.acm.tol_1.shp")):
    for cmpfile in glob.glob(join(in_cmpdir, "*poly_shapefile.acm.tol_1.shp")):
        referfile = "ENS"+cmpfile.split("\\")[-1].split(".")[0]+".geojson"
        in_refpath = os.path.join(in_refdir, referfile)
        in_cmppath = join(in_cmpdir, cmpfile)

        filename, AP, AR, F1, precision, recall, AP_S, AR_S,  AP_M, AR_M,  AP_L, AR_L = score(
            in_refpath, in_cmppath)
        onlyfile = path_leaf(filename)
        dictscore[onlyfile] = [AP, AR, F1]+[vals for k,
                                            vals in precision.items()]+[vals2 for k2, vals2 in recall.items()]+[AP_S, AR_S,  AP_M, AR_M,  AP_L, AR_L]

    split_name = "cocometric"
    stats_filepath = join(in_cmpdir, "{}.csv".format(split_name))
    print(stats_filepath)
    # https://stackoverflow.com/questions/3348460/csv-file-written-with-python-has-blank-lines-between-each-row,add newline
    # otherwise there are blank lines between rows
    stats_file = open(stats_filepath, "w", newline='')
    fnames = ["name", "AP", "AR", "F1", "P50", "P55", "P60", "P65", "P70", "P75", "P80", "P85",
              "P90", "P95", "R50", "R55", "R60", "R65", "R70", "R75", "R80", "R85", "R90", "R95", "APS", "ARS",
              "APM", "ARM", "APL", "ARL"]
    writer = csv.DictWriter(stats_file, fieldnames=fnames)
    writer.writeheader()
    # for filename, averagescore in filenames,averagescores:
    Av_AP = 0
    Av_AR = 0
    Av_F1 = 0
    Av_P50 = 0
    Av_P55 = 0
    Av_P60 = 0
    Av_P65 = 0
    Av_P70 = 0
    Av_P75 = 0
    Av_P80 = 0
    Av_P85 = 0
    Av_P90 = 0
    Av_P95 = 0
    Av_P65 = 0
    Av_R50 = 0
    Av_R55 = 0
    Av_R60 = 0
    Av_R65 = 0
    Av_R70 = 0
    Av_R75 = 0
    Av_R80 = 0
    Av_R85 = 0
    Av_R90 = 0
    Av_R95 = 0
    Av_APS = 0
    Av_ARS = 0
    Av_APM = 0
    Av_ARM = 0
    Av_APL = 0
    Av_ARL = 0

    for filename in dictscore.keys():
        print(filename)
        writer.writerow({
            "name": filename,
            "AP": dictscore[filename][0],
            "AR": dictscore[filename][1],
            "F1": dictscore[filename][2],
            "P50": dictscore[filename][3],
            "P55": dictscore[filename][4],
            "P60": dictscore[filename][5],
            "P65": dictscore[filename][6],
            "P70": dictscore[filename][7],
            "P75": dictscore[filename][8],
            "P80": dictscore[filename][9],
            "P85": dictscore[filename][10],
            "P90": dictscore[filename][11],
            "P95": dictscore[filename][12],
            "R50": dictscore[filename][13],
            "R55": dictscore[filename][14],
            "R60": dictscore[filename][15],
            "R65": dictscore[filename][16],
            "R70": dictscore[filename][17],
            "R75": dictscore[filename][18],
            "R80": dictscore[filename][19],
            "R85": dictscore[filename][20],
            "R90": dictscore[filename][21],
            "R95": dictscore[filename][22],
            "APS": dictscore[filename][23],
            "ARS": dictscore[filename][24],
            "APM": dictscore[filename][25],
            "ARM": dictscore[filename][26],
            "APL": dictscore[filename][27],
            "ARL": dictscore[filename][28],
        })
        Av_AP += dictscore[filename][0]
        Av_AR += dictscore[filename][1]
        Av_F1 += dictscore[filename][2]
        Av_P50 += dictscore[filename][3]
        Av_P55 += dictscore[filename][4]
        Av_P60 += dictscore[filename][5]
        Av_P65 += dictscore[filename][6]
        Av_P70 += dictscore[filename][7]
        Av_P75 += dictscore[filename][8]
        Av_P80 += dictscore[filename][9]
        Av_P85 += dictscore[filename][10]
        Av_P90 += dictscore[filename][11]
        Av_P95 += dictscore[filename][12]
        Av_R50 += dictscore[filename][13]
        Av_R55 += dictscore[filename][14]
        Av_R60 += dictscore[filename][15]
        Av_R65 += dictscore[filename][16]
        Av_R70 += dictscore[filename][17]
        Av_R75 += dictscore[filename][18]
        Av_R80 += dictscore[filename][19]
        Av_R85 += dictscore[filename][20]
        Av_R90 += dictscore[filename][21]
        Av_R95 += dictscore[filename][22]
        Av_APS += dictscore[filename][23]
        Av_ARS += dictscore[filename][24]
        Av_APM += dictscore[filename][25]
        Av_ARM += dictscore[filename][26]
        Av_APL += dictscore[filename][27]
        Av_ARL += dictscore[filename][28]

    num_files = len(dictscore.keys())
    writer.writerow({
        "name": "average",
        "AP": Av_AP/num_files,
        "AR": Av_AR/num_files,
        "F1": Av_F1/num_files,
        "P50": Av_P50/num_files,
        "P55": Av_P55/num_files,
        "P60": Av_P60/num_files,
        "P65": Av_P65/num_files,
        "P70": Av_P70/num_files,
        "P75": Av_P75/num_files,
        "P80": Av_P80/num_files,
        "P85": Av_P85/num_files,
        "P90": Av_P90/num_files,
        "P95": Av_P95/num_files,
        "R50": Av_R50/num_files,
        "R55": Av_R55/num_files,
        "R60": Av_R60/num_files,
        "R65": Av_R65/num_files,
        "R70": Av_R70/num_files,
        "R75": Av_R75/num_files,
        "R80": Av_R80/num_files,
        "R85": Av_R85/num_files,
        "R90": Av_R90/num_files,
        "R95": Av_R95/num_files,
        "APS": Av_APS/num_files,
        "ARS": Av_ARS/num_files,
        "APM": Av_APM/num_files,
        "ARM": Av_ARM/num_files,
        "APL": Av_APL/num_files,
        "ARL": Av_ARL/num_files,
    })

    stats_file.close()


if __name__ == "__main__":
    # score()
    ComputeMetrics()
