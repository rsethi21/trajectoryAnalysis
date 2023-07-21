import concurrent.futures
import json
import os
import re
import sys
import time
from itertools import repeat
import argparse

import numpy as np
import pytraj as pt

## specific to the protein

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="file with trajectory information", required=True)
parser.add_argument("-o", "--outpath", help="output folder path", required=True)
parser.add_argument("-l", "--length", help="length of protein (in number of atoms)", required=True)
parser.add_argument("-n", "--number", help="number of structures", required=True)
parser.add_argument("-w", "--workers", help="workers for multithreading", required=False, default=4)
parser.add_argument("-i", "--ion", help="mask of ions for salt mask", required=False, default=None)
parser.add_argument("-t", "--hydration", help="mask of hydration for water mask", required=False, default=None)


d = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
}


def LoadTraj(caseToProcess):
    ## load

    traj = pt.iterload(caseToProcess)

    ## get structure info

    mask = ":MET@BB"
    fr_i = traj[0, mask]
    l = pt.get_coordinates(fr_i)
    nStruct = np.shape(l)[1]

    return nStruct, traj


def GetTrajData(traj, nStruct, protLen, ionMask, waterMask):

    timeseriesContainer = dict()

    start = ((nStruct-1) * protLen) + 1
    fin = start + protLen - 1
    mask = "@%d-%d" % (start, fin)  ############ find out what the d's mean

    com = pt.center_of_mass(traj, mask=mask)

    ## compute radgyr

    data = pt.radgyr(traj, mask=mask)

    ## compute rmsd to ref first frame

    rmsf = pt.rmsf(traj, mask=mask)

    ## compute number of sodium ions in solvation layers

    if ionMask != None:
        sodiumshell = pt.watershell(traj, solute_mask=mask, solvent_mask=ionMask)
        timeseriesContainer["Salt"] = sodiumshell.values.tolist()

    ## compute number of watermolecules in solvation layers
    
    if waterMask != None:
        watershell = pt.watershell(traj, solute_mask=mask, solvent_mask=waterMask)
        timeseriesContainer["Hydration"] = watershell.values.tolist()

    timeseriesContainer["RgSeries"] = data.tolist()
    timeseriesContainer["RMSF"] = rmsf.tolist()
    timeseriesContainer["COM"] = com.tolist()

    return timeseriesContainer


def doit(nStruct, protLen, inputFile, outpath, ionMask, waterMask, workers):
    
    dataSeries = os.path.join(outpath, "trajMetrics.json")

    print("Generating data from trajs")
    nStructPossible, traj = LoadTraj(caseToProcess)
    nStruct = np.min([nStruct, nStructPossible])
    print(f"Processing {nStruct} out of {nStructPossible}")

    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        indices = [i for i in range(1, nStruct + 1)]
        results = executor.map(GetTrajData, repeat(traj), indices)

    seriesData = [result for result in results]

    with open(dataSeries, "w") as file:
        json.dump(seriesData, file)


if __name__ == "__main__":

    args = parser.parse_args()        

    start = time.perf_counter()
    doit(args.number, args.length, args.file, args.outpath, args.ion, args.hydration, args.workers)
    end = time.perf_counter()

    print(f"Time to run = {end-start}")

