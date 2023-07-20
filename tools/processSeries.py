import concurrent.futures
import json
import os
import re
import sys
import time
from itertools import repeat

import numpy as np
import pytraj as pt

## specific to the protein

protLen = 202
mask = ":MET@BB"  ## i can't get resid 1 to work correctly, so this is a workaround
nStruct = 100

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

    if "pdb" not in caseToProcess:
        trajName = caseToProcess + ".xtc"
        protName = caseToProcess + ".pdb"
        traj = pt.iterload(trajName, protName)
    else:
        traj = pt.iterload(caseToProcess)

    ## get structure info

    fr_i = traj[0, mask]
    l = pt.get_coordinates(fr_i)
    nStruct = np.shape(l)[1]
    print("Finding %d structures" % nStruct)

    return nStruct, traj


def GetTrajData(traj, nStruct):

    copies = []
    timeseriesCopies = []

    start = ((nStruct-1) * protLen) + 1
    fin = start + protLen - 1
    mask = "@%d-%d" % (start, fin)  ############ find out what the d's mean

    com = pt.center_of_mass(traj, mask=mask)

    ## compute radgyr

    data = pt.radgyr(traj, mask=mask)

    ## compute rmsd to ref first frame

    rmsf = pt.rmsf(traj, mask=mask)

    ## compute number of sodium ions in solvation layers

    sodiumshell = pt.watershell(traj, solute_mask=mask, solvent_mask=":ION")

    ## compute number of watermolecules in solvation layers

    watershell = pt.watershell(traj, solute_mask=mask, solvent_mask=":W")

    timeseriesContainer = dict()
    timeseriesContainer["RgSeries"] = data.tolist()
    timeseriesContainer["RMSF"] = rmsd.tolist()
    timeseriesContainer["Salt"] = sodiumshell.values.tolist()
    timeseriesContainer["Hydration"] = watershell.values.tolist()
    timeseriesContainer["COM"] = com.tolist()

    return timeseriesContainer


def doit(nStruct, mode=None, case=None):
    if "trajs3" in case:
        caseToProcess = os.path.join(case, "system_reduced_all.pdb")
        dataSeries = "./traj3Metrics.json"
    elif "trajs7" in case:
        caseToProcess = os.path.join(case, "system_reduced_all_MDtraj.pdb")
        dataSeries = "./data/traj7Metrics.json"
    else:
        raise RuntimeError("dunno this case")

    ## inputs

    if mode is "generation":
        print("Generating data from trajs")
        nStructPossible, traj = LoadTraj(caseToProcess)
        nStruct = np.min([nStruct, nStructPossible])
        print("Processing %d" % nStruct)

        ## saving series data

        with concurrent.futures.ProcessPoolExecutor(max_workers=15) as executor:
            indices = [i for i in range(1, nStruct + 1)]
            results = executor.map(GetTrajData, repeat(traj), indices)

        seriesData = [result for result in results]

        with open(dataSeries, "w") as file:
            json.dump(seriesData, file)

    else:
        raise RuntimeError("mode not understood")


def helpmsg():
    scriptName = sys.argv[0]
    msg = """
Purpose: 
 
Usage:
"""
    msg += "  %s -validation" % (scriptName)
    msg += """

Notes:

"""
    return msg


## main code running script when run from commandline

mode = "generation"
if __name__ == "__main__":
    msg = helpmsg()
    remap = "none"

    if len(sys.argv) < 2:
        raise RuntimeError(msg)

    # Loops over each argument in the command line
    for i, arg in enumerate(sys.argv):
        # print(arg)
        # calls 'doit' with the next argument following the argument '-validation'
        if arg == "-generation":
            mode = "generation"
        if arg == "-postprocess":
            mode = "postprocess"
        if arg == "-nstruct":
            nStruct = int(sys.argv[i + 1])
        if arg == "-case":
            arg1 = sys.argv[i + 1]
            doit(nStruct, mode=mode, case=arg1)

    raise RuntimeError("Arguments not understood")
