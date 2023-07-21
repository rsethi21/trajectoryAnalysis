import concurrent.futures
import json
import os
import re
import sys
import time
from itertools import repeat
import argparse

import pickle
import numpy as np
import pandas as pd
import pytraj as pt

## specific to the protein
parser = argparse.ArgumentParser()
parser.add_argument("-l", "--length", help="protein length in number of atoms", required=True)
parser.add_argument("-n", "--number", help="number of protein structures desired (must be less than or equal to max)", required=True)
parser.add_argument("-f", "--file", help="file to trajectory information", required=True)
parser.add_argument("-d", "--distance", help="distance to count as contact", required=False, default=3.4)
parser.add_argument("-o", "--outpath", help="path to folder for output data", required=True)

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


def GetTrajData(traj, nStruct, protLen, path, distance):
    
    start = 1
    fin = protLen*nStruct
    mask = f"@{start}-{fin}&@BB"  ############ find out what the d's mean
    mask2 = f"@{start}-{fin}&@SC1"
    mask3 = f"@{start}-{fin}&@SC2"
    mask4 = f"@{start}-{fin}&@SC3"

    backbone = pt.native_contacts(traj, mask=mask, distance=distance, options="byresidue map mapout matrix series")
    sideChain1 = pt.native_contacts(traj, mask=mask2, distance=distance, options="byresidue map mapout matrix series")
    sideChain2 = pt.native_contacts(traj, mask=mask3, distance=distance, options="byresidue map mapout matrix series")
    sideChain3 = pt.native_contacts(traj, mask=mask4, distance=distance, options="byresidue map mapout matrix series")

    pt.io.to_pickle(backbone, os.path.join(path, f"{nStruct}_backbone.pkl"))
    pt.io.to_pickle(sideChain1, os.path.join(path, f"{nStruct}_sidechain1.pkl"))
    pt.io.to_pickle(sideChain2, os.path.join(path, f"{nStruct}_sidechain2.pkl"))
    pt.io.to_pickle(sideChain3, os.path.join(path, f"{nStruct}_sidechain3.pkl"))

def doit(nStruct, protLen, case, outpath, distance):
    
    print("Generating data from trajs")
    nStructPossible, traj = LoadTraj(caseToProcess)
    nStruct = np.min([nStruct, nStructPossible])
    print(f"Processing {nStruct} out of {nStructPossible}")

    GetTrajData(traj, nStruct, protLen, outpath, distance)

if __name__ == "__main__":
    
    args = parser.parse_args()

    protLen = args.length
    nStruct = args.number

    start = time.perf_counter()
    doit(nStruct, protLen, args.file, args.outpath, args.distance)
    end = time.perf_counter()
    
    print(f"Time to run = {end-start}")

