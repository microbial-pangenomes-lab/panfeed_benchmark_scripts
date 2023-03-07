import pandas as pd
from collections import defaultdict
from math import ceil
import shutil
import random
import os
import io


def check_bmds(path, bmsets, check, gffpath):
    #checks whether there are more lineages than there are strains in the smallest subset
    print("checking benchmark sets for length")
    lineagedict = defaultdict(list)
    linfile = open(path, "r")
    
    for idx, line in enumerate(linfile):
        strain = line.split("\t")[0]
        lineage = line.split("\t")[1].rstrip("\n")
        print(strain, lineage)
        lineagedict[f"{lineage}"].append(f"{strain}")   
    
    if check == True:
        for ds in bmsets:
            if ds < len(lineagedict):
                bmsets.remove(ds)
    gffs = set()
    
    for file in os.listdir(gffpath):
        if file.endswith(".fasta" or "fasta"):
            gffs.add(file.split(".")[0])
    
    for lin in lineagedict:
        for strain in lineagedict[lin]:
            if strain not in gffs:
                print(f"{strain} is not in gff directory")
    
    print("done with bm set checking")
    return lineagedict, bmsets, idx


def stroi_prep(bmnr, lindict, strainnr, outpath, gffpath):
    print(f"preparing stroi file for dataset {bmnr}")
    os.mkdir(f"{outpath}{bmnr}")
    
    linlist = []
    #calculates the number of how many strains a lineage contributes to the overall number of strains
    for lin in lindict:
        num2pick = len(lindict[lin]) / strainnr * bmnr
        linlist.append(ceil(num2pick))
        
    lins = list(lindict.keys())
    gfflist = []
    strainlist = []
    print(lins)
    print(linlist)
    print(len(lins), len(linlist))
    
    for lin in range(len(linlist)):
        print("lineage", lins[lin])

        #take random sample of strains in a lineage
        if len(lindict[lins[lin]]) < linlist[lin]:
            strains = random.sample(lindict[lins[lin]], len(lindict[lins[lin]]))
        else:
            strains = random.sample(lindict[lins[lin]], linlist[lin])
        
        for strain in strains:
        #pick a random strain from the random lineage sample to add to the strains of interest (is done once for every lineage)
            print("strain", strain)
            gfflist.append(f"{strain}.fasta")
            strainlist.append(f"{strain}")
                
    print("done with all lineages")
    #since we are rounding (ceil()) the number of strains to pick for each lineage, there are a larger number of strains than we actually want
    #we therefore delete the excess numbers by randomly picking strains (that are not strains of interest), regardless of their lineage.
    while len(strainlist) > bmnr:
        popstrain = random.choice(strainlist)
        

        strainlist.remove(popstrain)
        gfflist.remove(popstrain + ".fasta")

    
    for gff in gfflist:
        shutil.copy(f"{gffpath}{gff}", f"{outpath}{bmnr}/")

    stroimem = io.StringIO()
    stroimem.write("ID\tPath\n")
    for gff in gfflist:
        stroimem.write(f"{gff.rstrip('.fasta')}\t{outpath}{bmnr}/{gff}\n")
    
    ucin = open(f"{outpath}{bmnr}/uc_in.txt", "w")
    ucin.write(stroimem.getvalue())
    ucin.close()
    
    return strainlist