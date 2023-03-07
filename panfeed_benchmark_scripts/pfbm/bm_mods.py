import pandas as pd
from collections import defaultdict
from math import ceil
import shutil
import random
import os
import io


def check_bmds(path, bmsets, check, gffpath):
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
        if file.endswith(".gff" or "gff"):
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
    stroi = []
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
            gfflist.append(f"{strain}.gff")
            strainlist.append(f"{strain}")
                
        stroi.append(f"{random.choice(strains)}")
        print("done with picking for one lineage")
    print("done with all lineages")
    #since we are rounding (ceil()) the number of strains to pick for each lineage, there are a larger number of strains than we actually want
    #we therefore delete the excess numbers by randomly picking strains (that are not strains of interest), regardless of their lineage.
    while len(strainlist) > bmnr:
        popstrain = random.choice(strainlist)
        
        if popstrain not in stroi:
            strainlist.remove(popstrain)
            gfflist.remove(popstrain + ".gff")
        else:
            continue
    
    for gff in gfflist:
        shutil.copy(f"{gffpath}{gff}", f"{outpath}{bmnr}/")

    stroimem = io.StringIO()
    for strain in stroi:
        stroimem.write(f"{strain}\n")
    
    stroi_out = open(f"{outpath}{bmnr}/stroi.txt", "w")
    stroi_out.write(stroimem.getvalue())
    stroi_out.close()
    
    return strainlist


def presab_prep(bmnr, strainlist, outp, prespath):
    
    presab = pd.read_csv(f"{prespath}", sep=",", index_col=0, header=0, low_memory=False)
    cols = presab.columns
    
    for col in range(2,len(cols)):
        if cols[col] not in strainlist:
            presab.pop(cols[col])
            
    presab.to_csv(f"{outp}{bmnr}/gene_presence_absence.csv")



























