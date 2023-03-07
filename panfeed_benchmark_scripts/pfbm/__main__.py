from collections import defaultdict
import shutil
import os
from io import StringIO
import argparse

from .bm_mods import check_bmds, stroi_prep, presab_prep
import subprocess



def get_options():
    description = "benchmarking script for panfeed"
    
    parser = argparse.ArgumentParser(description=description)
    
    parser.add_argument("-lin", "--lineages",
                        help="TSV file containing strains and their lineage (no header)")
    
    parser.add_argument("-rep", "--repetitions", type = int,
                        default=1,
                        help="repetitions of benchmark set")
    
    parser.add_argument("-bmset", "--benchmark_sets", nargs="+", type = int,
                        default=[50,100,300,450,912],
                        help="benchmark-set sizes to use")
    
    parser.add_argument("--check_set", 
                        action="store_true",
                        default=False)
    
    parser.add_argument("-outpath", "--output_path",
                        default="./",
                        help="path for output")
    
    parser.add_argument("-presab", "--presence_absence",
                        help="gene_presence_absence.csv from panaroo")
    
    parser.add_argument("-gff", "--gff_directory",
                        help="path with GFFs")
    
    return parser.parse_args()
    
def main():
    
    args = get_options()
    reps = args.repetitions
    
    for rep in range(reps):
        
        bmsets = args.benchmark_sets
        linpath = args.lineages
        check = args.check_set
        outp = args.output_path
        prespath = args.presence_absence
        gffpath = args.gff_directory
        
        lindict, bmsets, strainnr = check_bmds(linpath, bmsets, check, gffpath)
        
        procdict = {}

        #iterate through subsets and start process for each
        for ds in bmsets:
            strainlist = stroi_prep(ds, lindict, strainnr, outp, gffpath)
            presab_prep(ds, strainlist, outp, prespath)

            procdict[f"p{ds}"] = subprocess.Popen(f"/usr/bin/time -v mprof run --include-children -o {outp}{ds}mem{rep}.txt python /PATH/TO/panfeed-runner.py -g {outp}{ds} --cores 4 -p {outp}{ds}/gene_presence_absence.csv -o {outp}{ds}/panfeed_out --targets {outp}{ds}/stroi.txt", 
                                                  shell = True, 
                                                  stdout = subprocess.PIPE,
                                                  stderr = subprocess.STDOUT)

        for ds in bmsets:
            dstring = f"{outp}{ds}time{rep}.txt"
            open(dstring, "w").close()
            outout = open(dstring, "w")
            outbuffer = StringIO()
            #wait for process to finish before starting the next one
            outbuffer.write(procdict[f"p{ds}"].communicate(timeout=None)[0].decode("utf-8"))
            outout.write(outbuffer.getvalue())
            outout.flush()
            outout.close()
            shutil.rmtree(f"{outp}{ds}/")

        tamdict = defaultdict(list)
        #iterate through subsets and gather memory and time consumption data for each run
        for ds in bmsets:

            with open(f"{outp}{ds}time{rep}.txt", "r") as timefile:
                for line in timefile:
                    if "CPU" in line:
                        tamdict[f"{ds}"].append(line.split(":")[1].rstrip("\n"))
                    
                    elif "wall clock" in line:
                        tamdict[f"{ds}"].append(line.split("):")[1].rstrip("\n"))

            with open(f"{outp}{ds}mem{rep}.txt", "r") as memfile:
                next(memfile)
                maxmem = 0
                for line in memfile:
                    splitline = line.split(" ")
                    if maxmem == 0:
                        maxmem = int(float(splitline[1]))
                    elif maxmem < int(float(splitline[1])):
                        maxmem = int(float(splitline[1]))

            tamdict[f"{ds}"].append(maxmem)

        #for each subset durin a repetition, append data to final file
        outtsv = open(f"{outp}benchfile.tsv", "a")
        
        if os.stat(f"{outp}benchfile.tsv").st_size == 0:
            outtsv.write("dataset\tcpu_usage\ttime_spent\tmemory_footprint\n")

        appendbuffer = StringIO()
        for ds in bmsets:
            appendbuffer.write(f"{ds}\t{tamdict[str(ds)][0]}\t{tamdict[str(ds)][1]}\t{tamdict[str(ds)][2]}\n")

            outtsv.write(appendbuffer.getvalue())
            outtsv.flush()
            
            appendbuffer = StringIO()
        outtsv.close()
    
if __name__ == "__main__":
    main()
    
    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        