import os
from io import StringIO as sio
import random
import subprocess


fastpath = "/PATH/TO/FASTA/DIRECTORY"

datasets = [250,450,650,912]

reps = 3

strainlist = os.listdir(fastpath)

buffer = sio()

bigbuffer = sio()

for rep in range(reps):
	for ds in datasets:
		strains = random.sample(strainlist, ds)
		with open(f"input{ds}_{rep}.txt", "w") as inp:
			for strain in strains:
				strstrain = strain.split(".")[0]
				buffer.write(f"{strstrain}\t{fastpath}{strain}\n")
			inp.write(buffer.getvalue())
			inp.flush()
		buffer = sio()
		with open("useless.txt", "w") as useless:
			process = subprocess.Popen(f"""
			 							conda run -n fsm-lite /usr/bin/time -o log{ds}_{rep}.txt -v fsm-lite -l input{ds}_{rep}.txt -t tmp
			 							""",
			 							shell = True,
			 							stderr = subprocess.STDOUT,
			 							stdout = subprocess.PIPE)

			useless.write(process.communicate(timeout=None)[0].decode("utf-8"))
		
		os.remove("useless.txt")

