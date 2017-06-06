import re

saveLine = False

with open("log.lammps",'r') as infile:  # open  this file and call it infile
	with open("output.txt",'w') as outfile:
		for line in infile:  # for each line in infile
			# Ending condition
			if re.match("^Loop time.*\n",line):
				saveLine = False
			if saveLine:	
				if len(line.split()) == 11: # look for the lines with 11 items
					outfile.write(line.split()[4] + "\n")
			# Starting condition		
			if re.match("^Step Temp Press Pxy v_visc f_vave f_pave v_xyrate PotEng KinEng TotEng\s*\n",line):
				saveLine = True
		
