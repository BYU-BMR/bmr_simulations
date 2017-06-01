import sys,re

# Find the file to open
for arg in sys.argv:
	if re.match(".*lammpstrj",arg):
		inFileName = arg

mol_id_dict = {}
frame = 0

# Open the file
with open(inFileName) as inFile:
	for line in inFile:
		if re.match("^\d+\s+\d+\s+\d+\s+[\.\d]+\s+[\.\d]+\s+[\.\d]+\s*\n",line):
			fragments = line.split()
			if fragments[0] == '1':
				frame += 1
			if fragments[0] in mol_id_dict:
				if mol_id_dict[fragments[0]] != fragments[1]:
					print("Frame",frame,"Atom",fragments[0],"changed from",
						mol_id_dict[fragments[0]],"to",fragments[1])
			mol_id_dict[fragments[0]] = fragments[1]

print("Total Frames:",frame)
			

