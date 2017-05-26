


import os, re
import sys

#print('Number of arguments:', len(sys.argv), 'arguments.')
#print('Argument List:', str(sys.argv))

# class FrameExtractor():

# 	def()

def extractLastFrame(filename):
	numFrames = 0
	prog = re.compile("^ITEM: TIMESTEP.*")
	with open(filename) as f:
		for line in f:
			if prog.match(line):
				numFrames += 1
	print("numFrames:",numFrames)

	inFile = open(filename)
	outFile = open('last_frame.lammpstrj','w+')
	curFrame = 0
	targetFrame = numFrames
	writeLine = False
	for line in inFile:
		if prog.match(line):
			curFrame += 1
		#if curFrame%40 == 0:
		if curFrame == targetFrame:
			writeLine = True
		else:
			writeLine = False
		if writeLine:
			outFile.write(line)

	outFile.close()
	inFile.close()



# Find all the lammpstrj files in the directory
prog = re.compile(".*lammpstrj")
path = str(sys.argv[1])
files = os.listdir(path)
files = list(filter(prog.match,files))

if len(files) < 1:
	print("There are no .lammpstrj files in this directory.  Cannot do anything.")
elif len(files) > 1:
	print("Multiple lammpstrj files found in this directory. Proceeding with", files[0])
else:
	print("Reading file:", files[0])

filename = files[0]
filepath = path + '/' + filename

extractLastFrame(filepath)


