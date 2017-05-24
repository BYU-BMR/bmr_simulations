import os
import errno
import re
from shutil import copyfile

baseName = "shear_dependent"
sizesToRun = [0.2,0.25,0.3,0.35,0.4]
filesToCopy = ['test.in','test.py','submit.sh']

def create_directory(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def editInputFile(path,scale):
	f = open(path, 'r')    # pass an appropriate path of the required file
	lines = f.readlines()
	for i,line in enumerate(lines):
		if re.match(".*dump.*dump_id.*all.*lammpstrj.*",line):
			lines[i] = "dump               dump_id all custom 10000 coating_" + str(scale).replace('.','_') + ".lammpstrj id type xs ys zs\n"
	f.close()   # close the file and reopen in write mode to enable writing to file; you can also open in append mode and use "seek", but you will have some unwanted old data if the new data is shorter in length.

	# Write the file
	f = open(path, 'w')
	f.writelines(lines)
	f.close()

def editPythonFile(path,scale):
	f = open(path, 'r')    # pass an appropriate path of the required file
	lines = f.readlines()
	for i,line in enumerate(lines):
		if re.match(".*scale\s*=.*",line):
			lines[i] = "    scale = " + str(scale) + "\n"
	f.close()   # close the file and reopen in write mode to enable writing to file; you can also open in append mode and use "seek", but you will have some unwanted old data if the new data is shorter in length.

	# Write the file
	f = open(path, 'w')
	f.writelines(lines)
	f.close()

def editSubmitFile(path,scale):
	f = open(path, 'r')    # pass an appropriate path of the required file
	lines = f.readlines()
	for i,line in enumerate(lines):
		if re.match('.*#SBATCH -J ".*"   # job name.*',line):
			lines[i] = '#SBATCH -J "' + str(scale) + '_coating"   # job name\n'
	f.close()   # close the file and reopen in write mode to enable writing to file; you can also open in append mode and use "seek", but you will have some unwanted old data if the new data is shorter in length.

	# Write the file
	f = open(path, 'w')
	f.writelines(lines)
	f.close()


for size in sizesToRun:
	path = baseName + '/size_' + str(size).replace('.','_')
	create_directory(path)
	print("Created directory " + path)
	for filename in filesToCopy:
		dest = path + '/' + filename
		copyfile(filename,dest)
		print("Copied file " + filename + " to " + dest)

	# Do specific modifications here
	print("Modifying files in " + path)
	editInputFile(path + '/' + filesToCopy[0],str(size))
	editPythonFile(path + '/' + filesToCopy[1],str(size))
	editSubmitFile(path + '/' + filesToCopy[2],str(size))

	print("Changing directory to " + path)
	os.chdir(path)
	os.system('sbatch submit.sh')
	print("Returning to starting directory")
	os.chdir("../..")



