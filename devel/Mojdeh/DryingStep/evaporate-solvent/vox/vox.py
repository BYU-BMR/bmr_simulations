import re 
import numpy as np
from math import sqrt

cbd_type,active_type,solvent_type,wall_type = 1,2,3,4

class Atom:
	def __init__(self, x, y, z, Type):
		self.x = x
		self.y = y
		self.z = z
		self.Type = Type


readinFile = "start.lammpstrj"
outputFile = "voxel.data"
#handles
with open(readinFile,"r") as inFile: 
	with open(outputFile, "w") as outFile:
		# Find out how many frames there are
		NumFrame = 0
		for line in inFile:
			if re.match('ITEM: TIMESTEP.*',line):
				NumFrame += 1
		#print (NumFrame)

		# go to the last frame
		CurrentFrame = 0
		FoundOnLastFrame = False
		atomList = []
		position = {}
		fields = ['type', 'xs', 'ys', 'zs']

		# we need to tell it to go back to the beginning of the file
		inFile.seek(0) # go to the top of the file
		for line in inFile:
			if re.match('ITEM: TIMESTEP.*',line):
				CurrentFrame += 1
			if CurrentFrame == NumFrame:
				if re.match('ITEM: BOX.*',line):
					xedge = next(inFile).split()
					xmin = float(xedge[0])
					xmax = float(xedge[1])
					#print (xmin)
					#print (xmax)
					yedge = next(inFile).split()

					ymin = float(yedge [0])
					ymax = float(yedge [1])
					#print (ymin)
					#print (ymax)

					zedge = next(inFile).split()
					zmin = float(zedge [0])
					zmax = float(zedge [1])

					#print (zmin)
					#print (zmax)


			if CurrentFrame == NumFrame:
				if re.match('ITEM: ATOMS.*',line):
					#SAVE ALL THE LINES AFTER
					FoundOnLastFrame = True
					for i,item in enumerate(line.split()):
						for field in fields:
							if item == field:
								position[field] = i-2		

				elif FoundOnLastFrame:
					items = line.split()
					atomtype = int(items[position['type']])
					x = float(items[position['xs']])
					y = float(items[position['ys']])
					z = float(items[position['zs']])
					if atomtype in [active_type,cbd_type,solvent_type]:
						atom = Atom(x,y,z,atomtype)
						atomList.append(atom)


		# Find the limits of what we need to convrt to voxels
		'''xmin = float('inf')
		xmax = float('-inf')
		ymin = float('inf')
		ymax = float('-inf')
		zmin = float('inf')
		zmax = float('-inf')

		for atom in atomList:
			xmin = min(atom.x,xmin)
			ymin = min(atom.y,ymin)
			zmin = min(atom.z,zmin)

			xmax = max(atom.x,xmax)
			ymax = max(atom.y,ymax)
			zmax = max(atom.z,zmax)'''

		# Defining the size of the grid

		resolution = 1 #each voxel is 1 micrometer cubed, we need to decrease that late
		radius = 1
		def CheckVoxelType(x,y,z,Type):
			for atom in atomList: # we are checking if it's active
				if atom.Type == Type:
					dx = abs(x-atom.x) 
					dy = abs(y-atom.y) 
					dz = abs(z-atom.z)
					if dx < radius and  dy < radius and dz < radius:
						length = sqrt(dx**2+dy**2+dz**2)
						if length < radius:
							return True
			return False


		LinesToWrite = []
		atomCount = 0
		atomCount1 = 0
		atomCount2 = 0
		activeCounter = 0
		cbdCounter = 0
		solventCounter = 0
		poreCounter = 0

		for x in np.arange(xmin+0.5*resolution,xmax,resolution):
			for y in np.arange(ymin+0.5*resolution,ymax,resolution):
				for z in np.arange(zmin+0.5*resolution,zmax,resolution):
					if CheckVoxelType(x,y,z,active_type):
						#print ('found active voxel')
						atomline = "%d %f %f %f %f \n" % (active_type,x,y,z,resolution/2)
						LinesToWrite.append(atomline)
						atomCount += 1
						activeCounter += 1
					elif CheckVoxelType(x,y,z,cbd_type):
						#print ('found cbd voxel')
						atomline = "%d %f %f %f %f \n" % (cbd_type,x,y,z, resolution/2)
						LinesToWrite.append(atomline)
						atomCount1 += 1
						cbdCounter += 1
					elif CheckVoxelType(x,y,z,solvent_type):
						atomline = "%d %f %f %f %f \n" % (solvent_type,x,y,z, resolution/2)
						LinesToWrite.append(atomline)
						atomCount2 += 1
						solventCounter += 1
					else:
						poreCounter += 1
						
		print ("active voxels are", activeCounter, "voxels")
		print ("cbd voxels are", cbdCounter, "voxels")
		print ("solvent voxels are", solventCounter, "voxels")
		print ("pore voxels are", poreCounter, "voxels")
		#wrtie lammps trajectory file
		outFile.write("ITEM: TIMESTEP\n")
		outFile.write("0\n")
		outFile.write("ITEM: NUMBER OF ATOMS\n")
		outFile.write("%d\n" % (atomCount))
		outFile.write("ITEM: BOX BOUNDS pp pp pp\n")
		outFile.write("%f %f\n" % (xmin-radius, xmax+radius))
		outFile.write("%f %f\n" % (ymin-radius, ymax+radius))
		outFile.write("%f %f\n" % (zmin-radius, zmax+radius))
		outFile.write("ITEM: ATOMS type x y z radius\n")

		for line in LinesToWrite:
			outFile.write(line)
			






					





						

