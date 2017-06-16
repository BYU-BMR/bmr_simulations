import re 
import numpy as np
from math import sqrt

cbd_type,active_type,wall_type = 1,2,3

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
		fields = ['type', 'x', 'y', 'z']

		

		# Find the limits of what we need to convrt to voxels
		xmin = float('inf')
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
			zmax = max(atom.z,zmax)

		# Defining the size of the grid

		resolution = 0.4 #each voxel is 1 micrometer cubed, we need to decrease that late
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
		activeCounter = 0


						#print ('found pore voxel')
		print ("active voxels are", activeCounter, "voxels")
		#wrtie lammps trajectory file
		outFile.write("ITEM: TIMESTEP\n")
		outFile.write("0\n")
		outFile.write("ITEM: NUMBER OF ATOMS\n")
		






					





						

