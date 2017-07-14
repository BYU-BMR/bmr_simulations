import re 
import numpy as np
import math as math

cbd_type,active_type,wall_type = 1,2,4
cbd_type_new,active_type_new,pore = 0, 140, 255

porosity_of_cbd = 0.55

slice_horizontal = True

class Atom:
	def __init__(self, x, y, z, Type):
		self.x = x
		self.y = y
		self.z = z
		self.Type = Type


readinFile = "less_repel_6.lammpstrj"
outputFile = "slice_of_voxel.data"
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
		zlist = []
		position = {}
		fields = ['type', 'x', 'y', 'z']

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
					x = float(items[position['x']])
					y = float(items[position['y']])
					z = float(items[position['z']])
					if atomtype in [active_type,cbd_type]:
						zlist.append(z)
						atom = Atom(x,y,z,atomtype)
						atomList.append(atom)

		
		# Find the limits of what we need to convrt to voxels
		'''xmin = float('inf')
		xmax = float('-inf')
		ymin = float('inf')
		ymax = float('-inf')
		zmin = float('inf')
		zmax = float('-inf')'''
	
		'''for atom in atomList:
			xmin = min(atom.x,xmin)
			ymin = min(atom.y,ymin)
			zmin = min(atom.z,zmin)

			xmax = max(atom.x,xmax)
			ymax = max(atom.y,ymax)
			zmax = max(atom.z,zmax) + 1'''

	

		
		# Defining the size of the grid
		print ("old zmax = ",zmax)
		print ("old zmin = ",zmin)
		resolution = 0.4 #each voxel is 1 micrometer cubed, we need to decrease that late
		Aradius = 1.0
		Aradius2 = Aradius*Aradius
		Cradius = 1.00
		Cradius2 = Cradius*Cradius
		zmax = max(zlist) + Aradius
		zmin = min(zlist) - Aradius
		print ("new zmax = ",zmax)
		print ("new zmin = ",zmin)

		xlen = xmax - xmin
		ylen = ymax - ymin
		zlen = zmax - zmin

		'''NX = math.floor((xlen - 0.2)/0.4)
		NY = math.floor((ylen - 0.2)/0.4)	
		NZ = math.floor((zlen - 0.2)/0.4)

		print ("NX = ",NX)
		print ("NY = ",NY)
		print ("NZ = ",NZ)

		outFile.write("NX = %f NY = %f NZ = %f\n" % (NX,NY,NZ))'''


		def CheckVoxelType(x,y,z,Type,radius,radiussqrd):
			for atom in atomList: # we are checking if it's active
				if atom.Type == Type:
					
					dx = (x-atom.x)
					if  dx > xlen*0.5:  # the following if statements allow this simulation to adhere to "minimum image convention"... google it :)
						dx = dx - xlen
					elif dx <= -xlen*0.5:
						dx = dx + xlen
					if dx*dx < radiussqrd:

						dy = (y-atom.y)
						if  dy > ylen*0.5:
							dy = dy - ylen
						elif dy <= -ylen*0.5:
							dy = dy + ylen 
						if dy*dy < radiussqrd:

							dz = (z-atom.z)
							if  dz > zlen*0.5:
								dz = dz - zlen
							elif dz <= -zlen*0.5:
								dz = dz + zlen
							if dz*dz < radiussqrd:

								length = (dx*dx+dy*dy+dz*dz)
								if length < radiussqrd:
									return True
			return False


		LinesToWrite = []
		atomCount = 0
		activeCounter = 0
		carbonCounter = 0
		poreCounter = 0

		if slice_horizontal == True:
			zmax = zmax - 0.4*zlen
			zmin = zmin + 0.4*zlen


		for x in np.arange(xmin+0.5*resolution,xmax,resolution):
			for y in np.arange(ymin+0.5*resolution,ymax,resolution):
				for z in np.arange(zmin+0.5*resolution,zmax,resolution):
					if CheckVoxelType(x,y,z,active_type,Aradius,Aradius2):
						#print ('found active voxel')
						#atomline = "%d %f %f %f %f \n" % (active_type_new,x,y,z,resolution/2) #use this line to see in Ovito
						atomline = "%f %f %f %d \n" % (x,y,z,active_type_new)
						LinesToWrite.append(atomline)
						atomCount += 1
						activeCounter += 1
					elif CheckVoxelType(x,y,z,cbd_type,Cradius,Cradius2):
						#print ('found cbd voxel')
						#atomline = "%d %f %f %f %f \n" % (cbd_type_new,x,y,z, resolution/2) #use this line to see in Ovito
						atomline = "%f %f %f %d \n" % (x,y,z,cbd_type_new)
						LinesToWrite.append(atomline)
						atomCount += 1
						carbonCounter += 1
					else:
						poreCounter += 1
						#pass # use just this line to see in Ovito
						#atomline = "%d %f %f %f %f \n" % (pour,x,y,z, resolution/2)
						atomCount += 1
						atomline = "%f %f %f %d \n" % (x,y,z,pore)
						LinesToWrite.append(atomline)			

		print ("There are", activeCounter, "active voxels")
		print ("There are", carbonCounter, "carbon voxels")
		print ("There are", poreCounter, "pore voxels")

		volact = activeCounter/(activeCounter+carbonCounter+poreCounter)
		volcbd = carbonCounter/(activeCounter+carbonCounter+poreCounter)
		volpore = poreCounter/(activeCounter+carbonCounter+poreCounter)
		Porosity = (poreCounter + porosity_of_cbd*carbonCounter)/(atomCount) 

		print ("Volume Fraction of Active = ",volact)
		print ("Volume Fraction of Carbon = ",volcbd)
		print ("Volume Fraction of Pore = ",volpore)
		print ("Porosity = ",Porosity)


		print ("Writing new file, ",outFile)
		#wrtie lammps trajectory file
		outFile.write("ITEM: TIMESTEP\n")
		outFile.write("0\n")
		outFile.write("ITEM: NUMBER OF ATOMS\n")
		outFile.write("%d\n" % (atomCount))
		outFile.write("ITEM: BOX BOUNDS pp pp pp\n")
		outFile.write("%f %f\n" % (xmin, xmax))
		outFile.write("%f %f\n" % (ymin, ymax))
		outFile.write("%f %f\n" % (zmin, zmax))
		#outFile.write("ITEM: ATOMS type x y z radius\n") #use this line to see in Ovito
		outFile.write("ITEM: ATOMS x y z type\n")
		for line in LinesToWrite:
			outFile.write(line)
			






					





						

