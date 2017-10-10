import numpy as np 
import math, random, copy, re
from shutil import copyfile 

class DatafileGenerator():
	newFileName = "nmp_collapse_overcompensation"
	positionLines = []
	velocityLines = []
	xline = []

	# atom data string containing atom_id atom_type rho ev cp x y z
	#atomstr = "%d %d 1030 0.0 1.0 %f %f %f\n"
	atomstr = "%d %d %f %f %f\n" # use this string instead if previewing in ovito

	solvent_type,wall_type = 1,2

	m_solvent = 0.19418
	m_wall = 0.1

	dia = 0.1

	ID = 0

	xlo = 0
	xhi = 4.0010000000000003e+00
	ylo = 0
	yhi = 8.0009999999999994e+00
	zlo = 0
	zhi = 4.0010000000000003e+00


	def generateDataFile(self):
		self.setUpSimulation()
		self.writeFile()
		print ("generateDataFile is working")

	def setUpSimulation(self):
		xlen = abs(self.xlo-self.xhi)
		ylen = abs(self.yhi-self.ylo)
		zlen = abs(self.zhi-self.zlo)

		radius = self.dia/2
		wall_radius = self.dia/3

		vertex1 = (0,0)
		vertex2 = (0.4*xlen,0.4*zlen)

		vertexa = (0,ylen-radius/2,0)
		vertexb = (xlen,ylen,zlen)

		print ("setUpSimulation is working")

		self.FccVertexes(vertex1,vertex2,self.yhi*.5,1)
		self.DrawWallVtx(vertexa,vertexb,wall_radius,'XZ')

	def appendLine(self,atomType,x,y,z):
		self.ID += 1
		self.positionLines.append(self.atomstr % (self.ID, atomType, x, y, z))


	def FccVertexes(self,vertexa,vertexb,yhi,atomType):
		xlo = vertexa[0]
		xhi = vertexb[0]
		zlo = vertexa[1]
		zhi = vertexb[1]
		ylo = self.ylo

		print ("FccVertexes is working")

		self.FCC(xlo,xhi,ylo,yhi,zlo,zhi,atomType)


	def FCC(self,xlo,xhi,ylo,yhi,zlo,zhi,atomType):

		for k in np.arange(zlo-zhi,5*zhi,self.dia/2):
			for j in np.arange(ylo-yhi,5*yhi,self.dia/2):
				for i in np.arange(xlo-xhi,5*xhi,self.dia/2):
					x = (2*i + (( j + k )%2))
					y = (np.sqrt(3)*( j + (1/3)*(k%2)))
					z = ((2*np.sqrt(6)/3)*k)
					if z > zhi or x > xhi or y > yhi:
						break
					if z > zlo:
						if x > xlo:
							if y > ylo:
								self.appendLine(atomType,x,y,z)
		print ("FCC is working")

	def DrawWallVtx(self,Vertex1,Vertex2,radius,plane):
		xlo = Vertex1[0]
		xhi = Vertex2[0]
		ylo = Vertex1[1]
		yhi = Vertex2[1]
		zlo = Vertex1[2]
		zhi = Vertex2[2]

		if plane == 'XY':
			self.DrawWallXY(xlo,xhi,ylo,yhi,zlo,zhi,radius)
		if plane == 'XZ':
			self.DrawWallXZ(xlo,xhi,ylo,yhi,zlo,zhi,radius)
		if plane == 'YZ':
			self.DrawWallYZ(xlo,xhi,ylo,yhi,zlo,zhi,radius)

	def DrawWallXZ(self,xlo,xhi,ylo,yhi,zlo,zhi,radius):
		atomType = self.wall_type
		beta = 0
		dia = 2*radius
		for k in np.arange(zlo,zhi,radius*np.sqrt(3)):
			beta += 1
			for i in np.arange(xlo,xhi,dia):
				alpha = int((k-zlo)/(radius*np.sqrt(3)))
				x = (-1)**beta*(radius/2) + i
				y = ylo + (yhi-ylo)/2
				z = k
				self.appendLine(atomType,x,y,z)

	def writeFile(self):
		print("writeFile is working")
		print("Writing new file: %s" % self.newFileName)

		xlo = self.xlo
		xhi = self.xhi
		ylo = self.ylo
		yhi = self.yhi
		zlo = self.zlo
		zhi = self.zhi

		with open(self.newFileName,"w") as outFile:

			outFile.write("\n")
			outFile.write("%d atoms\n" % self.ID)
			outFile.write("\n")
			outFile.write("%d atom types\n" % 2)
			outFile.write("\n")
			outFile.write("%f %f xlo xhi\n" % (xlo,xhi))   
			outFile.write("%f %f ylo yhi\n" % (ylo,yhi))
			outFile.write("%f %f zlo zhi\n" % (zlo,zhi))

			outFile.write("\n")
			outFile.write("Masses\n")
			outFile.write("\n")
			outFile.write("1 %f\n" % self.m_solvent)
			outFile.write("2 %f\n" % self.m_wall)

			outFile.write("\n")
			outFile.write("Atoms\n")
			outFile.write("\n")

			for line in self.positionLines:
				outFile.write(line)
			for line in self.xline:
				outFile.write(line)

if __name__ == "__main__":
	generator = DatafileGenerator()
	generator.generateDataFile()
