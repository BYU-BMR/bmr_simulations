import numpy as np
#from PIL import Image
import random, math, copy

class DatafileGenerator():
    newFileName = "test_height.data"
    positionLines = []

    # data string containing molecule ID, type, dia, rho, x, y, z, 0 0 0 
    cbdStr = "%d %d %d 0.93 0.0 1.0 %f %f %f\n" #string for cbd particles
    #cbdStr = "%d %d %d %f %f %f\n" #string for cbd particles

    m_cbd = 3.72
  
    dia = 2.0

    cbd_type,solvent_type,active_type,wall_type = 1,2,3,4

    scale = 1.00

    x0 = 0.0
    x1 = scale*10.0*dia
    y0 = 0.0
    y1 = scale*10.0*dia
    z0 = 0
    z1 = scale*10.0*dia

    ID = 0
    molID = 0

    def generateDataFile(self):
        self.setUpSimulation()
        self.writeFile()

    def setUpSimulation(self):
        xlen = abs(self.x1-self.x0)
        ylen = abs(self.y1-self.y0)
        zlen = abs(self.z1-self.z0)
              
        self.placeparticles(xlen,ylen,zlen)
              
    def placeparticles(self,xlen,ylen,zlen):
        x = 2
        y = 1
        z = 1
        self.molID += 1
        self.appendLine(self.cbd_type,x,y,z)

        x = 20
        self.molID += 1
        self.appendLine(self.active_type,x,y,z)

        x = 1
        z = 10
        self.molID += 1
        self.appendLine(self.cbd_type,x,y,z)

        x = 20
        y = 20
        z = 10
        self.molID += 1
        self.appendLine(self.active_type,x,y,z)


    def appendLine(self,atomType,x,y,z):
        self.ID += 1
        self.positionLines.append(self.cbdStr % (self.ID, self.molID, atomType, x, y, z))

    def writeFile(self):
        print("Writing new file: %s" % self.newFileName)

        xmin = self.x0
        xmax = self.x1
        ymin = self.y0
        ymax = self.y1
        zmin = self.z0
        zmax = self.z1

        

        with open(self.newFileName,"w") as outFile:

            outFile.write("\n")
            outFile.write("%d atoms\n" % self.ID)
            outFile.write("\n")
            outFile.write("%d atom types\n" % 8)
            outFile.write("\n")
            outFile.write("%f %f xlo xhi\n" % (xmin,xmax))   
            outFile.write("%f %f ylo yhi\n" % (ymin,ymax))
            outFile.write("%f %f zlo zhi\n" % (zmin,zmax))

            outFile.write("\n")
            outFile.write("Masses\n")
            outFile.write("\n")
            outFile.write("1 %f\n" % self.m_cbd)
            outFile.write("2 %f\n" % self.m_cbd)
            outFile.write("3 %f\n" % self.m_cbd)
            outFile.write("4 %f\n" % self.m_cbd)
            outFile.write("5 %f\n" % self.m_cbd)
            outFile.write("6 %f\n" % self.m_cbd)
            outFile.write("7 %f\n" % self.m_cbd)
            outFile.write("8 %f\n" % self.m_cbd)

            outFile.write("\n")
            outFile.write("Atoms\n")
            outFile.write("\n")

            for line in self.positionLines:
                outFile.write(line)
    
if __name__ == "__main__":
    generator = DatafileGenerator()
    generator.generateDataFile()