# Python code to extract equilibrium drop coordinates and add the dr blade
import numpy as np
#from PIL import Image
import random, math, copy

class DatafileGenerator():
    newFileName = "test_tall1.data"
    positionLines = []

    # data string containing molecule ID, type, dia, rho, x, y, z, 0 0 0 
    cbdStr = "%d %d %d 0.93 0.0 1.0 %f %f %f\n" #string for cbd particles
    solStr = "%d %d %d 1.028 0.0 1.0 %f %f %f\n" #string for cbd particles
    actStr = "%d %d %d 4.79 0.0 1.0 %f %f %f\n" #string for cbd particles
    #cbdStr = "%d %d %d %f %f %f\n" #string for cbd particles
    #solStr = "%d %d %d %f %f %f\n" #string for cbd particles
    #actStr = "%d %d %d %f %f %f\n" #string for cbd particles

    m_cbd = .9739
    m_sol = 1.0765
    m_act = 30.10
    dia = 2.0
    act_dia = 8

    cbd_type,solvent_type,active_type = 1,2,3

    scale = 1.00

    solventcount = 0
    cbdcount = 0

    x0 = 0.0
    x1 = scale*10.0*dia
    y0 = 0.0
    y1 = scale*10.0*dia
    z0 = 0
    z1 = scale*400.0*dia

    ID = 0
    molID = 0

    def generateDataFile(self):
        self.setUpSimulation()
        self.writeFile()

    def setUpSimulation(self):
        xlen = abs(self.x1-self.x0)
        ylen = abs(self.y1-self.y0)
        zlen = abs(self.z1-self.z0)
              
        vertex1 = (0,zlen*.1)
        vertex2 = (xlen,zlen*.9)

        self.fillCubeWithRandomMixVtxs(vertex1,vertex2)

    def fillCubeWithRandomMixVtxs(self,v1,v2):
        (x1,z1) = v1
        (x2,z2) = v2
        self.fillCubeWithRandomMix(x1,self.y0,z1,x2,self.y1,z2)

    def fillCubeWithRandomMix(self,x,y,z,x2,y2,z2):
        xl = min(x,x2)
        xh = max(x,x2)
        yl = min(y,y2)
        yh = max(y,y2)
        zl = min(z,z2)
        zh = max(z,z2)
        radius = self.act_dia/2
        for yi in np.arange(yl,yh-self.dia,self.dia*2):
            for zi in np.arange(zl,zh-self.dia,self.dia*2):
                for xi in np.arange(xl,xh-self.dia*2,self.dia*2):
                    val = random.randint(0,83)
                    if val >= 16 and val < 18:
                        self.molID += 1
                        atom_type = self.active_type
                        self.appendLine(atom_type,xi+self.dia*1/2,yi+self.dia/2,zi+self.dia/2)

                    elif val >= 18 and val < 32:
                        atom_type = self.cbd_type
                        self.appendLine(atom_type,xi+self.dia*1/2,yi+self.dia/2,zi+self.dia/2)
                    else:
                        atom_type = self.solvent_type
                        self.appendLine(atom_type,xi+self.dia*1/2,yi+self.dia/2,zi+self.dia/2)


    def appendLine(self,atomType,x,y,z):
        self.ID += 1
        if atomType == self.cbd_type:
            self.positionLines.append(self.cbdStr % (self.ID, self.molID, atomType, x, y, z))
        elif atomType == self.solvent_type:
            self.positionLines.append(self.solStr % (self.ID, self.molID, atomType, x, y, z))
        else:
            self.positionLines.append(self.actStr % (self.ID, self.molID, atomType, x, y, z))

    
    def writeFile(self):
        print("Number of CBD particles =",self.cbdcount)
        print("Number of Solvent particles =",self.solventcount)
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
            outFile.write("2 %f\n" % self.m_sol)
            outFile.write("3 %f\n" % self.m_act)
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