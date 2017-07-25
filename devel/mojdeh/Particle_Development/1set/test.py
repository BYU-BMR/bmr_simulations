# Python code to extract equilibrium drop coordinates and add the dr blade
import numpy as np
#from PIL import Image
import random, math, copy

class DatafileGenerator():
    newFileName = "test.data"
    positionLines = []

    # data string containing molecule ID, type, dia, rho, x, y, z, 0 0 0 
    cbdStr = "%d %d %d 0.93 0.0 1.0 %f %f %f\n" #string for cbd particles
    #actStr = "%d %d %d 0.93 0.0 1.0 %f %f %f\n" #string for active particles
    #solStr = "%d %d %d 0.93 0.0 1.0 %f %f %f\n" #string for solvent particles
    #wallStr = "%d %d %d 0.93 0.0 1.0 %f %f %f\n" #string for wall particles
    #cbdStr = "%d %d %d %f %f %f\n" #string for cbd particles
    #solStr = "%d %d %d %f %f %f\n" #string for solvent particles
    #actStr = "%d %d %d %f %f %f\n" #string for active particles
    #wallStr = "%d %d %d %f %f %f\n" #string for wall particles

    m_cbd = 3.896
    m_act = 1284
    m_sol = 4.306
    #m_wall = 20.94

    cbd_dia = 2.0
    act_dia = 2.0
    sol_dia = 2.0
    #wall_dia = 2.0

    cbd_type,active_type,solvent_type = 1,2,3

    scale = 1.00

    solventcount = 0
    cbdcount = 0

    x0 = 0.0
    x1 = scale*40.0*cbd_dia
    y0 = 0.0
    y1 = scale*20.0*cbd_dia
    z0 = 0.0
    z1 = scale*20.0*cbd_dia

    ID = 0
    molID = 0

    def generateDataFile(self):
        self.setUpSimulation()
        self.writeFile()

    def setUpSimulation(self):
        xlen = abs(self.x1-self.x0)
        ylen = abs(self.y1-self.y0)
        zlen = abs(self.z1-self.z0)
              
        self.placeparticles1(xlen,ylen,zlen,self.cbd_dia,self.cbd_type,self.m_cbd)
        #self.placeparticles2(xlen,ylen,zlen,self.act_dia,self.active_type,self.m_act)
        #self.placeparticles3(xlen,ylen,zlen,self.sol_dia,self.solvent_type,self.m_sol)

        

        
        #self.placeparticlesclosely(xlen,ylen,zlen,self.sol_dia,self.solvent_type)
              
    def placeparticles1(self,xlen,ylen,zlen,dia,Type,mass):
        x = (1/3)*xlen - 10
        y = .5*ylen
        z = .5*zlen - 15
        self.molID += 1
        self.appendLine(self.cbd_type,x,y,z)

        x = (2/3)*xlen - 10
        self.molID += 1
        self.appendLine(self.cbd_type,x,y,z)

        x = (3/3)*xlen - 10
        self.molID += 1
        self.appendLine(self.cbd_type,x,y,z)

    '''def placeparticles2(self,xlen,ylen,zlen,dia,Type,mass):
        x = (1/3)*xlen - 10
        y = 0.5*ylen
        z = 0.75*zlen - 14
        self.molID += 1
        self.appendLine(self.active_type,x,y,z)

        x = (2/3)*xlen - 10
        self.molID += 1
        self.appendLine(self.active_type,x,y,z)

        x = (3/3)*xlen - 10
        self.molID += 1
        self.appendLine(self.active_type,x,y,z)

    def placeparticles3(self,xlen,ylen,zlen,dia,Type,mass):
        x = (1/3)*xlen - 10
        y = .5*ylen
        z = 1*zlen - 13
        self.molID += 1
        self.appendLine(self.solvent_type,x,y,z)

        x = (2/3)*xlen - 10
        self.molID += 1
        self.appendLine(self.solvent_type,x,y,z)

        x = (3/3)*xlen - 10
        self.molID += 1
        self.appendLine(self.solvent_type,x,y,z) '''

    def placeparticlesclosely(self,xlen,ylen,zlen,dia,Type):
        x = .5*xlen - .5*dia
        y = .5*ylen
        z = .5*zlen
        self.molID += 1
        self.appendLine(Type,x,y,z)

        x = .5*xlen
        self.molID += 1
        self.appendLine(Type,x,y,z)

        x = .5*xlen + dia
        self.molID += 1
        self.appendLine(Type,x,y,z)


    def appendLine(self,atomType,x,y,z):
        self.ID += 1
        if atomType == self.cbd_type:
            self.positionLines.append(self.cbdStr % (self.ID, self.molID, atomType, x, y, z))
        elif atomType == self.solvent_type:
            self.positionLines.append(self.actStr % (self.ID, self.molID, atomType, x, y, z))
        else: #atomType == self.active_type:
            self.positionLines.append(self.solStr % (self.ID, self.molID, atomType, x, y, z))
        #else:
            #self.positionLines.append(self.wallStr % (self.ID, self.molID, atomType, x, y, z))

    
    def writeFile(self):
        #print("Number of CBD particles =",self.cbdcount)
        #print("Number of Solvent particles =",self.solventcount)
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
            outFile.write("2 %f\n" % self.m_act)
            outFile.write("3 %f\n" % self.m_sol)
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