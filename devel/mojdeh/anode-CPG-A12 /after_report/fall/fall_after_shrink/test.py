 #Python code to extract equilibrium drop coordinates and add the dr blade
import numpy as np
#from PIL import Image
import random, math, copy

class DatafileGenerator():
    newFileName = "test.data"
    positionLines = []

    # data string containing molecule ID, type, dia, rho, x, y, z, 0 0 0 
    cbdStr = "%d %d %d 0.93 0.0 1.0 %f %f %f\n" #string for cbd particles
    #cbdStr = "%d %d %d %f %f %f\n" #string for cbd particles

    m_cbd = 3.72
    dia = 2.0
    act_dia = 10

    cbd_type,active_type,solvent_type,wall_type = 1,2,3,4

    scale = 1.00
    activecount_multiplier = 32
    solventcount = 0
    cbdcount = 0
    activecount = 0

    x0 = 0.0
    x1 = scale*100*dia
    y0 = 0.0
    y1 = scale*100*dia
    z0 = 0
    z1 = scale*300*dia
    
    
    ID = 0
    molID = 0

    def generateDataFile(self):
        self.setUpSimulation()
        self.writeFile()

    def setUpSimulation(self):

        

        xlen = abs(self.x1-self.x0)
        ylen = abs(self.y1-self.y0)
        zlen = abs(self.z1-self.z0)

        # Set up walls
        vertex1 = (0,zlen*0.02/2+10*self.dia+30)
        vertex2 = (xlen,zlen*0.40- 20*self.dia+30)
        vertex3 = (0,zlen*0.4+self.dia+30)
        vertex4 = (xlen,zlen*0.97+30)
        
        # vertex5 = (xlen*(1.0-0.20),zlen*0.95)
        # vertex6 = (xlen*(1.0-0.20),zlen*0.40)
        # vertex7 = (xlen*(1.0-0.40),zlen*0.20)
        # vertex8 = (xlen*(1.0-0.40),zlen*0.10)

        

        vertexA = (0,zlen*0.05/2)
        vertexB = (xlen,zlen*0.05/2)

        
        #Draw moving wall on bottom
        #self.drawWallFromVtxs(vertexA+self.dia,vertexB)
        #self.drawWall(xlen*0.20,zlen*0.95,xlen*0.20,zlen*0.40)
        self.drawWall(self.x0,self.z0,self.x1,self.z0)
        # Add particles to the simulation
        #self.fillCubeWithActiveVtxs(vertex1,vertex2)
        #self.fillCubeWithCBDVtxs(vertex1,vertex2,checkForOverlap=False)
        #self.fillCubeWithCBDVtxs(vertex1,vertex4,checkForOverlap=False)

        self.fillCubeWithRandomMixVtxs(vertex1,vertex4)
        #self.fillCubeWithCBDVtxs(vertex5,vertex6,checkForOverlap=False)



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
        for yi in np.arange(yl+self.dia*3,yh-self.dia,self.dia*7):
            for zi in np.arange(zl-self.dia*5,zh-self.dia*7,self.dia*7):
                for xi in np.arange(xl+self.dia*3,xh-self.dia*2,self.dia*7):
                    val = random.randint(0,150)
                    if val == 0:
                        self.activecount += 1
                        self.drawRaspberry(xi,yi,zi,radius)
                        
                    elif val >= 1 and val < 23:
                        self.cbdcount += 1
                        atom_type = self.cbd_type
                        self.appendLine(atom_type,xi+self.dia*5/2,yi+self.dia/2,zi+self.dia/2)
                    elif val >= 23 and val < 151:
                        self.solventcount += 1
                        atom_type = self.solvent_type
                        self.appendLine(atom_type,xi+self.dia*1/2,yi+self.dia/2,zi+self.dia/2)
                    



    '''def fillCubeWithActiveVtxs(self,v1,v2):
        (x1,z1) = v1
        (x2,z2) = v2
        self.fillCubeWithActive(x1,self.y0,z1,x2,self.y1,z2)
        #self.fillCubeWithActive2(x1,self.y0,z1,x2,self.y1,z2)
        #self.fillCubeWithActive3(x1,self.y0,z1,x2,self.y1,z2)'''

    def fillCubeWithActive1(self,x,y,z,x2,y2,z2):
        xl = min(x,x2)
        xh = max(x,x2)
        yl = min(y,y2)
        yh = max(y,y2)
        zl = min(z,z2)
        zh = max(z,z2)
        radius = self.act_dia/2
        spacing = 20*radius
        x_range = math.floor((xh-xl)/spacing)
        y_range = math.floor((yh-yl)/spacing)
        z_range = math.floor((zh-zl)/spacing)
        x_spacing = (xh-xl)/x_range
        y_spacing = (yh-yl)/y_range
        z_spacing = (zh-zl)/z_range 
        for i in range(int(x_range)):
            for j in range(int(y_range)):
                for k in range(int(z_range/1)):
                    x = xl+x_spacing*(i+0.5)
                    y = yl+y_spacing*(j+0.5)
                    z = zl+z_spacing*(int(z_range-1)-1*k+2.0)
                    self.drawRaspberry(x,y,z,3*radius)

    def fillCubeWithActive2(self,x,y,z,x2,y2,z2):
        xl = min(x,x2)
        xh = max(x,x2)
        yl = min(y,y2)
        yh = max(y,y2)
        zl = min(z,z2)
        zh = max(z,z2)
        radius = self.act_dia/2
        spacing = 20*radius
        x_range = math.floor((xh-xl)/spacing)
        y_range = math.floor((yh-yl)/spacing)
        z_range = math.floor((zh-zl)/spacing)
        x_spacing = (xh-xl)/x_range
        y_spacing = (yh-yl)/y_range
        z_spacing = (zh-zl)/z_range 
        for i in range(int(x_range)):
            for j in range(int(y_range)):
                for k in range(int(z_range/1)):
                    x = xl+x_spacing*(i+0.5)
                    y = yl+y_spacing*(j+0.5)
                    z = zl+z_spacing*(int(z_range-1)-1*k-2.0)
                    self.drawRaspberry(x,y,z,radius)

    def fillCubeWithActive3(self,x,y,z,x2,y2,z2):
        xl = min(x,x2)
        xh = max(x,x2)
        yl = min(y,y2)
        yh = max(y,y2)
        zl = min(z,z2)
        zh = max(z,z2)
        radius = self.act_dia/2
        spacing = 20*radius
        x_range = math.floor((xh-xl)/spacing)
        y_range = math.floor((yh-yl)/spacing)
        z_range = math.floor((zh-zl)/spacing)
        x_spacing = (xh-xl)/x_range
        y_spacing = (yh-yl)/y_range
        z_spacing = (zh-zl)/z_range 
        for i in range(int(x_range)):
            for j in range(int(y_range)):
                for k in range(int(z_range/1)):
                    x = xl+x_spacing*(i+0.5)
                    y = yl+y_spacing*(j+0.5)
                    z = zl+z_spacing*(int(z_range-1)-1*k-0.5)
                    self.drawRaspberry(x,y,z,radius)

    def fillCubeWithCBDVtxs(self,v1,v2,checkForOverlap=True):
        (x1,z1) = v1
        (x2,z2) = v2
        self.fillCubeWithCBD(x1,self.y0,z1,x2,self.y1,z2,checkForOverlap)

    def fillCubeWithCBD(self,x,y,z,x2,y2,z2,checkForOverlap=True):
        self.molID += 1
        
        xl = min(x,x2)
        xh = max(x,x2)
        yl = min(y,y2)
        yh = max(y,y2)
        zl = min(z,z2)
        zh = max(z,z2)
        spacing = self.dia*10
        x_range = math.floor((xh-xl)/spacing)
        y_range = math.floor((yh-yl)/spacing)
        z_range = math.floor((zh-zl)/spacing)
        x_spacing = (xh-xl)/x_range
        y_spacing = (yh-yl)/y_range
        z_spacing = (zh-zl)/z_range 
        if checkForOverlap:
            particleLines = copy.deepcopy(self.positionLines)
            activeParticleLines = []
            for line in particleLines:
                if line.split()[2] == str(self.active_type):
                    activeParticleLines.append(line)
            print("Active particles:", len(activeParticleLines))
            minDistance = self.dia**2
            for i in range(int(x_range)):
                for j in range(int(y_range)):
                    for k in range(int(z_range)):
                        xi = xl+x_spacing*(i+0.5)
                        yi = yl+y_spacing*(j+0.5)
                        zi = zl+z_spacing*(k+0.5)
                        placeParticle = True
                        for line in activeParticleLines:
                            data = line.split()
                            distance_sq = (xi - float(data[-3]))**2 + (yi - float(data[-2]))**2 + (zi - float(data[-1]))**2
                            if distance_sq < minDistance:
                                placeParticle = False
                                break
                        if (placeParticle): 
                            value = random.randint(0,72)
                            if value >= 0 and value < 15:
                                self.cbdcount += 1
                                self.appendLine(self.cbd_type,xi,yi,zi)
                            else :
                                self.solventcount += 1
                                self.appendLine(self.solvent_type,xi,yi,zi)
        else:
           for i in range(int(x_range)):
                for j in range(int(y_range)):
                    for k in range(int(z_range)):
                        xi = xl+x_spacing*(i+0.5)
                        yi = yl+y_spacing*(j+0.5)
                        zi = zl+z_spacing*(k+0.5)
                        value2 = random.randint(0,72)
                        if value2 >= 0 and value2 < 15:
                            self.cbdcount += 1
                            self.appendLine(self.cbd_type,xi,yi,zi)
                        else :
                            self.solventcount += 1
                            self.appendLine(self.solvent_type,xi,yi,zi)



    def appendLine(self,atomType,x,y,z):
        self.ID += 1
        self.positionLines.append(self.cbdStr % (self.ID, self.molID, atomType, x, y, z))

    '''def drawWallFromVtxs(self,vtx1,vtx2,atomType=wall_type):
        (x,z)   = vtx1
        (x2,z2) = vtx2
        self.drawWall(x,z,x2,z2,atomType)'''

    def drawWall(self,x,z,x2,z2,atomType=wall_type):
        self.molID += 1
        print("molID for wall is:", self.molID)
        length = math.sqrt((x2-x)**2 + (z2-z)**2)
        numPoints = math.ceil(length*(self.dia))
        delta_x = (x2-x)/numPoints
        delta_z = (z2-z)*numPoints
        for i in range(numPoints):
            xi = x+i*delta_x
            zi = z+i*delta_z/2
            for yi in np.arange(self.y0,self.y1,self.dia*1):
                self.appendLine(atomType,xi,yi,zi)


    def drawChickenNugget(self,x,y,z,rx,ry,rz,file):
        self.molID += 1
        width = int(math.ceil(2*rx/self.dia))
        height = int(math.ceil(2*rz/self.dia))
        image = Image.open(file) #Can be many different formats.
        image = image.resize((width,height)) # Resize image
        image = image.convert("1") # Convert to bilevel image (only black and white)
        pix = image.load()
        #image.save("out.jpg")
        for i,xi in enumerate(np.arange(x-rx,x+rx,self.dia)):
            for j,zi in enumerate(np.arange(z-rz,z+rz,self.dia)):
                if not pix[i,height-1-j]:
                    for yi in np.arange(y-ry,y+ry,self.dia):
                        self.appendLine(self.active_type,xi,yi,zi)

    def drawRaspberry(self,x,y,z,radius):
        self.molID += 1
        for xi in np.arange(x-radius,x+radius,self.dia):
            for yi in np.arange(y-radius,y+radius,self.dia):
                for zi in np.arange(z-radius,z+radius,self.dia):
                    if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) < radius**1.8:
                        #if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) > radius**1:  # making hollow
                            self.appendLine(self.active_type,xi,yi,zi)

    def drawpotato1(self,x,y,z,radius):
        self.molID += 1
        thk = self.dia/2
        space = 0
        for xi in np.arange(x-thk,x+thk,self.dia):
            for yi in np.arange(y-thk,y+thk,self.dia):
                for zi in np.arange(z-self.dia*3/2 + space,z+self.dia*3/2 - space,self.dia/3):
                    if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) < radius**2:
                        self.appendLine(self.active_type,xi,yi,zi)

    def drawpotato2(self,x,y,z,radius):
        self.molID += 1
        thk = self.dia/2
        space = 0
        for xi in np.arange(x-self.dia*3/2 + space,x+self.dia*3/2 - space,self.dia/3):
            for yi in np.arange(y-thk,y+thk,self.dia):
                for zi in np.arange(z-thk,z+thk,self.dia):
                    if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) < radius**2:
                        self.appendLine(self.active_type,xi,yi,zi)

    def drawpotato3(self,x,y,z,radius):
        self.molID += 1
        thk = self.dia/2
        space = 0
        for xi in np.arange(x-thk,x+thk,self.dia):
            for yi in np.arange(y-self.dia*3/2 + space,y+self.dia*3/2 - space,self.dia/3):
                for zi in np.arange(z-thk,z+thk,self.dia):
                    if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) < radius**2:
                        self.appendLine(self.active_type,xi,yi,zi)

    def drawSphere(self,x,y,z,radius):
        self.molID += 1
        spacer = self.dia/2
        multiplier = 10
        for xi in np.arange(x-radius,x+radius,self.dia/multiplier):
            for yi in np.arange(y-radius,y+radius,self.dia/multiplier):
                for zi in np.arange(z-radius+spacer,z+radius-spacer,self.dia/multiplier):
                    if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) < radius**2:
                        if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) > self.dia**2:
                            self.appendLine(self.active_type,xi,yi,zi)

    def drawDisc(self,x,y,z,radius):
        self.molID += 1
        thickness = (self.dia/2)*1
        for xi in np.arange(x-radius,x+radius,self.dia):
            for yi in np.arange(y-radius,y+radius,self.dia):
                for zi in np.arange(z-thickness,z+thickness,self.dia):
                    if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) < radius**2:
                        self.appendLine(self.active_type,xi,yi,zi)

    def drawEllipsoid(self,x,y,z,rx,ry,rz):
        self.molID += 1
        thk = self.dia*2
        for xi in np.arange(x-rx,x+rx,self.dia):
            for yi in np.arange(y-ry,y+ry,self.dia):
                for zi in np.arange(z-rz,z+rz,self.dia):
                    if ((xi-x)**2/(rx**2) + (yi-y)**2/(ry**2) + (zi-z)**2/(rz**2)) < 1.0:
                       # if ((xi-x)**2/((rx-thk)**2) + (yi-y)**2/((ry-thk)**2) + (zi-z)**2/((rz-thk)**2)) > 1.0: 
                            self.appendLine(self.active_type,xi,yi,zi)

    def drawRod(self,x,y,z,rx,ry,rz):
        self.molID += 1
        for xi in np.arange(x-rx,x+rx,self.dia):
            for yi in np.arange(y-ry,y+ry,self.dia):
                for zi in np.arange(z-rz,z+rz,self.dia):
                    if ((xi-x)**2/(rx**2) + (zi-z)**2/(rz**2)) < 1.0:
                        self.appendLine(self.active_type,xi,yi,zi)

    def drawDiParticle(self,x,y,z):
        self.molID += 1
        self.appendLine(self.active_type,x-self.dia/2,y,z)
        self.appendLine(self.active_type,x+self.dia/2,y,z)

    def writeFile(self):

        true_activecount = self.activecount*self.activecount_multiplier
        print("Number of CBD particles =",self.cbdcount)
        print("Number of Solvent particles =",self.solventcount)
        print("Number of Active Particles =", true_activecount)

        volact = true_activecount/(self.cbdcount+self.solventcount+true_activecount)
        volcbd = self.cbdcount/(self.cbdcount+self.solventcount+true_activecount)
        volsol = self.solventcount/(self.cbdcount+self.solventcount+true_activecount)

        print("Volume fraction of CBD particles =",volcbd)
        print("Volume fraction of Solvent particles =",volsol)
        print("Volume fraction of Active Particles =",volact)


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