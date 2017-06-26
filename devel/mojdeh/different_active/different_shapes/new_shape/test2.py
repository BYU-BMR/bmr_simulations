# Python code to extract equilibrium drop coordinates and add the dr blade
import numpy as np

import random, math, copy

class DatafileGenerator():
    newFileName = "test_box2.data"
    positionLines = []

    # data string containing molecule ID, type, dia, rho, x, y, z, 0 0 0 
    #cbdStr = "%d %d %d 0.93 0.0 1.0 %f %f %f\n" #string for cbd particles
    cbdStr = "%d %d %d %f %f %f\n" #string for cbd particles

    m_cbd = 3.72
    dia = 1.14

    cbd_type,active_type,wall_type = 1,2,3

    scale = 1.0

    x0 = 0.0
    x1 = scale*7*50.0*dia
    y0 = 0.0
    y1 = scale*5.0*dia
    z0 = 0
    z1 = scale*3*50.0*dia

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
        # vertex1 = (xlen*0.20,zlen*0.95)
        # vertex2 = (xlen*0.20,zlen*0.40)
        # vertex3 = (xlen*0.40,zlen*0.20)
        # vertex4 = (xlen*0.40,zlen*0.10)
        
        # vertex5 = (xlen*(1.0-0.20),zlen*0.95)
        # vertex6 = (xlen*(1.0-0.20),zlen*0.40)
        # vertex7 = (xlen*(1.0-0.40),zlen*0.20)
        # vertex8 = (xlen*(1.0-0.40),zlen*0.10)

        vertex1 = (xlen*0.04,zlen*2)
        vertex2 = (xlen*0.04,zlen*0.45)
        vertex3 = (xlen*0.12,zlen*0.25)
        vertex4 = (xlen*0.12,zlen*0.20)
        
        vertex5 = (xlen*0.48,zlen*2)
        vertex6 = (xlen*0.48,zlen*0.45)
        vertex7 = (xlen*0.40,zlen*0.25)
        vertex8 = (xlen*0.40,zlen*0.20)

        vertex9 = (0,zlen*0.02)
        vertex10 = (xlen,zlen*0.02)

        # Draw top wall
        self.drawWallFromVtxs(vertex1,vertex5)
        # Draw bottom wall
        #self.drawWallFromVtxs(vertex4,vertex8)

        # Draw moving wall on bottom
        self.drawWallFromVtxs(vertex9,vertex10)

        # Draw left walls
        self.drawWallFromVtxs(vertex1,vertex2)
        self.drawWallFromVtxs(vertex2,vertex3)
        self.drawWallFromVtxs(vertex3,vertex4)
        # Draw right walls
        self.drawWallFromVtxs(vertex5,vertex6)
        self.drawWallFromVtxs(vertex6,vertex7)
        self.drawWallFromVtxs(vertex7,vertex8)
        
        
        radius = 4
        spacing = 4*radius
        x_range = math.floor((vertex1[0]+(vertex5[0]-vertex1[0]))/spacing)
        z_range = math.floor((vertex2[1]+(vertex1[1]-vertex2[1]))/spacing)
        for i in range(int(x_range)):
            for j in range(int(z_range)):
            	for k in range(int(z_range)):
	                x = (vertex1[0]+(vertex5[0]-vertex1[0])/x_range*(i+0.5))
	                y = 5
	                z = (vertex2[1]+(vertex1[1]-vertex2[1])/z_range*(j+0.5))
	                self.drawSphere(x,y,z,radius)

        #self.fillCubeWithCBDVtxs(vertex1,vertex6)


        #self.drawWall(xlen*0.20,zlen*0.95,xlen*0.20,zlen*0.40)
        #self.drawWall(xlen*0.80,zlen*0.95,xlen*0.80,zlen*0.40)
        

        self.drawSphere(xlen*0.45,ylen*0.5,zlen*0.45,ylen*0.45)

        self.drawSphere(80,100,100,10)
        self.drawSphere(120,100,100,10)
        # self.drawEllipsoid(100,100,120,40,40,10)
        #self.drawChickenNugget(70,5,120,40,1,40,"ACTIVE PARTICLE.png")
        # self.drawWall(150,150,150,50)


    def fillCubeWithCBDVtxs(self,v1,v2):
        (x1,z1) = v1
        (x2,z2) = v2
        self.fillCubeWithCBD(x1,self.y0,z1,x2,self.y1,z2)

    def fillCubeWithCBD(self,x,y,z,x2,y2,z2):
        xl = min(x,x2) + 2*self.dia
        xh = max(x,x2) - self.dia
        yl = min(y,y2) + self.dia
        yh = max(y,y2) - self.dia
        zl = min(z,z2) + self.dia
        zh = max(z,z2) - self.dia
        particleLines = copy.deepcopy(self.positionLines)
        activeParticleLines = []
        for line in particleLines:
            if line.split()[2] == str(self.active_type):
                activeParticleLines.append(line)
        print("Active particles:", len(activeParticleLines))
        minDistance = self.dia**2
        for xi in np.arange(xl,xh,self.dia):
            for yi in np.arange(yl,yh,self.dia):
                for zi in np.arange(zl,zh,self.dia):
                    placeParticle = True
                    for line in activeParticleLines:
                        data = line.split()
                        distance_sq = (xi - float(data[-3]))**2 + (yi - float(data[-2]))**2 + (zi - float(data[-1]))**2
                        if distance_sq < minDistance:
                            placeParticle = False
                            break
                    if (placeParticle): self.appendLine(self.cbd_type,xi,yi,zi)

    def appendLine(self,atomType,x,y,z):
        self.ID += 1
        self.positionLines.append(self.cbdStr % (self.ID, self.molID, atomType, x, y, z))

    def drawWallFromVtxs(self,vtx1,vtx2):
        (x,z)   = vtx1
        (x2,z2) = vtx2
        self.drawWall(x,z,x2,z2)

    def drawWall(self,x,z,x2,z2):
        self.molID += 1
        print("molID for wall is:", self.molID)
        length = math.sqrt((x2-x)**2 + (z2-z)**2)
        numPoints = math.ceil(length/(self.dia*5/8))
        delta_x = (x2-x)/numPoints
        delta_z = (z2-z)/numPoints
        for i in range(numPoints):
            xi = x+i*delta_x
            zi = z+i*delta_z
            for yi in np.arange(self.y0,self.y1,self.dia*5/8):
                self.appendLine(self.wall_type,xi,yi,zi)

    def drawSphere(self,x,y,z,radius):
            self.molID += 1
            for xi in np.arange(x-radius,x+radius,self.dia):
                for yi in np.arange(y-radius,y+radius,self.dia):
                    for zi in np.arange(z-radius,z+radius,self.dia):
                        if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) < radius**2:
                            self.appendLine(self.active_type,xi,yi,zi)

    '''def drawChickenNugget(self,x,y,z,rx,ry,rz,file):
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


   def drawSphere(self,x,y,z,radius):
        self.molID += 1
        for xi in np.arange(x-radius,x+radius,self.dia):
            for yi in np.arange(y-radius,y+radius,self.dia):
                for zi in np.arange(z-radius,z+radius,self.dia):
                    if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) < radius**2:
                        self.appendLine(self.active_type,xi,yi,zi)

    def drawEllipsoid(self,x,y,z,rx,ry,rz):
        self.molID += 1
        thk = self.dia*2
        for xi in np.arange(x-rx,x+rx,self.dia):
            for yi in np.arange(y-ry,y+ry,self.dia):
                for zi in np.arange(z-rz,z+rz,self.dia):
                    if ((xi-x)**2/(rx**2) + (yi-y)**2/(ry**2) + (zi-z)**2/(rz**2)) < 1.0:
                       if ((xi-x)**2/((rx-thk)**2) + (yi-y)**2/((ry-thk)**2) + (zi-z)**2/((rz-thk)**2)) > 1.0: 
                            self.appendLine(self.active_type,xi,yi,zi)

    #def drawDiParticle(self,x,y,z):
        self.molID += 1
        self.appendLine(self.active_type,x-self.dia/2,y,z)
        self.appendLine(self.active_type,x+self.dia/2,y,z)

    #def differentparticles(self,x,y,z):
        self.molID += 1
        self.appendLine(self.active_type,x-self.dia/2,y+1,z)
        self.appendLine(self.active_type,x+self.dia/2,y-1,z)'''

        


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