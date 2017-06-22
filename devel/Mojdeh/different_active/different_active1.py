# Python code to extract equilibrium drop coordinates and add the dr blade
import numpy as np
#from PIL import Image
import random, math, copy

class DatafileGenerator():
    newFileName = "active_spheres6.data"
    positionLines = []

    # data string containing molecule ID, type, dia, rho, x, y, z, 0 0 0 
    #cbdStr = "%d %d %d 0.93 0.0 1.0 %f %f %f\n" #string for cbd particles
    cbdStr = "%d %d %d %f %f %f\n" #string for cbd particles

    m_cbd = 3.72
    dia = 2.0
    act_dia = 6

    activeShapes = ['sphere','ellipsoid','rod','single','diparticle']
    activeShape = activeShapes[4]
    if activeShape == 'sphere':
        active_dia = 8.0
    elif activeShape == 'ellipsoid':
        active_dia = 18.0
    elif activeShape == 'rod':
        active_dia = 12.0
    elif activeShape == 'single':
        active_dia = 3.0
    elif activeShape == 'diparticle':
        active_dia = 4.0


    cbd_type,solvent_type,active_type,wall_type = 1,2,3,4

    scale = 0.50

    shrink = .50

    x0 = 0.0
    x1 = scale*3*50.0*dia
    y0 = 0.0
    y1 = scale*1*8.0*dia
    z0 = 0
    z1 = scale*4*50.0*dia*shrink

    ID = 0
    molID = 0

    def generateDataFile(self):
        self.setUpSimulation()
        self.writeFile()

    def setUpSimulation(self):
        xlen = abs(self.x1-self.x0)
        ylen = abs(self.y1-self.y0)
        zlen = abs((self.z1-self.z0)/self.shrink)

        # Set up walls
        # vertex1 = (xlen*0.20,zlen*0.95)
        # vertex2 = (xlen*0.20,zlen*0.40)
        # vertex3 = (xlen*0.40,zlen*0.20)
        # vertex4 = (xlen*0.40,zlen*0.10)
        
        # vertex5 = (xlen*(1.0-0.20),zlen*0.95)
        # vertex6 = (xlen*(1.0-0.20),zlen*0.40)
        # vertex7 = (xlen*(1.0-0.40),zlen*0.20)
        # vertex8 = (xlen*(1.0-0.40),zlen*0.10)

        opening_width = zlen*0.075
        print("opening_width:",opening_width)
        side_thickness = zlen*(7/4)*0.04
        xi = zlen*(7/4)*0.04
        xf = xi + side_thickness
        vertex1 = (xi,zlen*(self.shrink - .05))
        vertex2 = (xi,zlen*0.35/2)
        vertex3 = (xf,zlen*0.25/2)
        vertex4 = (xf,zlen*0.05/2+opening_width)
        
        xi = xf + opening_width
        xf = xi + side_thickness
        vertex5 = (xf,zlen*(self.shrink - .05))
        vertex6 = (xf,zlen*0.35/2)
        vertex7 = (xi,zlen*0.25/2)
        vertex8 = (xi,zlen*0.05/2+opening_width)

        vertex9 = (0,zlen*0.05/2)
        vertex10 = (xlen,zlen*0.05/2)

        
        # Draw moving wall on bottom
        self.drawWallFromVtxs(vertex9,vertex10)

        # Only include the active particles if the opening_width is wide enough
        ##self.fillCubeWithRandomMixVtxs(vertex1,vertex6)
        # if opening_width > self.active_dia:
        #     self.fillCubeWithMixVtxs(vertex1,vertex6)
        # else:
        #     self.fillCubeWithCBDVtxs(vertex1,vertex6,checkForOverlap=False)

        # Fill the left half of the box with Active and the right side with CBD
        # middle_x_vertex = ((vertex1[0]+vertex5[0])/2,vertex1[1])
        # self.fillCubeWithActiveVtxs(middle_x_vertex,vertex2)
        # self.fillCubeWithCBDVtxs(middle_x_vertex,vertex6,checkForOverlap=False)
        
        # # Fill top half of box with Active and the bottom half with CBD
        # middle_z_vertex = (vertex1[0],(vertex1[1]+vertex2[1])/2)
        self.fillCubeWithActiveVtxs(vertex1,vertex6)
        self.fillCubeWithCBDVtxs(vertex1,vertex6,checkForOverlap=True)

         # Draw top wall
        self.drawWallFromVtxs(vertex1,vertex5,atomType=5)
        # Draw bottom wall
        self.drawWallFromVtxs(vertex4,(vertex8[0]+self.dia,vertex8[1]),atomType=6)

        # Draw left walls
        self.drawWallFromVtxs(vertex1,vertex2)
        self.drawWallFromVtxs(vertex2,vertex3)
        self.drawWallFromVtxs(vertex3,vertex4)
        # Draw right walls
        self.drawWallFromVtxs(vertex5,vertex6)
        self.drawWallFromVtxs(vertex6,vertex7)
        self.drawWallFromVtxs(vertex7,vertex8)
        
        # self.fillCubeWithCBDVtxs(vertex1,vertex6)


        #self.drawWall(xlen*0.20,zlen*0.95,xlen*0.20,zlen*0.40)
        #self.drawWall(xlen*0.80,zlen*0.95,xlen*0.80,zlen*0.40)
        

        # self.drawSphere(xlen*0.45,ylen*0.5,zlen*0.45,ylen*0.45)

        # self.drawSphere(80,100,100,10)
        # self.drawSphere(120,100,100,10)
        # self.drawEllipsoid(100,100,120,40,40,10)
        # self.drawChickenNugget(50,50,50,40,1,40,"test.jpg")
        # self.drawWall(150,150,150,50)

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
        for yi in np.arange(yl,yh-self.dia,self.dia):
            for zi in np.arange(zl,zh-self.dia,self.dia):
                for xi in np.arange(xl,xh-self.dia*2,self.dia*2):
                    val = random.randint(0,99)
                    if val >= 0 and val < 18:
                        atom_type = self.active_type
                        self.molID += 1
                    elif val >= 18 and val < 32:
                        atom_type = self.cbd_type
                    else:
                        atom_type = self.solvent_type
                    self.appendLine(atom_type,xi+self.dia*1/2,yi+self.dia/2,zi+self.dia/2)
                    self.appendLine(atom_type,xi+self.dia*3/2,yi+self.dia/2,zi+self.dia/2)



    def fillCubeWithMixVtxs(self,v1,v2):
        (x1,z1) = v1
        (x2,z2) = v2
        self.fillCubeWithMix(x1,self.y0,z1,x2,self.y1,z2)

    def fillCubeWithMix(self,x,y,z,x2,y2,z2):
        xl = min(x,x2)
        xh = max(x,x2)
        yl = min(y,y2)
        yh = max(y,y2)
        zl = min(z,z2)
        zh = max(z,z2)
        radius = self.active_dia/2.0
        spacing = 2.0*radius
        x_act_range = math.floor((xh-xl)/spacing)
        y_act_range = math.floor((yh-yl)/(spacing/2))
        z_act_range = math.floor((zh-zl)/(spacing/2))
        x_cbd_range = math.floor((spacing)/self.dia)
        y_cbd_range = math.floor((yh-yl)/self.dia)
        z_cbd_range = math.floor((zh-zl)/self.dia)
        x_act_spacing = (xh-xl)/x_act_range
        y_act_spacing = (yh-yl)/y_act_range
        z_act_spacing = (zh-zl)/z_act_range
        x_cbd_spacing = (spacing)/x_cbd_range
        y_cbd_spacing = (yh-yl)/y_cbd_range
        z_cbd_spacing = (zh-zl)/z_cbd_range
        for i in range(int(x_act_range)):
            if (i % 2) == 1:
                for j in range(int(y_act_range)):
                    for k in range(int(z_act_range)):
                        x = xl+x_act_spacing*(i+0.5)
                        y = yl+y_act_spacing*(j+0.5)
                        z = zl+z_act_spacing*(int(z_act_range-1)-k+0.5)
                        if self.activeShape == 'sphere':
                            self.drawSphere(x,y,z,radius)
                        elif self.activeShape == 'ellipsoid':
                            self.drawEllipsoid(x,y,z,radius/2,radius,radius)  
                        elif self.activeShape == 'rod':
                            self.drawRod(x,y,z,radius/2,radius,radius/2)  
                        elif self.activeShape == 'single':
                            self.appendLine(self.active_type,x,y,z)
                        elif self.activeShape == 'diparticle':
                            self.drawDiParticle(x,y,z)    
            else:
                xi = xl+x_act_spacing*(i)
                for i_cbd in range(int(x_cbd_range)):
                    for j in range(int(y_cbd_range)):
                        for k in range(int(z_cbd_range)):
                            x = xi+x_cbd_spacing*(i_cbd+0.5)
                            y = yl+y_cbd_spacing*(j+0.5)
                            z = zl+z_cbd_spacing*(k+0.5)
                            self.appendLine(self.cbd_type,x,y,z)


    def fillCubeWithActiveVtxs(self,v1,v2):
        (x1,z1) = v1
        (x2,z2) = v2
        self.fillCubeWithActive(x1,self.y0,z1,x2,self.y1,z2)

    def fillCubeWithActive(self,x,y,z,x2,y2,z2):
        xl = min(x,x2)
        xh = max(x,x2)
        yl = min(y,y2)
        yh = max(y,y2)
        zl = min(z,z2)
        zh = max(z,z2)
        radius = self.act_dia/2
        spacing = 2*radius
        x_range = math.floor((xh-xl)/spacing)
        y_range = math.floor((yh-yl)/spacing)
        z_range = math.floor((zh-zl)/spacing)
        x_spacing = (xh-xl)/x_range
        y_spacing = (yh-yl)/y_range
        z_spacing = (zh-zl)/z_range 
        for i in range(int(x_range)):
            for j in range(int(y_range)):
                for k in range(int(z_range)):
                    x = xl+x_spacing*(i+0.5)
                    y = yl+y_spacing*(j+0.5)
                    z = zl+z_spacing*(int(z_range-1)-k+0.5)
                    self.drawSphere(x,y,z,radius)

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
        spacing = self.dia
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
                        if (placeParticle): self.appendLine(self.cbd_type,xi,yi,zi)
        else:
           for i in range(int(x_range)):
                for j in range(int(y_range)):
                    for k in range(int(z_range)):
                        xi = xl+x_spacing*(i+0.5)
                        yi = yl+y_spacing*(j+0.5)
                        zi = zl+z_spacing*(k+0.5)
                        self.appendLine(self.cbd_type,xi,yi,zi)

    def appendLine(self,atomType,x,y,z):
        self.ID += 1
        self.positionLines.append(self.cbdStr % (self.ID, self.molID, atomType, x, y, z))

    def drawWallFromVtxs(self,vtx1,vtx2,atomType=wall_type):
        (x,z)   = vtx1
        (x2,z2) = vtx2
        self.drawWall(x,z,x2,z2,atomType)

    def drawWall(self,x,z,x2,z2,atomType=wall_type):
        self.molID += 1
        print("molID for wall is:", self.molID)
        length = math.sqrt((x2-x)**2 + (z2-z)**2)
        numPoints = math.ceil(length/(self.dia*1))
        delta_x = (x2-x)/numPoints
        delta_z = (z2-z)/numPoints
        for i in range(numPoints):
            xi = x+i*delta_x
            zi = z+i*delta_z
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