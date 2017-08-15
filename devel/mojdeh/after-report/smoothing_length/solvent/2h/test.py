# Python code to extract equilibrium drop coordinates and add the dr blade
import numpy as np
#from PIL import Image
import random, math, copy

class DatafileGenerator():
    newFileName = "slvnt_clmn.data"
    positionLines = []

   # data string containing atom ID, molecule ID, type, rho, e, cv, x, y, z 
    cbdStr = "%d %d %d 0.93 0.0 1.0 %f %f %f\n" #string for cbd particles
    solStr = "%d %d %d 1.028 0.0 1.0 %f %f %f\n" #string for solvent particles
    actStr = "%d %d %d 4.79 0.0 1.0 %f %f %f\n" #string for active particles
    wallStr = "%d %d %d 2.0 0.0 1.0 %f %f %f\n" #string for wall particles
    #cbdStr = "%d %d %d %f %f %f\n" #string for cbd particles
    #solStr = "%d %d %d %f %f %f\n" #string for solvent particles
    #actStr = "%d %d %d %f %f %f\n" #string for active particles
    #wallStr = "%d %d %d %f %f %f\n" #string for wall particles

    m_cbd = 3.896
    m_sol = 4.306
    m_act = 20.064
    m_wall = 8.378

    dia = 2.0
    cbd_dia = 2.0
    sol_dia = 2.0
    act_dia = 8.0
    wall_dia = 2.0

    cbd_type,solvent_type,active_type,wall_type = 1,2,3,4

    active_thickness_factor = 1

    True_activecount = 0
    activecount = 0
    solventcount = 0
    cbdcount = 0

    '''activeShapes = ['sphere','ellipsoid','rod','single','diparticle']
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
        active_dia = 4.0'''

    scale = 1

    shrink = 1

    x0 = 0.0
    x1 = scale*1*20.0*dia
    y0 = 0.0
    y1 = scale*1*20.0*dia
    z0 = 0
    z1 = scale*1*50.0*dia*shrink

    ID = 0
    molID = 1

    def generateDataFile(self):
        self.setUpSimulation()
        self.writeFile()

    def setUpSimulation(self):
        xlen = abs(self.x1-self.x0)
        ylen = abs(self.y1-self.y0)
        zlen = abs((self.z1-self.z0)/self.shrink)

        vertexalpha = (.4*xlen,zlen*.05 + self.dia)
        vertexbeta =  (.6*xlen,zlen*.9)

        # Set up walls
        # vertex1 = (xlen*0.20,zlen*0.95)
        # vertex2 = (xlen*0.20,zlen*0.40)
        # vertex3 = (xlen*0.40,zlen*0.20)
        # vertex4 = (xlen*0.40,zlen*0.10)
        
        # vertex5 = (xlen*(1.0-0.20),zlen*0.95)
        # vertex6 = (xlen*(1.0-0.20),zlen*0.40)
        # vertex7 = (xlen*(1.0-0.40),zlen*0.20)
        # vertex8 = (xlen*(1.0-0.40),zlen*0.10)

        '''opening_width = 2*self.act_dia
        print("opening_width:",opening_width)
        side_thickness = zlen*(7/4)*0.04
        xi = zlen*(7/4)*0.04
        xf = xi + side_thickness

        vertex1 = (xi,zlen*(self.shrink - .05)+self.act_dia)
        vertex2 = (xi,zlen*0.35/2)
        vertex3 = (xf,zlen*0.25/2)
        vertex4 = (xf,zlen*0.05/2+opening_width)
    
        ai = xi
        xi = xf + opening_width
        xf = xi + side_thickness

        vertex5 = (xf,zlen*(self.shrink - .05)+self.act_dia)
        vertex6 = (xf,zlen*0.35/2)
        vertex7 = (xi,zlen*0.25/2)
        vertex8 = (xi,zlen*0.05/2+opening_width)
       
        bi = xf'''
        vertex9 = (0,zlen*0.05/2)
        vertex10 = (xlen,zlen*0.05/2)

        vertex1 = (0,zlen*0.98)
        vertex2 = (xlen,zlen*0.98)

        '''ci = math.floor(bi-ai) + 1
        di = opening_width

        vertex11 = (vertex1[0]+ci,vertex1[1])
        vertex41 = (vertex4[0]+di,vertex4[1])
        vertex51 = (vertex5[0]+ci,vertex5[1])
        vertex81 = (vertex8[0]+di,vertex8[1])

        vertex20 = (vertex2[0],0)
        vertex24 = (vertex2[0],vertex4[1]-self.wall_dia)
        vertex60 = (vertex6[0],0)
        vertex68 = (vertex6[0],vertex8[1]-self.wall_dia)'''

        
        
        # Draw moving wall on bottom
        self.drawWallFromVtxs(vertex9,vertex10)
        #self.drawWallFromVtxs2(vertex1,vertex2,atomType=5)
        #self.fillCubeWithRandomMixVtxs(vertex1,vertex6)
        self.fillCubeWithCBDVtxs2(vertexalpha,vertexbeta)
        #self.placeparticle2(.875*xlen,.5*ylen,.5*zlen)
        # Only include the active particles if the opening_width is wide enough
        #self.fillCubeWithActiveVtxs(vertex1,vertex6)
        #self.fillCubeWithCBDVtxs(vertex1,vertex6,checkForOverlap=True)

         # Draw top wall
        #self.drawWallFromVtxs(vertex11,vertex51,atomType=5)
        # Draw bottom wall
        #self.drawWallFromVtxs(vertex41,vertex81,atomType=6)

        # Draw Guards
        #self.drawWallFromVtxs(vertex20,vertex2,atomType=7)
        #self.drawWallFromVtxs(vertex60,vertex6,atomType=7)

        # Draw left walls
        #self.drawWallFromVtxs(vertex1,vertex2)
        #self.drawWallFromVtxs(vertex2,vertex3)
        #self.drawWallFromVtxs(vertex3,vertex4)
        # Draw right walls
        #self.drawWallFromVtxs(vertex5,vertex6)
        #self.drawWallFromVtxs(vertex6,vertex7)
        #self.drawWallFromVtxs(vertex7,vertex8)
        
        # self.fillCubeWithCBDVtxs(vertex1,vertex6)


        #self.drawWall(xlen*0.20,zlen*0.95,xlen*0.20,zlen*0.40)
        #self.drawWall(xlen*0.80,zlen*0.95,xlen*0.80,zlen*0.40)
        

        # self.drawSphere(xlen*0.45,ylen*0.5,zlen*0.45,ylen*0.45)

        # self.drawSphere(80,100,100,10)
        # self.drawSphere(120,100,100,10)
        # self.drawEllipsoid(100,100,120,40,40,10)
        # self.drawChickenNugget(50,50,50,40,1,40,"test.jpg")
        # self.drawWall(150,150,150,50)

    def placeparticle2(self,x,y,z):
        self.drawRaspberry(x,y,z,self.act_dia/2)

    def fillCubeWithRandomMixVtxs(self,v1,v2):
        (x1,z1) = v1
        (x2,z2) = v2

        self.fillCubeWithRandomMix(x1,self.y0,z1,x2,self.y1,z2)
       

    def fillCubeWithRandomMix(self,x,y,z,x2,y2,z2):
        cbd_to_active_ratio = math.floor(self.act_dia/self.cbd_dia)
        xl = min(x,x2)
        xh = max(x,x2)
        yl = min(y,y2)
        yh = max(y,y2)
        zl = min(z,z2)
        zh = max(z,z2)
        for yi in np.arange(yl+self.act_dia/2+self.cbd_dia/2,yh-self.act_dia/2-self.cbd_dia/2,self.act_dia):
            for zi in np.arange(zl+self.act_dia/2+self.cbd_dia/2,zh-self.act_dia/2-self.cbd_dia/2,self.act_dia):
                for xi in np.arange(xl+self.act_dia/2+self.cbd_dia/2,xh-self.act_dia/2-self.cbd_dia/2,self.act_dia):
                    val = random.randint(0,99)
                    if val >= 0 and val < 21:
                        atom_type = self.active_type
                        self.molID += 1
                        #self.appendLine(atom_type,xi,yi,zi)
                        self.drawRaspberry(xi,yi,zi,self.act_dia/2)
                    else:
                        for y in np.arange(yi-self.act_dia/2 + self.cbd_dia/2,yi + self.act_dia/2,self.cbd_dia):
                            for z in np.arange(zi-self.act_dia/2 + self.cbd_dia/2,zi + self.act_dia/2,self.cbd_dia):
                                for x in np.arange(xi-self.act_dia/2 + self.cbd_dia/2,xi + self.act_dia/2,self.cbd_dia):
                                    value = random.randint(0,79)
                                    if value >=0 and value < 66:
                                        atom_type = self.solvent_type
                                        self.appendLine(atom_type,x,y,z)
                                    else:
                                        atom_type = self.cbd_type
                                        self.appendLine(atom_type,x,y,z)

    def fillCubeWithCBDVtxs2(self,v1,v2):
        (x1,z1) = v1
        (x2,z2) = v2
        (y1,y2) = (0.25*self.y1,0.75*self.y1)
        self.fillCubeWithCBD2(x1,y1,z1,x2,y2,z2)
       

    def fillCubeWithCBD2(self,x,y,z,x2,y2,z2):
        xl = min(x,x2)
        xh = max(x,x2)
        yl = min(y,y2)
        yh = max(y,y2)
        zl = min(z,z2)
        zh = max(z,z2)
        for yi in np.arange(yl,yh,self.cbd_dia):
            for zi in np.arange(zl,zh,self.cbd_dia):
                for xi in np.arange(xl,xh,self.cbd_dia):
                    self.molID += 1
                    atom_type = self.solvent_type
                    self.appendLine(atom_type,xi,yi,zi)   


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
                for k in range(int(z_range/2)):
                    x = xl+x_spacing*(i+0.5)
                    y = yl+y_spacing*(j+0.5)
                    z = zl+z_spacing*(int(z_range-1)-k+0.5)
                    chance = random.randint(0,2)
                    if chance == 0:
                        self.drawdisc1(x,y,z,radius)
                    elif chance == 1:
                        self.drawdisc2(x,y,z,radius)
                    else:
                        self.drawdisc3(x,y,z,radius)

    '''def fillCubeWithCBDVtxs(self,v1,v2,checkForOverlap=True):
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
            print("Active particles:", len(activeParticleLines)*(3/4))
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
                            self.appendLine(self.solvent_type,xi,yi,zi)'''

    def appendLine(self,atomType,x,y,z):
        self.ID += 1
        if atomType == self.cbd_type:
            self.positionLines.append(self.cbdStr % (self.ID, self.molID, atomType, x, y, z))
        elif atomType == self.solvent_type:
            self.positionLines.append(self.solStr % (self.ID, self.molID, atomType, x, y, z))
        elif atomType == self.active_type:
            self.positionLines.append(self.actStr % (self.ID, self.molID, atomType, x, y, z))
        else:
            self.positionLines.append(self.wallStr % (self.ID, self.molID, atomType, x, y, z))

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

    def drawWallFromVtxs2(self,vtx1,vtx2,atomType=wall_type):
        (x,z)   = vtx1
        (x2,z2) = vtx2
        self.drawWall2(x,z,x2,z2,atomType)

    def drawWall2(self,x,z,x2,z2,atomType=wall_type):
        self.molID += 1
        print("molID for wall is:", self.molID)
        length = math.sqrt((x2-x)**2 + (z2-z)**2)
        numPoints = math.ceil(length/(self.dia*.5))
        delta_x = (x2-x)/numPoints
        delta_z = (z2-z)/numPoints
        for i in range(numPoints):
            xi = x+i*delta_x
            zi = z+i*delta_z
            for yi in np.arange(self.y0,self.y1,self.dia*1):
                self.appendLine(atomType,xi,yi,zi)

    def drawdisc1(self,atomType,x,y,z):
        r_act = self.act_dia/2
        atf = self.active_thickness_factor
        pi2 = 2*np.pi
        for ri in np.arange(0,r_act,atf*1.5):
            for thetai in np.arange(0,pi2,atf/4):
                X = x + ri*np.cos(thetai)
                Y = y + ri*np.sin(thetai)
                Z = z
                self.True_activecount += 1
                self.appendLine(atomType,X,Y,Z)

    def drawdisc2(self,atomType,x,y,z):
        r_act = self.act_dia/2
        atf = self.active_thickness_factor
        pi2 = 2*np.pi
        for ri in np.arange(0,r_act,atf*1.5):
            for thetai in np.arange(0,pi2,atf/4):
                X = x + ri*np.cos(thetai)
                Y = y 
                Z = z + ri*np.sin(thetai)
                self.True_activecount += 1
                self.appendLine(atomType,X,Y,Z)

    def drawdisc3(self,atomType,x,y,z):
        r_act = self.act_dia/2
        atf = self.active_thickness_factor
        pi2 = 2*np.pi
        for ri in np.arange(0,r_act,atf*1.5):
            for thetai in np.arange(0,pi2,atf/4):
                X = x 
                Y = y + ri*np.cos(thetai)
                Z = z + ri*np.sin(thetai)
                self.True_activecount += 1
                self.appendLine(atomType,X,Y,Z)

    def drawRaspberry(self,x,y,z,radius):
        self.molID += 1
        for xi in np.arange(x-radius+self.dia/4,x+radius,self.dia*(3/4)):
            for yi in np.arange(y-radius+self.dia/4,y+radius,self.dia):
                for zi in np.arange(z-radius+self.dia/4,z+radius,self.dia):
                    if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) < radius**2:
                        #if ((xi-x)**2 + (yi-y)**2 + (zi-z)**2) > radius**1:  # making hollow
                            self.appendLine(self.active_type,xi,yi,zi)


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
            outFile.write("2 %f\n" % self.m_sol)
            outFile.write("3 %f\n" % self.m_act)
            outFile.write("4 %f\n" % self.m_wall)
            outFile.write("5 %f\n" % self.m_wall)
            outFile.write("6 %f\n" % self.m_wall)
            outFile.write("7 %f\n" % self.m_wall)
            outFile.write("8 %f\n" % self.m_wall)

            outFile.write("\n")
            outFile.write("Atoms\n")
            outFile.write("\n")

            for line in self.positionLines:
                outFile.write(line)
    
if __name__ == "__main__":
    generator = DatafileGenerator()
    generator.generateDataFile()