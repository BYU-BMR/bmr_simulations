# Python code to extract equilibrium drop coordinates and add the dr blade
import numpy as np
import random

newFileName = "different_particles.data"


positionLines = []
# data string containing molecule ID, type, dia, rho, x, y, z, 0 0 0 
cbdStr = "%d %d 0.93 0.0 1.0 %f %f %f\n" #string for cbd particles
# cbdStr = "%d %d %f %f %f\n" #string for cbd particles

m_cbd = 3.896 # results in dia = 2 for cbd particles
rho_cbd = 0.93

m_active = m_cbd


vol = m_cbd/rho_cbd
# dia = (3.0*vol/(4.0*3.1415))**(1.0/3.0)*2.0
dia_cbd = 2.0
dia_act = 4.0

solid_d = 3
x0 = 0.0
x1 = x0 + int(10.0*dia_cbd)
x2 = x1 + int(10.0*dia_cbd)
y0 = 0.0
y1 = int(10.0*dia_cbd)
z0 = 0.0
z1 = int(10.0*dia_cbd)


ID = 0

# Add CBD particles to grid
for x in np.arange(x0+1,x1-1,dia_cbd):
    for y in np.arange(y0+1,y1-1,dia_cbd):
        for z in np.arange(z0+1,z1-1,dia_cbd):
            xc = x + random.randrange(-100,100)/1000
            yc = y + random.randrange(-100,100)/1000
            zc = z + random.randrange(-100,100)/1000
            ID += 1
            positionLines.append(cbdStr % (ID, 1, xc, yc, zc))


# Add Active particles to grid
for x in np.arange(x1+1,x2-1,dia_act):
    for y in np.arange(y0+1,y1-1,dia_act):
        for z in np.arange(z0+1,z1-1,dia_act):
            xc = x + random.randrange(-100,100)/1000
            yc = y + random.randrange(-100,100)/1000
            zc = z + random.randrange(-100,100)/1000
            ID += 1
            positionLines.append(cbdStr % (ID, 2, xc, yc, zc))



print("Writing new file: %s" % newFileName)

xmin = x0
xmax = x2
ymin = y0
ymax = y1
zmin = z0
zmax = z1

with open(newFileName,"w") as outFile:

    outFile.write("\n")
    outFile.write("%d atoms\n" % ID)
    outFile.write("\n")
    outFile.write("%d atom types\n" % 2)
    outFile.write("\n")
    outFile.write("%f %f xlo xhi\n" % (xmin,xmax))   
    outFile.write("%f %f ylo yhi\n" % (ymin,ymax))
    outFile.write("%f %f zlo zhi\n" % (zmin,zmax))

    outFile.write("\n")
    outFile.write("Masses\n")
    outFile.write("\n")
    outFile.write("1 %f\n" % m_cbd)
    outFile.write("2 %f\n" % m_active)
    

    outFile.write("\n")
#    outFile.write("Atoms # meso\n")
    outFile.write("Atoms\n")
    outFile.write("\n")

    for line in positionLines:
        outFile.write(line)
    
