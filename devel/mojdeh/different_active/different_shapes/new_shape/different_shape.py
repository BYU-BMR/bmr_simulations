# Python code to extract equilibrium drop coordinates and add the dr blade
import numpy as np
import random

newFileName = "active_particle.data"


positionLines = []
# data string containing molecule ID, type, dia, rho, x, y, z, 0 0 0 
#cbdStr = "%d %d 0.93 0.0 1.0 %f %f %f\n" #string for cbd particles
cbdStr = "%d %d %f %f %f\n" #string for cbd particles

m_cbd = 3.896 # results in dia = 2 for cbd particles
rho_cbd = 0.93
print (cbdStr)

vol = m_cbd/rho_cbd
dia = (3.0*vol/(4.0*3.1415))**(1.0/3.0)*2.0
dia = 2.0

solid_d = 3
x0 = 1.0
x1 = int(10.0*dia)
x2 = x1 + solid_d
y0 = 0.0
y1 = int(10.0*dia)
z0 = 0.0
z1 = 1.0
z2 = int(10.0*dia)
z3 = z2 + 1.0

ID = 0

'''# Add CBD particles to grid
for x in np.arange(x0+1,x1-1,dia):
    for y in np.arange(y0+1,y1-1,dia):
        for z in np.arange(z1+1,z2-1,dia):
            xc = x + random.randrange(-100,100)/1000
            yc = y + random.randrange(-100,100)/1000
            zc = z + random.randrange(-100,100)/1000
            ID += 1
            positionLines.append(cbdStr % (ID, 1, xc, yc, zc))'''

# Add active particles to conglomerate
'''for x in np.arange(x1,x2,dia/2.0):
    for y in np.arange(y0,y0+2*solid_d,dia/2.0):
        for z in np.arange(z1,z1+3*solid_d,dia/2.0):
            ID += 1
            positionLines.append(cbdStr % (ID, 2, x, y, z))'''

for x in np.arange(x0+1,x1-5.5*solid_d,dia/10):
    for y in np.arange(y1-solid_d,y1+solid_d,dia/10):
        for z in np.arange(z0,z1,dia):
            ID += 1
            positionLines.append(cbdStr % (ID, 2, x, y, z))

'''for x in np.arange(x0,x1,dia/2.0):
    for y in np.arange(y1-solid_d,y1,dia/2.0):
        for z in np.arange(z3-solid_d,z3,dia/2.0):
            ID += 1
            positionLines.append(cbdStr % (ID, 4, x, y, z))

for x in np.arange(x1,x2,dia/2.0):
    for y in np.arange(y0,y0+solid_d,dia/2.0):
        for z in np.arange(z3-solid_d,z3,dia/2.0):
            ID += 1
            positionLines.append(cbdStr % (ID, 5, x, y, z))'''
            
    
'''# Add top and bottom plates for shearing
for x in np.arange(x0,x2,dia):
    for y in np.arange(y0,y1,dia):
        for z in np.arange(z0,z1,dia/4.0):
            ID += 1
            positionLines.append(cbdStr % (ID, 6, x, y, z))

for x in np.arange(x0,x2,dia):
    for y in np.arange(y0,y1,dia):
        for z in np.arange(z2,z3,dia/4.0):
            ID += 1
            positionLines.append(cbdStr % (ID, 7, x, y, z))'''

print("Writing new file: %s" % newFileName)

xmin = x0
xmax = x2
ymin = y0
ymax = y1
zmin = z0
zmax = z3

with open(newFileName,"w") as outFile:

    outFile.write("\n")
    outFile.write("%d atoms\n" % ID)
    outFile.write("\n")
    outFile.write("%d atom types\n" % 7)
    outFile.write("\n")
    outFile.write("%f %f xlo xhi\n" % (xmin,xmax))   
    outFile.write("%f %f ylo yhi\n" % (ymin,ymax))
    outFile.write("%f %f zlo zhi\n" % (zmin,zmax))

    outFile.write("\n")
    outFile.write("Masses\n")
    outFile.write("\n")
    outFile.write("1 %f\n" % m_cbd)
    outFile.write("2 %f\n" % m_cbd)
    outFile.write("3 %f\n" % m_cbd)
    outFile.write("4 %f\n" % m_cbd)
    outFile.write("5 %f\n" % m_cbd)
    outFile.write("6 %f\n" % m_cbd)
    outFile.write("7 %f\n" % m_cbd)

    outFile.write("\n")
#    outFile.write("Atoms # meso\n")
    outFile.write("Atoms\n")
    outFile.write("\n")

    for line in positionLines:
        outFile.write(line)
    
