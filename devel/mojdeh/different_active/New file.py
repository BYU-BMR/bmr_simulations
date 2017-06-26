# Python code to extract equilibrium drop coordinates and add the dr blade
import numpy as np
import random

newFileName = "tutorial.data"


positionLines = []
# data string containing molecule ID, type, dia, rho, x, y, z, 0 0 0 
#cbdStr = "%d %d 0.93 0.0 1.0 %f %f %f\n" #string for cbd particles
cbdStr = "%d %d %f %f %f\n" #string for cbd particles

m_cbd = 3.896 # results in dia = 2 for cbd particles
rho_cbd = 0.93


vol = m_cbd/rho_cbd
dia = (3.0*vol/(4.0*3.1415))**(1.0/3.0)*2.0
dia = 2.0

solid_d = 3
x0 = 0.0
x1 = int(10.0*dia)
x2 = x1 + solid_d
y0 = 0.0
y1 = int(10.0*dia)
z0 = 0.0
z1 = 1.0
z2 = int(10.0*dia)
z3 = z2 + 1.0

ID = 0



# Add active particles to conglomerate
for x in np.arange(x1,x2,dia/30):
    for y in np.arange(y0,y0+solid_d,dia/30):
        for z in np.arange(z1,z1+solid_d,dia/30):
            ID += 1
            positionLines.append(cbdStr % (ID, 2, x, y, z))

    

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
    
