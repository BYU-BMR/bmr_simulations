import numpy as np
import random

dataFileName = "visc.data"

positionLines = []
# data string containing molecule ID, type, dia, rho, x, y, z, 0 0 0 
cbdStr = "%d %d 0.93 0.0 1.0 %f %f %f\n" #string for cbd particles- this one will work for the final simulation
#cbdStr = "%d %d %f %f %f\n" #string for cbd particles-this one will work for obserivng the general shape by .data file

m_cbd = 3.896 # results in dia = 2 for cbd particles
rho_cbd = 0.93

m_active = m_cbd

vol = m_cbd/rho_cbd
dia = (3.0*vol/(4.0*3.1415))**(1.0/3.0)*2.0
dia_act = 2.0
dia_cbd = 2.4

solid_d = 3

xmin = 0.0
xmax = xmin + int(20.0*dia_cbd)
ymin = 0.0
ymax = ymin + int(10.0*dia_cbd)
zmin = 0.0
zmax = zmin + int(23.0*dia_cbd)

#z_wall_low = zmin + 3*dia_cbd
#z_wall_high = zmax - 3*dia_cbd

ID = 0
# Add CBD particles to grid
for x in np.arange(xmin,xmax,dia_cbd):
    for y in np.arange(ymin,ymax,dia_cbd):
        for z in np.arange(zmin,zmax,dia_cbd):
            xc = x + random.randrange(-100,100)/1000
            yc = y + random.randrange(-100,100)/1000
            zc = z + random.randrange(-100,100)/1000
            ID += 1
            positionLines.append(cbdStr % (ID, 1, xc, yc, zc))

# az inja minevise data file ro

print("Writing data file: %s" % dataFileName)

with open(dataFileName,"w") as outFile:

    outFile.write("\n")
    outFile.write("%d atoms\n" % ID)
    outFile.write("\n")
    outFile.write("%d atom types\n" % 3)
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
    

    outFile.write("\n")
#    outFile.write("Atoms # meso\n")
    outFile.write("Atoms\n")
    outFile.write("\n")

    for line in positionLines:
        outFile.write(line)