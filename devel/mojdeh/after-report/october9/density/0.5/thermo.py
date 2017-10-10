# Python code to extract equilibrium drop coordinates and add the dr blade
import numpy as np
import random

newFileName = "random_configuration.data"


positionLines = []
# data string containing molecule ID, type, dia, rho, x, y, z, 0 0 0 
cbdStr = "%d %d %d 0.93 0.0 1.0 %f %f %f\n" #string for cbd particles
actStr = "%d %d %d 4.79 0.0 1.0 %f %f %f\n" #string for active particles
solStr = "%d %d %d 0.5 0.0 1.0 %f %f %f\n" #string for solvent particles
wallStr = "%d %d %d 5.0 0.0 1.0 %f %f %f\n" #string for wall particles
#cbdStr = "%d %d %d %f %f %f\n" #string for cbd particles
#solStr = "%d %d %d %f %f %f\n" #string for solvent particles
#actStr = "%d %d %d %f %f %f\n" #string for active particles
#wallStr = "%d %d %d %f %f %f\n" #string for wall particles

m_cbd = 3.896
#m_act = 20.064 #me before:     1284.11 
m_sol = 4.31      
#m_wall = 28.378    #me before: 2.61

dia = 2.0
#cbd_dia = 2.0
#act_dia = 6.0 #????is it correct?
sol_dia = 2.0
#wall_dia = 3.0

#cbd_type,active_type,solvent_type,wall_type = 1,2,3,4
solvent_type = 1
rho_sol = 1.028


vol = m_sol/rho_sol
# dia = (3.0*vol/(4.0*3.1415))**(1.0/3.0)*2.0
dia_cbd = 2.0
#dia_act = 2.0

solid_d = 3
x0 = 0.0
x1 = x0 + int(10.0*dia_cbd)
x2 = x1 
y0 = 0.0
y1 = int(10.0*dia_cbd)
z0 = 0.0
z1 = int(10.0*dia_cbd)


ID = 0

# Add Solvent particles to grid
for x in np.arange(x0+2.5*dia_cbd,x1-1,1.6*sol_dia):
    for y in np.arange(y0+1,y1-1,1.6*sol_dia):
        for z in np.arange(z0+1,z1-1,1.6*sol_dia):
            xc = x + random.randrange(-100,100)/100
            yc = y + random.randrange(-100,100)/100
            zc = z + random.randrange(-100,100)/100
            ID += 1
            positionLines.append(solStr % (ID, solvent_type, sol_dia, xc, yc, zc))


'''# Add Active particles to grid
for x in np.arange(x1+1,x2-1,dia_act):
    for y in np.arange(y0+1,y1-1,dia_act):
        for z in np.arange(z0+1,z1-1,dia_act):
            xc = x + random.randrange(-100,100)/1000
            yc = y + random.randrange(-100,100)/1000
            zc = z + random.randrange(-100,100)/1000
            ID += 1
            positionLines.append(cbdStr % (ID, 2, xc, yc, zc))'''



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
    #outFile.write("2 %f\n" % m_act)
    outFile.write("2 %f\n" % m_sol)
    #outFile.write("4 %f\n" % m_wall)
    #outFile.write("5 %f\n" % m_wall)
    #outFile.write("6 %f\n" % m_wall)
    #outFile.write("7 %f\n" % m_wall)
    #outFile.write("8 %f\n" % m_wall)
    

    outFile.write("\n")
#    outFile.write("Atoms # meso\n")
    outFile.write("Atoms\n")
    outFile.write("\n")

    for line in positionLines:
        outFile.write(line)
    
