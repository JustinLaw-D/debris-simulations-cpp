# contains methods to update Constants.h, specifically the atmospheric data

import os
import numpy as np
from datetime import date

G = 6.67430e-11 # gravitational constant (N*m^2/kg^2)
Me = 5.97219e24 # mass of Earth (kg)
Re = 6371.0 # radius of Earth (km)
SCOL = ";"
ENDL = "\n"

filepath = os.path.dirname(os.path.abspath(__file__)) # path to current folder
filepath += '/'

# read density model
atmfile=open(filepath + "atmospheric_data/cira-2012.dat","r")
header=atmfile.readline()
zmodel=[]
denmodelL=[]
denmodelM=[]
denmodelHL=[]
for line in atmfile:
  alt,low,med,highL,_=line.split()
  zmodel.append(float(alt))
  denmodelL.append(float(low))
  denmodelM.append(float(med))
  denmodelHL.append(float(highL))

atmfile.close()

zmodel=np.array(zmodel)*1000 # convert to m
denmodelL=np.array(denmodelL)
denmodelM=np.array(denmodelM)
denmodelHL=np.array(denmodelHL)

logdenL = np.log10(denmodelL)
logdenM = np.log10(denmodelM)
logdenHL = np.log10(denmodelHL)
logz = np.log10(zmodel)

# read solar cycle template (using F10.7 as the solar activity index)
f107file = open(filepath + "atmospheric_data/solar_cycle_table36_cira2012.dat","r")
header=f107file.readline()
f107_mo=[]
for line in f107file:
   mo,_,f107,_,_,_,_,_,_,_,_,_,_=line.split()
   f107_mo.append(float(f107))
f107file.close()
f107_mo=np.array(f107_mo) 

f107file = open(filepath + "atmospheric_data/solar_cycle_table36_cira2012.dat","r")
header=f107file.readline()
f107_mo=[]
for line in f107file:
   mo,_,f107,_,_,_,_,_,_,_,_,_,_=line.split()
   f107_mo.append(float(f107))
f107file.close()
f107_mo=np.array(f107_mo)

# write the new Constants.h

consts_file = open(filepath + "Constants.h", 'w')
consts_file.write("// file to hold the constants used by the program" + ENDL)
consts_file.write("// last updated : " + str(date.today()) + ENDL)
consts_file.write(ENDL)
consts_file.write("#include <cstddef>" + ENDL)
consts_file.write("#pragma once" + ENDL)
consts_file.write(ENDL)
consts_file.write("// physics constants" + ENDL)
consts_file.write("const double G = " + str(G) + SCOL + ENDL)
consts_file.write("const double Re = " + str(Re) + SCOL + ENDL)
consts_file.write("const double Me = " + str(Me) + SCOL + ENDL)
consts_file.write(ENDL)
consts_file.write("// atmospheric data" + ENDL)
logz_str, logdenL_str, logdenM_str, logdenHL_str = str(logz[0]), str(logdenL[0]), str(logdenM[0]), str(logdenHL[0])
f107_mo_str = str(f107_mo[0])
num_alts = str(len(logz))
num_months = str(len(f107_mo))

for i in range(1, len(logz)):
    logz_str += "," + str(logz[i])
    logdenL_str += "," + str(logdenL[i])
    logdenM_str += "," + str(logdenM[i])
    logdenHL_str += "," + str(logdenHL[i])
for i in range(1, len(f107_mo)):
    f107_mo_str += "," + str(f107_mo[i])

consts_file.write("const size_t num_alts = " + num_alts + SCOL + ENDL)
consts_file.write("const double logz[" + num_alts + "] = {" + logz_str + "}" + SCOL + ENDL)
consts_file.write("const double logdenL[" + num_alts + "] = {" + logdenL_str + "}" + SCOL + ENDL)
consts_file.write("const double logdenM[" + num_alts + "] = {" + logdenM_str + "}" + SCOL + ENDL)
consts_file.write("const double logdenHL[" + num_alts + "] = {" + logdenHL_str + "}" + SCOL + ENDL)
consts_file.write("const size_t num_months = " + num_months + SCOL + ENDL)
consts_file.write("const double f107_mo[" + num_months + "] = {" + f107_mo_str + "}" + SCOL + ENDL)
consts_file.close()