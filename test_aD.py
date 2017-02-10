from pylab import *
from LIV_postproc import *
za = 0.1 + arange(20)*0.1
alpha = 2.0

print LCDM_Cosmology()
for z in za:
    print "z = " + str(z)
    print "D_alpha = " + str(alphaDistance(z, LCDM_Cosmology(), alpha)) + "   DL = " + str(LuminosityDistance(z, LCDM_Cosmology()))
