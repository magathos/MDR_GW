#import matplotlib
#matplotlib.use('Agg')
from pylab import *
from LIV_postproc import *
za = 0.01 + arange(30)*0.01
alpha = 0.0

Da=[]
Da0=[]
DL=[]
print shape(za)

print LCDM_Cosmology().values()
for z in za:
    print "z = " + str(z)
    ad = alphaDistance(z, LCDM_Cosmology(), alpha)
    da0 = Dalpha0_LO(z, LCDM_Cosmology())
    ld = LuminosityDistance(z, LCDM_Cosmology())
    print "D_alpha = " + str(ad) + "   DL = " + str(ld)
    Da.append(ad)
    Da0.append(da0)
    DL.append(ld)

fig = figure()
ax = fig.add_subplot(111)
ax.set_title("$\\alpha = " + str(alpha) + "$")
ax.set_xlim((0.0, 0.3))
ax.set_ylabel('cosmological distance [m]')
ax.set_xlabel('redshift')
ax.scatter(za,array(Da), color='g', label="$D_{\\alpha=" + str(alpha) + "}$")
ax.scatter(za,array(Da0), color='r', label="$D_{\\alpha=" + str(alpha) + "}^{LO}$")
ax.scatter(za,array(DL), color='b', label='$D_L$')
ax.legend(loc='upper left', fancybox=True)

fig.savefig("/home/spxma1/public_html/LVC/D-z_alpha0_small.png")
