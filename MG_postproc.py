#!/bin/python

'''
Copyright (C) 2015 Michalis Agathos

Massive Dispersion Relation post-processing module for PE 

'''

import sys
import os
import argparse

from pylab import *
from scipy import stats
from scipy import integrate
from scipy.optimize import newton 
import random 


###########################################
#
#  Define hardcoded
#
###########################################

fig_width_pt = 3*246.0/1.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean # *2.0/3.0      # height in inches
fig_size =  [fig_width,fig_height]

plot_params = {'axes.grid' : True,
               'savefig.dpi' : 150,
               'axes.labelsize': 20, 
               'axes.titlesize': 20, 
               'font.size': 16, 
               'font.family': 'serif', 
               'legend.fontsize': 16, 
               'xtick.labelsize': 16, 
               'xtick.major.size': 8,
               'xtick.minor.size': 4,
               'ytick.labelsize': 16, 
               'axes.grid' : True,
               'text.usetex': False,
               'lines.markersize' : 4, 
               'figure.figsize': fig_size
               }
          
rcParams.update(plot_params)

lambdag_bounds = {"Solar system":(2.8e15,'r','-.'), "Binary Pulsars":(1.6e13,'g','--')}
C_SI = 2.99792458e8
hbar_eVs = 6.582e-16
h_eVs = hbar_eVs*2*pi
G_SI = 6.67384e-11
MT_SUN = 4.92549e-6
pc_SI = 3.0856775807e16

#########################################################
#
#  Define functions
#
#  Massive Graviton calculations are defined in
#  Will (1997) [arXiv:gr-qc/9709011]
#
#  Most cosmological distances can be
#  found in Hogg (2000) [arXiv:astro-ph/9905116] 
#  consistent with LAL definitions.
#
#########################################################

def getParams(paramfile):
  '''Populate recovery and injection parameters from files'''
  recparams = {}
  if os.path.isfile(paramfile):
    f = open(paramfile, 'r')
    l = f.readline().strip()
    f.close()
    parlist = l.split()
    recparams = dict(zip(parlist, arange(len(parlist))))
  return recparams

def GComptonWavelength(l_g, redshift):
    """
    Calculate Compton wavelength of the graviton in m as:
    l_g * sqrt((1 + (2+z)*(1+z+sqrt(1+z)))/(5*(1+z)^3))
    Valid for \Omega_0 = 1 and for all z.
    """
    return l_g*sqrt((1.0 + (2.0+redshift)*(1.0+redshift+sqrt(1.0+redshift)))/(5.0*(1.0+redshift)**3.0))

def GravitonMass(lambda_g):
    """
    Calculate graviton mass in eV as:    
    m c^2 = h c / \lambda_g
    Valid for \Omega_0 = 1 and for all z.
    """
#    return 1.23982e-6/(lambda_g*sqrt((1 + (2+redshift)*(1+redshift+sqrt(1+redshift)))/(5*(1+redshift)**3)))
    return 1.23982e-6/lambda_g

def LCDM_Cosmology(H0=70.0, OmegaM=0.3, OmegaL=0.7):
  
  cosmo_dict = {"H0":H0, "OmegaM":OmegaM, "OmegaL":OmegaL, "OmegaK":1.0-OmegaM-OmegaL}
  return cosmo_dict


def HubbleTime(H0):
  """
  Calculates the hubble time in seconds
  """
  return 1000.0*pc_SI/H0

def HubbleDistance(H0):
  """
  Calculates the Hubble distance in meters
  """
  return C_SI*1000.0*pc_SI/H0

def Lambda(cosmology):
  """
  Calculates the cosmological constant Lambda in m^-2
  """
  return 3.0*cosmology["OmegaL"]*(cosmology["H0"]/(1000.0*pc_SI))**2/C_SI**2

def LuminosityDistance(redshift, cosmology):
  """
  Calculate luminosity distance from redshift 
  and LambdaCDM parameters, in meters.
  """
  
  DL = (1.0 + redshift) * TransverseComovingDistance(redshift, cosmology)
  return DL

def MGDistance(redshift, cosmology):
  """
  Calculate the cosmological distance that
  is relevant for GW propagation with a 
  massive graviton dispersion relation.
  (Will 1997), in meters
  """
  print "ERROR: Not implemented yet!"
  exit(-1)
  return 0

def LOSComovingDistance(redshift, cosmology):
  """
  Calculates line-of-sight comoving distance in meters
  """
  DH = HubbleDistance(cosmology["H0"])
  OM = cosmology["OmegaM"]
  OK = cosmology["OmegaK"]
  OL = cosmology["OmegaL"]

  Hparam = lambda z:1.0/sqrt(OM*(1+z)**3 + OK*(1+z)**2 + OL)
  cosmoint = integrate.quad(Hparam, 0, redshift)
  DC = DH * cosmoint[0]
  return DC

def TransverseComovingDistance(redshift, cosmology):
  """
  Calculates the transverse comoving distance in meters
  """
  eps = 1e-7
  k = cosmology["OmegaK"]
  DH = HubbleDistance(cosmology["H0"])
  DC = LOSComovingDistance(redshift, cosmology)
  if (fabs(k)<eps):
    DM = DC
  elif (k>0.0):
    DM = DH*sinh(sqrt(k)*DC/DH)/sqrt(k)
  elif (k<0.0):
    DM = DH*sin(sqrt(-k)*DC/DH)/sqrt(-k)
  return DM

def AngularDiameterDistance(redshift, cosmology):
  """
  The angular diameter distance in meters
  """
  DM = TransverseComovingDistance(redshift, cosmology)
  return DM/(1.0+redshift)

def AngularDiameterDistance12(redshift1, redshift2, cosmology):
  """
  The angular diameter distance between two points
  at given redshifts z1 and z2, in meters.
  """

  DM1 = TransverseComovingDistance(redshift1, cosmology)
  DM2 = TransverseComovingDistance(redshift2, cosmology)
  DH = HubbleDistance(cosmology["H0"])
  DA12 = (DM2*sqrt(1.0 + cosmology["OmegaK"]*DM1**2/DH**2) - DM1*sqrt(1.0 + cosmology["OmegaK"]*DM2**2/DH**2))/(1.0 + redshift2)
  return DA12


def calculate_redshift(distance, cosmology=None):
    """
    Calculate the redshift from the luminosity distance measurement 
    given in meters. Returns an array of redshifts.
    """
    def find_z_root(redshift, dl, cosmo):
        return dl - LuminosityDistance(redshift, cosmo)
    
    if cosmology is None:
      cosmology = LCDM_Cosmology()
    z = array([newton(find_z_root, random.uniform(0.0,2.0), args = (d, cosmology)) for d in distance[:]])
    return z

def plotPosteriorHist(data, nbins=50, xmin=None, xmax=None, logplot=False, ax=None, xlabel=None, label="GW"):
  '''Plot posterior histogram'''
  if ax is None:
    fig = figure()
    ax = fig.add_subplot(111)
    if xlabel is not None:
      ax.set_xlabel(xlabel)
    ax.set_ylabel("Probability Density")
    if xmin is not None:
      ax.set_xlim(xmin=xmin)
    else:
      xmin = min(data)
    if xmax is not None:
      ax.set_xlim(xmax=xmax)
    else:
      xmax = max(data)
  
  sdata = sort(data)
  CI_50 = (sdata[floor(0.5*len(sdata))],sdata[ceil(0.5*len(sdata))])
  CI_68 = (sdata[floor(0.32*len(sdata))])
  CI_90 = (sdata[floor(0.1*len(sdata))])
  CI_95 = (sdata[floor(0.05*len(sdata))])
  CI_99 = (sdata[floor(0.01*len(sdata))])
#  CI_68 = (sdata[floor(0.16*len(sdata))], sdata[ceil(0.84*len(sdata))])
#  CI_95 = (sdata[floor(0.025*len(sdata))], sdata[ceil(0.975*len(sdata))])
#  CI_99 = (sdata[floor(0.005*len(sdata))], sdata[ceil(0.995*len(sdata))])
  # lowerbound = min(data)
  print "shape:", shape(data), " min:", min(data), " max:", max(data)
  print "median: ", CI_50
  print "68% CI: ", CI_68
  print "90% CI: ", CI_90
  print "95% CI: ", CI_95
  print "99% CI: ", CI_99

  if logplot is True:
    ax.set_xscale('log')  
    logbins = 10**linspace(log10(xmin), log10(xmax), nbins)
    norm = 1.0/(len(data)*(log10(xmax)-log10(xmin))/nbins)
    ax.hist(data, bins=logbins, histtype='stepfilled', label=label, weights=norm*ones(len(data)), alpha=0.6)
    #    ax.axvline(log10(CL_95), linestyle=':', color='b', label=str(log10(CL_95))+' @ $95\%$', alpha=0.5)
    #    ax.axvline(log10(CL_99), linestyle=':', color='b', label=str(log10(CL_99))+' @ $99\%$', alpha=0.75)
    #    ax.axvline(log10(lowerbound), linestyle=':', color='b', label='$\lambda_g < '+str(lowerbound)+'$')
  else:
    ax.hist(data, bins=nbins, histtype='stepfilled', label=label, normed=True, alpha=0.6)
    #    ax.axvline(CL_95, linestyle=':', color='b', label=str(CL_95)+' @ $95\%$', alpha=0.5)
    #    ax.axvline(CL_99, linestyle=':', color='b', label=str(CL_99)+' @ $99\%$', alpha=0.75)
    #    ax.axvline(lowerbound, linestyle=':', color='b', label='$\lambda_g < '+str(lowerbound)+'$')
  
  ax.grid(color='grey', linestyle='--', linewidth=0.8, alpha=0.8)
  ax.axvline(CI_90, linewidth=3, color='b', label='90%')
      
  return fig


def plotPosteriorKDE(data, nbins=50, xmin=None, xmax=None, logplot=False, weights=None, ax=None, label="GW"):
  '''Posterior plots with Gaussian KDE'''
  if ax is None:
    fig = figure()
    ax = fig.add_subplot(111)
    if xmin is not None:
      ax.set_xlim(xmin=xmin)
    else:
      xmin = min(data)
    if xmax is not None:
      ax.set_xlim(xmax=xmax)
    else:
      xmax = max(data)
  
  if weights is None:
    kde = stats.gaussian_kde(data)
  else:
    print "Not ready yet"
    exit(-1)
#    kde = kde_weighted.gaussian_kde(data, weights)

  if logplot is True:
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    x = pow(10, linspace(log10(xmin), log10(xmax), 1000))
    ax.semilogx(x, kde.evaluate(x), label=label)
  else:
    x = linspace(xmin, xmax, 1000)
    ax.plot(x, kde.evaluate(x), label=label)
    ax.fill_between(x, zeros(len(x)), kde.evaluate(x), facecolor="b", alpha=0.3)
  ax.grid(color='grey', linestyle='--', linewidth=0.3, alpha=0.8)
  ax.set_ylabel("Probability Density")

  kdemax = max(kde.evaluate(x))
  def kdesolve(x):
    return kde.evaluate(x) - kdemax/2.0
  bound = newton(kdesolve, -22)
  print kdemax, pow(10, bound)


  return fig, kde


if __name__ == "__main__":

###########################################
#
#  Parse arguments (argparser)
#
###########################################

  parser = argparse.ArgumentParser(description="This is a post-processing script for massive graviton parameters.")
#  parser.add_argument("-v", action="store_true", dest="verbose", help="Print detailed output",default=False)
  parser.add_argument("-i", "--input", type=str, dest="datafile", nargs='+', help="file(s) containing data points", metavar="FILE")
  parser.add_argument("-l", "--label", type=str, dest="labels", nargs='+', help="source-identifying string(s)", default=None)
  parser.add_argument("-p", "--params", type=str, dest="paramfile", nargs='+', help="file(s) containing param names per column", metavar="FILE", default=None)
  parser.add_argument("-o", "--output", type=str, dest="outputfolder", help="outputfolder", metavar="FOLDER",default=".")
#  parser.add_argument("--weightchain", action="store_true", dest="weightchain", help="plot KDE using weighted chain points instead of posteriors",default=False)
  parser.add_argument("--no2d", action="store_true", dest="no2d", help="Do not produce 2d summary plots", default=False)
  parser.add_argument("-e", "--eventplot", action="store_true", dest="eventplot", help="plot sample of posterior points (mpl v1.3)", default=False)


  args = parser.parse_args()
  
  datafile = args.datafile
  labels = args.labels
  paramfile = args.paramfile
  outfolder = args.outputfolder
  no2d = args.no2d
  eventplot = args.eventplot

  if not os.path.exists(outfolder):
    os.makedirs(outfolder)
  
  if labels:
    if len(labels)!=len(datafile):
      print "ERROR: need to give same number of datafiles and labels"
      sys.exit(-1)

  if paramfile:
    if len(paramfile)!=len(datafile):
      print "ERROR: numbers of datafiles and paramfiles do not match"
      sys.exit(-1)
      
  kdelist_lg=[]
  kdelist_lamg=[]
  kdelist_mg=[]
  
  for (dfile, lab) in zip(datafile, labels):
    data = genfromtxt(dfile, names=True)

    #  if paramfile is None:
    #    paramfile = datafile+'_params.txt'
    #  paramdict = getParams(paramfile)

    """Converting (log)distance posterior to meters"""
    if "logdistance" in data.dtype.names:
      distdata = exp(data["logdistance"]) * 1e6 * pc_SI
      print "logdist"
    elif "distance" in data.dtype.names:
      distdata = data["distance"] * 1e6 * pc_SI
      print "lindist"
    else:
      print "ERROR: No distance posterior! Exiting..."
      sys.exit(-1)

    """Calculating redshifts"""
    zdata = calculate_redshift(distdata)

    logldata = data["logl"]

    """Calculating posteriors for lambda_eff parameters"""
    if "lambda_eff" in data.dtype.names:
      lgdata = data["lambda_eff"]
      lamgdata = GComptonWavelength(lgdata, zdata)
      mgdata = GravitonMass(lamgdata)
      figl_hist = plotPosteriorHist(lgdata, xlabel="$l_g$", label=lab)
      figl_kde, kde_lg = plotPosteriorKDE(lgdata, xlabel="$l_g$", label=lab)
      figlam_hist = plotPosteriorHist(lamgdata, xlabel="$\lambda_eff$", label=lab)
      figlam_kde, kde_lamg = plotPosteriorKDE(lamgdata, xlabel="$\lambda_eff$", label=lab)
      figm_hist = plotPosteriorHist(mgdata, xlabel="$m_g$", label=lab)
      figm_kde, kde_mg = plotPosteriorKDE(mgdata, xlabel="$m_g$", label=lab)
    elif "loglambda_eff" in data.dtype.names:
      loglgdata = array(data["loglambda_eff"])
      lgdata = pow(10, loglgdata)
      lamgdata = GComptonWavelength(lgdata, zdata)
      loglamgdata = log10(lamgdata)
      savetxt(os.path.join(outfolder, "lamda_g_posteriors.dat"), lamgdata)
      savetxt(os.path.join(outfolder, "loglamda_g_posteriors.dat"), loglamgdata)
      mgdata = GravitonMass(lamgdata)
      logmgdata = log10(mgdata)
    
      print "Plotting lg"
      figl_hist = plotPosteriorHist(lgdata, logplot=True, xmin=1e13, xmax=1e20, label=lab)
      figl_hist.gca().set_xlabel("$l_g \, [m]$")
      figl_kde, kde_lg = plotPosteriorKDE(loglgdata, xmin=13, xmax=20, label=lab)
      figl_kde.gca().set_xlabel("$\log_{10}l_g \, [m]$")
      for lb in lambdag_bounds:
        (bound, lc, ls) = lambdag_bounds[lb]
        figl_hist.gca().axvline(bound, linestyle=ls, linewidth=3, color=lc, label=lb)
        figl_kde.gca().axvline(log10(bound), linestyle=ls, linewidth=3, color=lc, label=lb)
      
      print "Plotting lambdag"
    #    figlam_hist = plotPosteriorHist(lamgdata, logplot=True, xmin=1e13, xmax=1e20, label=lab)
      figlam_hist = plotPosteriorHist(lamgdata, logplot=True, nbins=100, xmin=1e13, xmax=max(lamgdata), label=lab)
      figlam_hist.gca().set_xlabel("$\lambda_eff \, [m]$")
      figlam_kde, kde_lamg = plotPosteriorKDE(loglamgdata, xmin=13, xmax=20, label=lab)
      figlam_kde.gca().set_xlabel("$\log_{10}\lambda_eff \, [m]$")
      for lb in lambdag_bounds:
        (bound, lc, ls) = lambdag_bounds[lb]
        figlam_hist.gca().axvline(bound, linestyle=ls, linewidth=3, color=lc, label=lb)
        figlam_kde.gca().axvline(log10(bound), linestyle=ls, linewidth=3, color=lc, label=lb)

      print "Plotting mg"
      figm_hist = plotPosteriorHist(mgdata, logplot=True, label=lab)
      figm_hist.gca().set_xlabel("$m_g \, [eV/c^2]$")
      figm_kde, kde_mg = plotPosteriorKDE(logmgdata, label=lab)
      figm_kde.gca().set_xlabel("$\log_{10}m_g \, [eV/c^2]$")
      for lb in lambdag_bounds:
        (bound, lc, ls) = lambdag_bounds[lb]
        figm_hist.gca().axvline(GravitonMass(bound), linestyle=ls, linewidth=3, color=lc, label=lb)
        figm_kde.gca().axvline(log10(GravitonMass(bound)), linestyle=ls, linewidth=3, color=lc, label=lb)

    figl_hist.gca().legend(loc='upper left', fancybox=True, framealpha=0.8)
    figl_kde.gca().legend(loc='upper left', fancybox=True, framealpha=0.8)
    figl_hist.savefig(os.path.join(outfolder, "l_g_hist_"+lab+".png"), bbox_inches='tight')
    figl_kde.savefig(os.path.join(outfolder, "l_g_KDE_"+lab+".png"), bbox_inches='tight')

    figlam_hist.gca().legend(loc='upper left', fancybox=True, framealpha=0.8)
    figlam_kde.gca().legend(loc='upper left', fancybox=True, framealpha=0.8)
    figlam_hist.savefig(os.path.join(outfolder, "lambda_eff_hist_"+lab+".png"), bbox_inches='tight')
    figlam_kde.savefig(os.path.join(outfolder, "lambda_eff_KDE_"+lab+".png"), bbox_inches='tight')
    
    figm_hist.gca().legend(loc='upper right', fancybox=True, framealpha=0.8)
    figm_kde.gca().legend(loc='upper right', fancybox=True, framealpha=0.8)
    figm_hist.savefig(os.path.join(outfolder, "m_g_hist_"+lab+".png"), bbox_inches='tight')
    figm_kde.savefig(os.path.join(outfolder, "m_g_KDE_"+lab+".png"), bbox_inches='tight')
    
    
    print "Plotting Likelihood scatter"
    figll = figure()
    axll = figll.add_subplot(111)
    axll.scatter(loglgdata, logldata)
    axll.set_xlabel("$\log\lambda_eff$")
    axll.set_ylabel("$\log L$")
    figll.savefig(os.path.join(outfolder, "lambda_L_"+lab+".png"), bbox_inches='tight')
    
    kdelist_lg.append(kde_lg)
    kdelist_lamg.append(kde_lamg)
    kdelist_mg.append(kde_mg)

  print "Combining posteriors..."
  
  fig_combined_lamg = figure()
  ax_combined_lamg = fig_combined_lamg.add_subplot(111) 
  fig_combined_lamgp = figure()
  ax_combined_lamgp = fig_combined_lamgp.add_subplot(111) 
  logplot=False
  xmin=13
  xmax=20
  ax_combined_lamg.set_xlim(xmin=xmin)
  ax_combined_lamg.set_xlim(xmax=xmax)
  ax_combined_lamgp.set_xlim(xmin=xmin)
  ax_combined_lamgp.set_xlim(xmax=xmax)
  x = linspace(xmin, xmax, 1000)

  lamgprior = lambda l,lmin,lmax: 1.0/(lmax-lmin)
  lamgpriormg = lambda l,lmin,lmax: (lmax-lmin)/(l*lmax*lmin)

  # Prior uniform in loglambda_eff
  y = lamgprior(x, xmin, xmax)
  # Prior uniform in m_g
  yp = lamgpriormg(pow(10.0,x), pow(10.0,xmin), pow(10.0,xmax))

  for ki in kdelist_lamg:
    y *= ki.evaluate(x)
    yp *= ki.evaluate(x)
    
  lnorm = trapz(y,x)
  y = y/lnorm
  lpnorm = trapz(yp,x)
  yp = yp/lpnorm

  ax_combined_lamg.plot(x, y, label="GW150914 + GW151226 (PDF)")
  ax_combined_lamgp.plot(x, yp, label="GW150914 + GW151226 (PDF)")
  ax_combined_lamg.fill_between(x, zeros(len(x)), y, facecolor="b", alpha=0.3)
  ax_combined_lamgp.fill_between(x, zeros(len(x)), yp, facecolor="b", alpha=0.3)
  ax_combined_lamg.grid(color='grey', linestyle='--', linewidth=0.3, alpha=0.8)
  ax_combined_lamgp.grid(color='grey', linestyle='--', linewidth=0.3, alpha=0.8)
  ax_combined_lamg.set_ylabel("PDF / CDF")
  ax_combined_lamgp.set_ylabel("PDF / CDF")
  ax_combined_lamg.set_xlabel("$\log_{10}\lambda_eff \, [m]$")
  ax_combined_lamgp.set_xlabel("$\log_{10}\lambda_eff \, [m]$")

  kdemax = max(y)
  from scipy.integrate import cumtrapz
  ycum = cumtrapz(y,x)
  ypcum = cumtrapz(yp,x)

  ax_combined_lamg.plot(x, hstack((array([0]),ycum)), color='k', label="CDF")
  ax_combined_lamgp.plot(x, hstack((array([0]),ypcum)), color='k', label="CDF")
  ax_combined_lamg.fill_between(x, zeros(len(x)), hstack((array([0]),ycum)), facecolor="k", alpha=0.2)
  ax_combined_lamgp.fill_between(x, zeros(len(x)), hstack((array([0]),ypcum)), facecolor="k", alpha=0.2)

  bounds = []
  for lb in lambdag_bounds:
    (bound, lc, ls) = lambdag_bounds[lb]
    fig_combined_lamg.gca().axvline(log10(bound), linestyle=ls, linewidth=3, color=lc, label=lb)
    fig_combined_lamgp.gca().axvline(log10(bound), linestyle=ls, linewidth=3, color=lc, label=lb)
    bounds.append(log10(bound))
  print where(ypcum<0.1)[0][-1], where(ypcum>0.1)[0][0], x[where(ypcum<0.1)[0][-1]], x[where(ypcum>0.1)[0][0]]
  combined_90bound = x[where(ypcum>0.1)[0][0]]
  bounds.append(combined_90bound)
  ax_combined_lamgp.errorbar(bounds, 0.5*ones(len(lambdag_bounds)+1), xerr=0.3, color='k', linewidth=3,  xlolims=True)
  fig_combined_lamgp.gca().axvline(combined_90bound, linewidth=3, color='darkblue', label="90% lower bound")

  fig_combined_lamg.gca().legend(loc='upper left', fancybox=True, framealpha=0.8)
  fig_combined_lamgp.gca().legend(loc='upper left', fancybox=True, framealpha=0.8)
  fig_combined_lamg.savefig(os.path.join(outfolder, "combined_lambda_eff.png"), bbox_inches='tight')
  fig_combined_lamgp.savefig(os.path.join(outfolder, "combined_lambda_eff_mgprior.png"), bbox_inches='tight')


  print "DONE!"
