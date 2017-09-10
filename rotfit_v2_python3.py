import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
import scipy.integrate as integrate
import scipy.optimize as opt
from tkinter import Tk
import locale
import sys

plt.style.use('classic')
plt.close()

# conversion factors
arcsec2rad = np.pi / 180.0 / 3600.0
mpc2kpc = 1000
kpc2pc = 1000
pc2m = 3.1e16
km2m = 1000
kg2solarm = 5e-31

epsilon = 0.1   #prevents division by 0

# definition of the rotation curve
def rotvel(r,rc,vinf):
	x = (r + epsilon) / rc
	aux = 1.0 - np.arctan(x) / x
	result = vinf * np.sqrt(aux)
	return(result)

# defintion of mass integrand
def massfunc(r, rc, p0):
   return(4*  np.pi * r**2 * p0 / (1+(r/rc)**2))

# copy results to clipboard function
def clip(self):
    locale.setlocale(locale.LC_ALL, '')
    r = Tk()
    r.withdraw()
    r.clipboard_clear()
    r.clipboard_append(locale.format("%.1f", distance / mpc2kpc, grouping=False) + '\t' + locale.format("%.0f", vinf, grouping=False) + '\t' + locale.format("%.2f", rc, grouping=False) + '\t' + locale.format("%.2f", p0, grouping=False) + '\t' + locale.format("%.2E", M, grouping=False))
    r.update() # now it stays on the clipboard after the window is closed
    r.destroy()

# data: radius in arcsec, velocity in km/s and error in km/s
while True:
    try:
        dist_data = np.loadtxt('../data/summary.data', dtype={'names':('name','distance'),'formats':('S10','float32')}) #returns data as ndarray
    except:
        print('Could not open summary.data. Check directories. ' + str(sys.exc_info()[0]))
        input('Press enter to quit.')
    print('Found ' + str(len(dist_data)) + ' entries in summary.data:')
    i = 0
    while i < len(dist_data):
        print(dist_data[i][0].decode('utf-8'))
        i = i + 1
    filename = input('Load entry: ')
    try:
        i = 0
        while dist_data[i][0] != filename.encode('utf-8'):       # search for the correct entry in summary.data   
                i = i + 1
        distance = dist_data[i][1] * mpc2kpc
        data = np.loadtxt('../data/' + filename + '.data')
        r = abs(data[:,0]) * arcsec2rad * distance  # radial distance
        v = abs(data[:,1])                          # orbital speed
        e = data[:,2]                               # orbital speed error
        break
    except:
        print('Could not open. Check input. ' + str(sys.exc_info()[0]))

# plot data
plt.errorbar(r,v,yerr=e,fmt='ro',label="data")

# fit to data, and print out distance, rc and vinf
guess = np.array([5.0, 160.0])                               
param, covar = opt.curve_fit(rotvel,r,v,guess,e)
rc = param[0]
vinf = param[1]
print('dist: ' + str(distance))     #Fit of the function rotvel to r und v with the error e
print('rc: ' + str(rc))            # and a guess for the parameters.
print('vinf: ' + str(vinf))        # Print calculated parameters (param) and covariance (covar).

# calculate p0
p0 = 1/(4* np.pi * 6.7e-11) * vinf**2/rc**2 * km2m**2 * pc2m * kg2solarm / kpc2pc**2   #conversion to Msol/pc^3

# calculate mass
I = integrate.quad(massfunc, 0, 10 * rc * kpc2pc, args=(rc * kpc2pc , p0))   # calculate mass inside the tidal radius
M = I[0]

# plot fit
rr = np.linspace(min(r),max(r),1001)
v_fit = rotvel(rr, rc, vinf)
plt.plot(rr,v_fit,'b-',label="best-fit rotation curve")

plt.xlabel(r'radius $r$ in $\mathrm{kpc}/h$')
plt.ylabel(r'velocity $\upsilon$ in $\mathrm{km}/\mathrm{s}$')

plt.legend(loc="lower right")

plt.title(filename)

plt.gca().set_position((.1, .3, .8, .6))   #sets the position of the axes to .1, .25 with the width .8 and height .6

plt.figtext(0.05,0.16,r'$d = ' + "{:.1f}".format(distance / mpc2kpc) + '$ $Mpc$')
plt.figtext(0.05,0.13,r'$v_\inf = ' + "{:.0f}".format(vinf) + '$ $km/s$')
plt.figtext(0.05,0.1,r'$r_c = ' + "{:.2f}".format(rc) + '$ $kpc/h$')
plt.figtext(0.05,0.07, r'$p_0 = ' + "{:.2f}".format(p0) + '$ $M_\odot/pc^3$')
plt.figtext(0.05,0.04, r'$M = ' + "{:.2E}".format(M) + '$ $M_\odot$')

axclip = plt.axes([.4, .05, .2, .1])    # set position and dimensions of button
bclip = widgets.Button(axclip, 'Copy Results')  # create button
bclip.on_clicked(clip)  # call function clip

plt.show()