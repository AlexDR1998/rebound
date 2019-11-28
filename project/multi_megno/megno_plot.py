import numpy as np
#import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#Graphing half of python_examples/megno/problem.py. Split so that graphing did not
#occur on computing cluster
#
#Comment/uncomment a and e arrays for wide/narrow ranges of values, change
#.txt files loaded for different plots

N =200                      # Grid size, increase this number to see more detail
#a = np.linspace(0.97,0.98,N)   # range of cruithne semi-major axis in AU
#e = np.linspace(0.51,0.52,N)   # range of cruithne eccentricity
a = np.linspace(0.8,1.2,N)   # wider range
e = np.linspace(0.3,0.7,N)  
megno = np.load("no_jupiter/wide/megno2.npy")
lyaptimescale = np.load("no_jupiter/wide/lyaptimescale2.npy")
#lyaptimescale = 1E6/(2*megno)

extent = [a.min(), a.max(), e.min(), e.max()]
plt.imshow(megno, vmin=1.8, vmax=20., aspect='auto', origin="lower", interpolation='nearest', cmap="plasma", extent=extent)
plt.colorbar(label="MEGNO")
plt.xlabel("Semi Major Axis (AU)")
plt.ylabel("Eccentricity")
plt.show()
plt.imshow(lyaptimescale, vmin=1e1, vmax=1e5, norm=LogNorm(), aspect='auto', origin="lower", interpolation='nearest', cmap="viridis", extent=extent)
plt.colorbar(label="Lyapunov Timescale (years)")
plt.xlabel("Semi Major Axis (AU)")
plt.ylabel("Eccentricity")
plt.show()
"""
f, axarr = plt.subplots(2,figsize=(10,10))
extent = [a.min(), a.max(), e.min(), e.max()]
for ax in axarr:
    ax.set_xlim(extent[0],extent[1])
    ax.set_ylim(extent[2],extent[3])
    ax.set_xlabel("$a_{\mathrm{Cruithne}}$ [AU]")
    ax.set_ylabel("$e_{\mathrm{Cruithne}}$")

# Plot MEGNO 
im1 = axarr[0].imshow(megno, vmin=1.8, vmax=4., aspect='auto', origin="lower", interpolation='nearest', cmap="RdYlGn_r", extent=extent)
cb1 = plt.colorbar(im1, ax=axarr[0])
cb1.solids.set_rasterized(True)
cb1.set_label("MEGNO $\\langle Y \\rangle$")

# Plot Lyapunov timescale
im2 = axarr[1].imshow(lyaptimescale, vmin=1e1, vmax=1e5, norm=LogNorm(), aspect='auto', origin="lower", interpolation='nearest', cmap="RdYlGn", extent=extent)
cb2 = plt.colorbar(im2, ax=axarr[1])
cb2.solids.set_rasterized(True)
cb2.set_label("Lyapunov timescale [years]")

plt.savefig("megno_wide_no_jupiter.pdf")
"""