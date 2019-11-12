#!/usr/bin/python
# This example integrates Jupiter and Saturn in the Solar system for a variety of initial conditions.
# Alongside the normal equations of motions, IAS15 is used to integrate the variational equations.
# These can be used to measure the Mean Exponential Growth of Nearby Orbits (MEGNO), a chaos indicator.
# This example script runs 12^2 simulations and plots the MEGNO value. Values close to <Y>=2 correspond 
# to regular quasi-periodic orbits. Higher values of <Y> correspond to chaotic orbits.

# Import matplotlib
import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
 
# Import the rebound module
import rebound
# Import other modules
import numpy as np
import multiprocessing
import warnings
import time
# Runs one simulation.
def simulation(par):
    cruithne_a, cruithne_e = par
    sim = rebound.Simulation() 
    sim.integrator = "whfast"
    sim.min_dt = 0.05
    sim.dt = 1
    sim.G = 1.4880826e-34 #Units of AU^3/kg/day^2
    ss_pos = np.array([[-3.013424684820004E-03 , 7.577400311699641E-03 , 1.522327948994377E-06],
                     [-1.744888637875314E-01 ,-4.244459073219732E-01 ,-1.957052267164663E-02],
                     [-5.832166811075774E-01 ,-4.227197431546666E-01 , 2.757923523005813E-02],
                     [ 9.926803557750922E-01 , 1.171071592931596E-01 ,-5.897765002577186E-06],
                     [ 8.228862161210234E-01 , 6.587996475726317E-01 ,-3.785334833271938E-01],
                     [-1.644137782803167E+00 , 2.530657764233049E-01 , 4.541188790127934E-02],
                     [-1.699642400881444E-01 ,-5.250743091389841E+00 , 2.557799438953158E-02],
                     [ 3.337284978938324E+00 ,-9.463581115929795E+00 , 3.169757971579321E-02],
                     [ 1.643185373589321E+01 , 1.110227970027843E+01 ,-1.716425425492940E-01],
                     [ 2.917750633262754E+01 ,-6.646446374100315E+00 ,-5.355542983258703E-01],
                     [ 1.269610583507422E+01 ,-3.141113199578617E+01 ,-3.112775015648088E-01]])
    ss_vel = np.array([[-8.492322341632166E-06 ,-9.353155515358904E-07 , 2.301541419735291E-07],
                     [ 2.048459325443986E-02 ,-8.995044409873299E-03 ,-2.614664569098718E-03],
                     [ 1.189797004491595E-02 ,-1.634153795167905E-02 ,-9.110754924475504E-04],
                     [-2.170315646338148E-03 , 1.703520105098581E-02 ,-5.463857304374388E-07],
                     [-1.381007577874730E-02 , 5.897333087115376E-03 , 2.754466583425123E-03],
                     [-1.556986377304836E-03 ,-1.264431146517457E-02 ,-2.267130777538514E-04],
                     [ 7.450947804402543E-03 , 1.166544750377484E-04 ,-1.671012875749180E-04],
                     [ 4.952121119936067E-03 , 1.839073137038874E-03 ,-2.293132844397104E-04],
                     [-2.230759206307085E-03 , 3.075630739324861E-03 , 4.027037883636828E-05],
                     [ 6.759177587143499E-04 , 3.079179855664010E-03 ,-7.882476544965271E-05],
                     [ 2.988188173523851E-03 , 5.096901737398172E-04 ,-9.289666940024388E-04]])

    ss_mass = np.array([1.988544e30,3.302e23,48.685e23,6.045476731e24,1.3e14,6.4185e23,1898.13e24,5.68319e26,86.8103e24,102.41e24,1.4639248e+22])

    # These parameters are only approximately those of Jupiter and Saturn.
    #Add particles
    for x in [0,1,2,3]:#,5,6,7,8,9,10]:
        p = rebound.Particle(m=ss_mass[x],x=ss_pos[x,0],y=ss_pos[x,1],z=ss_pos[x,2],vx=ss_vel[x,0],vy=ss_vel[x,1],vz=ss_vel[x,2])
        sim.add(p)

    #cruithne = rebound.Particle()
    #                       M=0.6989875253,
    #                   M=4.97689734325,
    #                   Omega=2.205628636,
    #                   omega=0.7616728033,
    cruithne = sim.add(m=ss_mass[4],
                       a=cruithne_a,
                       e=cruithne_e,
                       inc=0.34567242532, 
                       Omega=2.20306437947,
                       omega=0.76510954706,
                       f=3.98060732223)
    
    for x in [5,6,7,8,9,10]:
        p = rebound.Particle(m=ss_mass[x],x=ss_pos[x,0],y=ss_pos[x,1],z=ss_pos[x,2],vx=ss_vel[x,0],vy=ss_vel[x,1],vz=ss_vel[x,2])
        sim.add(p)
    


    #sun     = rebound.Particle(m=1.)
    #sim.add(sun)
    #jupiter = sim.add(primary=sun,m=0.000954, a=5.204, M=0.600, omega=0.257, e=0.048)
    #saturn  = sim.add(primary=sun,m=0.000285, a=saturn_a, M=0.871, omega=1.616, e=saturn_e)

    sim.move_to_com()
    sim.init_megno()
    # Hide warning messages (WHFast timestep too large)
    with warnings.catch_warnings(record=True) as w: 
        warnings.simplefilter("always")
        sim.integrate(1E6)

    return [sim.calculate_megno(),1./(sim.calculate_lyapunov()*2.*np.pi)] # returns MEGNO and Lypunov timescale in years


### Setup grid and run many simulations in parallel
# actual a=9.977049427315410E-01
# actual e=5.147985107081943E-01
N = 20                      # Grid size, increase this number to see more detail
a = np.linspace(0.97,0.98,N)   # range of cruithne semi-major axis in AU
e = np.linspace(0.51,0.52,N)   # range of cruithne eccentricity
#a = np.linspace(0.8,1.2,N)
#e = np.linspace(0.4,0.6,N)
parameters = []
for _e in e:
    for _a in a:
        parameters.append([_a,_e])
#print("Running initial simulation")
#simulation((1,0.5))

t1 = time.time()
# Run simulations in parallel
pool = rebound.InterruptiblePool()    # Number of threads default to the number of CPUs on the system
print("Running %d simulations on %d threads..." % (len(parameters), pool._processes))
res = np.nan_to_num(np.array(pool.map(simulation,parameters))) 
megno = np.clip(res[:,0].reshape((N,N)),1.8,4.)             # clip arrays to plot saturated 
lyaptimescale = np.clip(np.absolute(res[:,1].reshape((N,N))),1e1,1e5)
t2 = time.time()
print("Time taken: "+str(t2-t1))
### Create plot and save as pdf 

# Setup plots
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

plt.savefig("megno.pdf")

### Automatically open plot (OSX only)
from sys import platform as _platform
if _platform == "darwin":
    import os
    os.system("open megno.pdf")
