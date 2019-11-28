# Import the rebound module
import sys
#import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm
import rebound
import numpy as np
import time
from rebound.interruptible_pool import InterruptiblePool
import warnings

def simulation(par):
    integrator, run, trial = par
    sim = rebound.Simulation()
    k = 0.01720209895    
    Gfac = 1./k
    #Gfac = 1
    sim.dt = dt
    if integrator == "whfast-nocor":
        integrator = "whfast"
    else:
        sim.ri_whfast.corrector = 17
    sim.integrator = integrator 
    sim.ri_whfast.safe_mode = 0

    #massfac = 1.
    #sim.add(m=1.00000597682, x=-4.06428567034226e-3, y=-6.08813756435987e-3, z=-1.66162304225834e-6,      vx=+6.69048890636161e-6*Gfac, vy=-6.33922479583593e-6*Gfac, vz=-3.13202145590767e-9*Gfac)   # Sun
    #sim.add(m=massfac/1407.355,   x=+3.40546614227466e+0, y=+3.62978190075864e+0, z=+3.42386261766577e-2, vx=-5.59797969310664e-3*Gfac, vy=+5.51815399480116e-3*Gfac, vz=-2.66711392865591e-6*Gfac)   # Jupiter
    #sim.add(m=massfac/3501.6,     x=+6.60801554403466e+0, y=+6.38084674585064e+0, z=-1.36145963724542e-1, vx=-4.17354020307064e-3*Gfac, vy=+3.99723751748116e-3*Gfac, vz=+1.67206320571441e-5*Gfac)   # Saturn
    #sim.add(m=massfac/22869.,     x=+1.11636331405597e+1, y=+1.60373479057256e+1, z=+3.61783279369958e-1, vx=-3.25884806151064e-3*Gfac, vy=+2.06438412905916e-3*Gfac, vz=-2.17699042180559e-5*Gfac)   # Uranus
    #sim.add(m=massfac/19314.,     x=-3.01777243405203e+1, y=+1.91155314998064e+0, z=-1.53887595621042e-1, vx=-2.17471785045538e-4*Gfac, vy=-3.11361111025884e-3*Gfac, vz=+3.58344705491441e-5*Gfac)   # Neptune
    #sim.G = 1.4880826e-34 #Units of AU^3/kg/day^2
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
    for x in [0,1,2,3,4,5]:#,6,7,8]:#,9,10]:
    #for x in [0,6,7]:
        p = rebound.Particle(m=ss_mass[x]/ss_mass[0],x=ss_pos[x,0],y=ss_pos[x,1],z=ss_pos[x,2],vx=ss_vel[x,0]*Gfac,vy=ss_vel[x,1]*Gfac,vz=ss_vel[x,2]*Gfac)
        #p = rebound.Particle(m=ss_mass[x],x=ss_pos[x,0],y=ss_pos[x,1],z=ss_pos[x,2],vx=ss_vel[x,0],vy=ss_vel[x,1],vz=ss_vel[x,2])
        sim.add(p)


    N = sim.N
    particles = sim.particles
    np.random.seed(run)
    for p in particles:
        p.m *= 1.+1e-3*np.random.rand()
        p.x *= 1.+1e-3*np.random.rand()
        p.y *= 1.+1e-3*np.random.rand()
        p.z *= 1.+1e-3*np.random.rand()
        p.vx *= 1.+1e-3*np.random.rand()
        p.vy *= 1.+1e-3*np.random.rand()
        p.vz *= 1.+1e-3*np.random.rand()

    def move_to_heliocentric():
        particles[0].x  = 0.
        particles[0].y  = 0. 
        particles[0].z  = 0. 
        particles[0].vx = 0. 
        particles[0].vy = 0. 
        particles[0].vz = 0. 


    def energy():
        com_vx = 0.
        com_vy = 0.
        com_vz = 0.
        if integrator=="mercury" or integrator[0:7]=="swifter":
            mtot = 0.
            for p in particles:
                com_vx += p.vx*p.m 
                com_vy += p.vy*p.m 
                com_vz += p.vz*p.m 
                mtot += p.m
            com_vx /= mtot
            com_vy /= mtot
            com_vz /= mtot
        E_kin = 0.
        E_pot = 0.
        for i in xrange(N):
            #Kinetic energy per particle
            dvx = particles[i].vx - com_vx
            dvy = particles[i].vy - com_vy
            dvz = particles[i].vz - com_vz
            E_kin += 0.5*particles[i].m*(dvx*dvx + dvy*dvy + dvz*dvz)
            for j in xrange(i+1,N):
                #Potential energy felt between all particles
                dx = particles[i].x-particles[j].x
                dy = particles[i].y-particles[j].y
                dz = particles[i].z-particles[j].z
                r2 = dx*dx + dy*dy + dz*dz
                E_pot -= particles[i].m*particles[j].m/np.sqrt(r2)
                #E_pot -= sim.G*particles[i].m*particles[j].m/np.sqrt(r2)
        return E_kin+E_pot




    def momentum():
        com_vx = 0.
        com_vy = 0.
        com_vz = 0.
        if integrator=="mercury" or integrator[0:7]=="swifter":
            mtot = 0.
            for p in particles:
                com_vx += p.vx*p.m 
                com_vy += p.vy*p.m 
                com_vz += p.vz*p.m 
                mtot += p.m
            com_vx /= mtot
            com_vy /= mtot
            com_vz /= mtot
        lin_mom = np.array([0.,0.,0.])
        ang_mom = np.array([0.,0.,0.])
        for i in xrange(N):
            #Velocities of all particles
            dvx = particles[i].vx - com_vx
            dvy = particles[i].vy - com_vy
            dvz = particles[i].vz - com_vz
            #Linear momentum
            lin_mom_temp=particles[i].m*np.array([dvx,dvy,dvz])
            lin_mom+=lin_mom_temp
            
            #Angular momentum
            r_sun = np.array([particles[0].x,particles[0].y,particles[0].z])
            ang_mom+=np.cross(r_sun,lin_mom_temp)

            #E_kin += 0.5*particles[i].m*(dvx*dvx + dvy*dvy + dvz*dvz)
            
            #for j in xrange(i+1,N):
                #Potential energy felt between all particles
            #    dx = particles[i].x-particles[j].x
            #    dy = particles[i].y-particles[j].y
            #    dz = particles[i].z-particles[j].z
            #    r2 = dx*dx + dy*dy + dz*dz
            #    E_pot -= particles[i].m*particles[j].m/np.sqrt(r2)
                #E_pot -= sim.G*particles[i].m*particles[j].m/np.sqrt(r2)
        return lin_mom,ang_mom

    

    times = np.logspace(np.log10(orbit),np.log10(tmax),Ngrid)
    if integrator=="mercury" or integrator[0:7]=="swifter":
        move_to_heliocentric()
    else:
        sim.move_to_com()
    ei = energy()
    lm_i,am_i=momentum()

    mag_lm_i = np.linalg.norm(lm_i)
    mag_am_i = np.linalg.norm(am_i)
    es = []
    lms = []
    ams = []

    runtime = 0.
    start = time.time()
    # Capture warning messages (WHFast timestep too large)
    with warnings.catch_warnings(record=True) as w: 
        warnings.simplefilter("always")
        for t in times:
            sim.integrate(t,exact_finish_time=0)
            ef = energy()
            lm,am=momentum()
            e = np.fabs((ei-ef)/ei)+1.1e-16
            l = np.linalg.norm(lm-lm_i)#/mag_lm_i
            a = np.linalg.norm(am-am_i)#/mag_am_i
            es.append(e)
            lms.append(l)
            ams.append(a)
    
    integrator, run, trial = par
    print(integrator.ljust(13) + " %9.5fs"%(time.time()-start) + "\t Energy error: %e"  %( e) + "\t Linear momentum error: %e" %( l) +"\t Angular momentum error: %e" %( a))
    
    es = np.array(es)
    lms = np.array(lms)
    ams = np.array(ams)
    return [times, es,lms,ams]

Ngrid = 500
#3dt = 100.23
orbit = 11.8618*1.*np.pi
dt = orbit/3000.
tmax = orbit*1e2        # Maximum integration time.
#integrators = ["whfast-nocor", "whfast"]
integrators = ["whfast","ias15","janus"]
#integrators = ["mercury","swifter-whm","whfast-nocor", "whfast"]
colors = {
    'whfast-nocor': "#00AA00",
    'whfast':       "#FF0000",
    'mercury':      "#6E6E6E",
    'swifter-whm':  "#444444",
    'swifter-helio':"#AABBBB",
    'swifter-tu4':  "#FFAAAA",
    'ias15':        "g",
    'janus':        "#0000FF",
}
trials = 4
    
parameters = [(inte,i*trials+j,j) for i,inte in enumerate(integrators) for j in xrange(trials)]
if len(sys.argv)!=2:
    pool = InterruptiblePool()
    print("Running %d simulations" % (len(parameters)))
    res = np.array(pool.map(simulation,parameters)).reshape(len(integrators),trials,4,Ngrid)
    np.save("res.npy",res)
else:
    print("Loading %d simulations" % (len(parameters)))
    print(sys.argv[1])
    res = np.load(sys.argv[1])



f,axarr = plt.subplots(1,1,figsize=(13,4))
#extent=[res[:,:,0,:].min()/orbit, res[:,:,0,:].max()/orbit, 1e-16, 1e-9]
extent=[res[:,:,0,:].min()/orbit, res[:,:,0,:].max()/orbit, 1e-24, 1e-7]

axarr.set_xlim(extent[0], extent[1])
axarr.set_ylim(extent[2], extent[3])
axarr.set_xlabel(r"time [Jupiter orbits]")
axarr.set_ylabel(r"linear momentum error")
plt.xscale('log')
plt.yscale('log')
plt.grid(True)


res_mean = np.mean(res,axis=1)
for i in xrange(len(res)):
    for j in xrange(trials):
        #res_trial = res[i,j,:2,:]
        res_trial = res[i,j,[0,2],:]
        im1 = axarr.plot(res_trial[0]/orbit,res_trial[1], color=colors[integrators[i]],alpha=0.2)
    im1 = axarr.plot(res_mean[i][0]/orbit,res_mean[i][2], label=integrators[i].upper(),color=colors[integrators[i]], linewidth=2.0)

from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
#lgd = plt.legend(loc="upper center",  bbox_to_anchor=(0.5, -0.2),  prop = fontP,ncol=3,frameon=False, numpoints=1, scatterpoints=1 , handletextpad = 0.2, markerscale=2.)
plt.legend(loc="upper left")
plt.show()

#plt.savefig("longtermtest_angmom.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')

#from sys import platform as _platform
#if _platform == "darwin":
#    import os
#    os.system("open longtermtest.pdf")
