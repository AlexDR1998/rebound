import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
#plt.rcParams['axes.facecolor'] = 'black'
import sys


# --- Lots of plotting and graphs for Cruithne orbital elements
#
#
#

def cruithne_rotating_frame(sun,earth,cruithne):
	#Function to plot orbit of sun, earth, cruithne in a reference frame orbiting sun with earth.
	#i.e. earth appears stationary, see orbit of cruithne from earth's pov.


	#s = 10000
	s = 0
	#e = 8000
	e = 2000
	#e = 40000
	sun = sun[s:e]
	earth = earth[s:e]
	cruithne = cruithne[s:e]
	
	N = earth.shape[0]

	
	#Find angle earth is at
	thetas = np.arctan2(earth[:,1],earth[:,0])
	r = np.array([[np.cos(thetas) , np.sin(thetas)],
				  [-np.sin(thetas), np.cos(thetas)]])

	for x in range(N):
		#Apply rotation
		earth[x,0:2] = np.matmul(r[:,:,x],earth[x,0:2])
		cruithne[x,0:2] = np.matmul(r[:,:,x],cruithne[x,0:2])


	S = 36
	#Average over year
	#f = lambda x:np.ones(x)/float(x)
	#cruithne_smooth = cruithne[:]
	#cruithne_smooth[:,0] = np.convolve(cruithne[:,0],f(S),"same")
	#cruithne_smooth[:,1] = np.convolve(cruithne[:,1],f(S),"same")

	plt.plot(sun[:,0],sun[:,1],label="Sun")
	plt.plot(earth[:,0],earth[:,1],label="Earth")
	plt.plot(cruithne[s:e,0],cruithne[s:e,1],label="Cruithne")
	#plt.scatter(cruithne_smooth[::S,0],cruithne_smooth[::S,1])
	plt.legend()
	plt.xlabel("AU")
	plt.ylabel("AU")
	plt.show()



def fft(sun,earth,cruithne):
	#--- Plots fft of cruithne/earth positions. Obsolete, not really used
	earth_fft = np.fft.fftn(earth)
	cruithne_fft = np.fft.fftn(cruithne)
	plt.plot((np.abs(earth_fft[:,0])))
	plt.plot((np.abs(earth_fft[:,1])))
	#plt.plot(cruithne_fft[:,0].imag)
	

	plt.show()
	print(cruithne_fft.shape)

def earth_cruithne_distance(earth,cruithne):
	#--- Plots distance between earth and cruithne

	dist = np.linalg.norm(earth-cruithne,axis=1)
	print(dist.shape)
	
	#Need to filter out higher frequency data (yearly fluctuation)
	b,a = signal.butter(3,0.002)
	dist_smooth = signal.filtfilt(b,a,dist)

	plt.plot(dist)
	plt.show()

	peaks,_=signal.find_peaks(dist_smooth)
	plt.plot(dist_smooth)
	plt.plot(peaks,dist_smooth[peaks],"x",color="orange")
	print(np.diff(peaks)*100.0/356.0)
	plt.show()
	

	#t = np.arange(dist_smooth.shape[0])
	#freq = np.fft.fftfreq(dist_smooth.shape[0])
	ft = np.abs(np.fft.fft(dist))
	#plt.plot(freq,ft)
	plt.plot(ft)
	plt.show()


def cruithne_freq_corot(earth,cruithne):
	#--- Fft of cruithne in corotating reference frame with earth. Not really used in the end

	N = earth.shape[0]
	thetas = np.arctan2(earth[:,1],earth[:,0])
	r = np.array([[np.cos(thetas) , np.sin(thetas)],
				  [-np.sin(thetas), np.cos(thetas)]])

	for x in range(N):
		#Apply rotation
		earth[x,0:2] = np.matmul(r[:,:,x],earth[x,0:2])
		cruithne[x,0:2] = np.matmul(r[:,:,x],cruithne[x,0:2])

	#smooth data
	b,a = signal.butter(3,0.002)
	cruithne_smooth = signal.filtfilt(b,a,cruithne[:,0])

	#peak finding on smooth data
	peaks,_=signal.find_peaks(cruithne_smooth) 

	f = np.abs(np.fft.fft(cruithne[:,1]))
	f_smooth = np.abs(np.fft.fft(cruithne_smooth))
	freq = np.fft.fftfreq(N)




	plt.plot(cruithne[:,0])
	plt.show()
	plt.plot(cruithne_smooth)
	plt.plot(peaks,cruithne_smooth[peaks],"x",color="orange")
	print(np.diff(peaks)*100.0/356.0)
	plt.show()
	plt.plot(freq,f)
	plt.plot(freq,f_smooth)
	plt.show()


def earth_cruithne_freq(earth,cruithne):
	#--- Fft of earth and cruithne positions in ecliptic reference frame

	plt.plot(np.arange(earth.shape[0])*10/365,earth[:,0])
	peaks,_=signal.find_peaks(earth[:,0])
	print(np.mean(np.diff(peaks)*10))
	#plt.plot(cruithne[:,0])
	plt.show()
	freq = np.fft.fftfreq(earth.shape[0])
	ft_e = np.abs(np.fft.fft(earth[:,0]))
	ft_c = np.abs(np.fft.fft(cruithne[:,0]))
	#plt.plot(ft)

	plt.plot(freq,ft_e)
	plt.plot(freq,ft_c)
	plt.show()


def cruithne_semi_major():
	#--- Semi major axis of cruithne. Change .txt file to change which results

	data = np.loadtxt("cruithne_semi_major_axis_janus_rev_short.txt")
	xs = np.arange(data.shape[0])*10/365+2019

	b,a = signal.butter(3,0.002)
	data_smooth = signal.filtfilt(b,a,data)
	peaks,_=signal.find_peaks(data_smooth,height=1.0024)
	times = np.diff(peaks)*10/365
	print(np.mean(times[:3]))
	print(np.std(times[:3]))
	print(np.mean(times[:5]))
	print(np.std(times[:5]))
	plt.plot(xs,data)
	plt.plot(xs[peaks],data[peaks],"x",color="orange")
	plt.title("Fluctuation of semi major axis")
	plt.xlabel("Year")
	plt.ylabel("Semi major axis (AU)")
	#plt.axvline(7000,color="red")
	plt.show()

def cruithne_semi_major_multiplot():
	#--- Semi major axis for all three methods

	data_whf = np.loadtxt("cruithne_semi_major_axis_whfast_rev.txt")
	data_ias = np.loadtxt("cruithne_semi_major_axis_ias15_rev.txt")
	data_jan = np.loadtxt("cruithne_semi_major_axis_janus_rev_ord8.txt")
	xs_whf = np.arange(data_whf.shape[0])*10/365+2019
	xs_ias = np.arange(data_ias.shape[0])*10/365+2019
	xs_jan = np.arange(data_jan.shape[0])*10/365+2019
	plt.plot(xs_whf,data_whf,label="WHFast",color="red")
	plt.plot(xs_ias,data_ias,label="IAS15",color="green")
	plt.plot(xs_jan,data_jan,label="JANUS (8th order)",color="blue")
	#plt.title("Fluctuation of semi major axis")
	plt.xlabel("Year")
	plt.ylabel("Semi major axis (AU)")
	plt.axvline(250000*10/365+2019,color="orange")
	plt.legend()
	plt.show()



def single_reverse_comparison():
	#--- Compares forwards and backwards simulations around a velocity flip, for single orbital element

	data = np.loadtxt("cruithne_semi_major_axis_ias15_rev.txt")
	N = int(data.shape[0]/2)
	ys1 = np.zeros(N)
	ys2 = np.zeros(N)
	for x in range(N):
		ys1[x]=data[N-x]
		ys2[x]=data[N+x]

	plt.plot(ys1)
	plt.plot(ys2)
	plt.show()


def cruithne_reversed_sma():
	#--- Plots how time reversbile semi major axes are for all methods

	data_whf = np.loadtxt("cruithne_semi_major_axis_whfast_rev.txt")
	data_ias = np.loadtxt("cruithne_semi_major_axis_ias15_rev.txt")
	data_jan = np.loadtxt("cruithne_semi_major_axis_janus_rev_ord8.txt")
	ys_whf = np.zeros(int(data_whf.shape[0]/2))
	ys_ias = np.zeros(int(data_ias.shape[0]/2))
	ys_jan = np.zeros(int(data_jan.shape[0]/2))
	N = int(data_whf.shape[0]/2)

	for x in range(ys_whf.shape[0]):
		ys_whf[x]=data_whf[N-x]-data_whf[N+x]
		ys_ias[x]=data_ias[N-x]-data_ias[N+x]
		ys_jan[x]=data_jan[N-x]-data_jan[N+x]

	xs_whf = np.arange(ys_whf.shape[0])*10/365
	xs_ias = np.arange(ys_ias.shape[0])*10/365
	xs_jan = np.arange(ys_jan.shape[0])*10/365
	


	plt.plot(xs_ias,ys_ias,label="IAS15",color="green")
	plt.plot(xs_whf,ys_whf,label="WHFast",color="red")
	plt.plot(xs_jan,ys_jan,label="JANUS (8th order)",color="blue")
	plt.legend(loc="upper left")
	plt.xlabel("Years before/after velocity flip",size=20)
	plt.ylabel("Difference in semi major axis (AU)",size=20)
	plt.show()


def JANUS_order_comparison_reverse():
	#--- Shows which order of JANUS conserves time symmetry best

	data_4 = np.loadtxt("eccentricity_janus_rev_ord4.txt")
	data_6 = np.loadtxt("eccentricity_janus_rev_ord6.txt")
	data_8 = np.loadtxt("eccentricity_janus_rev_ord8.txt")
	data_10 = np.loadtxt("eccentricity_janus_rev3.txt")
	N = int(data_4.shape[0]/2)
	ys_4 = np.zeros(N)
	ys_6 = np.zeros(N)
	ys_8 = np.zeros(N)
	ys_10 = np.zeros(N)

	for x in range(N):
		ys_4[x] = data_4[N-x]-data_4[N+x]
		ys_6[x] = data_6[N-x]-data_6[N+x]
		ys_8[x] = data_8[N-x]-data_8[N+x]
		ys_10[x] = data_10[N-x]-data_10[N+x]

	xs = np.arange(N)*10/365
	plt.plot(xs,ys_4,label="4th order")
	plt.plot(xs,ys_6,label="6th order")
	plt.plot(xs,ys_8,label="8th order")
	plt.plot(xs,ys_10,label="10th order")
	plt.legend(loc="upper left")
	plt.xlabel("Years before/after velocity flip")
	plt.ylabel("Difference in Eccentricity")
	plt.show()

def cruithne_reversed_eccentricity():
	#--- Plots how time reversbile eccentriticies are for all methods

	data_whf = np.loadtxt("eccentricity_whfast_rev.txt")
	data_ias = np.loadtxt("eccentricity_ias15_rev.txt")
	data_jan = np.loadtxt("eccentricity_janus_rev_ord8.txt")
	ys_whf = np.zeros(int(data_whf.shape[0]/2))
	ys_ias = np.zeros(int(data_ias.shape[0]/2))
	ys_jan = np.zeros(int(data_jan.shape[0]/2))
	N = int(data_whf.shape[0]/2)

	for x in range(ys_whf.shape[0]):
		ys_whf[x]=data_whf[N-x]-data_whf[N+x]
		ys_ias[x]=data_ias[N-x]-data_ias[N+x]
		ys_jan[x]=data_jan[N-x]-data_jan[N+x]

	xs_whf = np.arange(ys_whf.shape[0])*10/365
	xs_ias = np.arange(ys_ias.shape[0])*10/365
	xs_jan = np.arange(ys_jan.shape[0])*10/365
	
	plt.plot(xs_ias,ys_ias,label="IAS15",color="green")
	plt.plot(xs_whf,ys_whf,label="WHFast",color="red")
	plt.plot(xs_jan,ys_jan,label="JANUS (8th order)",color="blue")
	plt.legend(loc="upper left")
	plt.xlabel("Years before/after velocity flip")
	plt.ylabel("Difference in Eccentricity")
	plt.show()

def cruithne_inc_multiplot():
	#--- Plots inclinations of all three methods

	data_whf = np.loadtxt("inclinations_whfast.txt")*180/np.pi
	data_ias = np.loadtxt("inclinations_ias15.txt")*180/np.pi
	data_jan = np.loadtxt("inclinations_janus.txt")*180/np.pi
	xs_whf = np.arange(data_whf.shape[0])*10/365+2019
	xs_ias = np.arange(data_ias.shape[0])*10/365+2019
	xs_jan = np.arange(data_jan.shape[0])*10/365+2019
	plt.plot(xs_whf,data_whf,label="WHFast",color="red")
	plt.plot(xs_ias,data_ias,label="IAS15",color="green")
	plt.plot(xs_jan,data_jan,label="JANUS",color="blue")
	#plt.title("Fluctuation of semi major axis")
	plt.xlabel("Year")
	plt.ylabel("Inclination (degrees)")
	#plt.axvline(7000,color="orange")
	plt.legend()
	plt.show()

def cruithne_argument_peri_multiplot():
	#--- Plots argument of periapsis for all methods

	data_whf = np.loadtxt("argument_peri_whfast.txt")*180/np.pi
	data_ias = np.loadtxt("argument_peri_ias15.txt")*180/np.pi
	data_jan = np.loadtxt("argument_peri_janus.txt")*180/np.pi
	xs_whf = np.arange(data_whf.shape[0])*10/365+2019
	xs_ias = np.arange(data_ias.shape[0])*10/365+2019
	xs_jan = np.arange(data_jan.shape[0])*10/365+2019
	plt.plot(xs_whf,data_whf,label="WHFast",color="red")
	plt.plot(xs_ias,data_ias,label="IAS15",color="green")
	plt.plot(xs_jan,data_jan,label="JANUS",color="blue")
	#plt.title("Fluctuation of semi major axis")
	plt.xlabel("Year")
	plt.ylabel("Argument of pericenter (degrees)")
	#plt.axvline(7000,color="orange")
	plt.legend()
	plt.show()

def cruithne_longitude_ascending_multiplot():
	#--- Plots longitude of ascending node for all methods

	data_whf = np.loadtxt("longitude_ascending_whfast.txt")*180/np.pi
	data_ias = np.loadtxt("longitude_ascending_ias15.txt")*180/np.pi
	data_jan = np.loadtxt("longitude_ascending_janus.txt")*180/np.pi
	xs_whf = np.arange(data_whf.shape[0])*10/365+2019
	xs_ias = np.arange(data_ias.shape[0])*10/365+2019
	xs_jan = np.arange(data_jan.shape[0])*10/365+2019
	plt.plot(xs_whf,data_whf,label="WHFast",color="red")
	plt.plot(xs_ias,data_ias,label="IAS15",color="green")
	plt.plot(xs_jan,data_jan,label="JANUS",color="blue")
	#plt.title("Fluctuation of semi major axis")
	plt.xlabel("Year")
	plt.ylabel("Longitude of Ascending Node (degrees)")
	#plt.axvline(7000,color="orange")
	plt.legend()
	plt.show()

def cruithne_eccentricity_multiplot():
	#--- Plots eccentricity for all methods

	data_whf = np.loadtxt("eccentricity_whfast_rev.txt")
	data_ias = np.loadtxt("eccentricity_ias15_rev.txt")
	data_jan = np.loadtxt("eccentricity_janus_rev_ord8.txt")
	xs_whf = np.arange(data_whf.shape[0])*10/365+2019
	xs_ias = np.arange(data_ias.shape[0])*10/365+2019
	xs_jan = np.arange(data_jan.shape[0])*10/365+2019
	plt.plot(xs_whf,data_whf,label="WHFast",color="red")
	plt.plot(xs_ias,data_ias,label="IAS15",color="green")
	plt.plot(xs_jan,data_jan,label="JANUS (8th order)",color="blue")
	#plt.title("Fluctuation of semi major axis")
	plt.xlabel("Year")
	plt.ylabel("Eccentricity")
	plt.axvline(250000*10/365+2019,color="orange")
	plt.legend()
	plt.show()



def cruithne_mean_anomaly_multiplot():
	#--- Plots mean anomaly for all methods

	data_whf = np.loadtxt("mean_anom_whfast.txt")*180/np.pi
	data_ias = np.loadtxt("mean_anom_ias15.txt")*180/np.pi
	data_jan = np.loadtxt("mean_anom_janus.txt")*180/np.pi
	xs_whf = np.arange(data_whf.shape[0])*10/365+2019
	xs_ias = np.arange(data_ias.shape[0])*10/365+2019
	xs_jan = np.arange(data_jan.shape[0])*10/365+2019
	plt.plot(xs_whf,data_whf,label="WHFast",color="red")
	plt.plot(xs_ias,data_ias,label="IAS15",color="green")
	plt.plot(xs_jan,data_jan,label="JANUS",color="blue")
	#plt.title("Fluctuation of semi major axis")
	plt.xlabel("Year")
	plt.ylabel("Mean Anomaly (degrees)")
	#plt.axvline(7000,color="orange")
	plt.legend()
	plt.show()

def cruithne_semi_major_fft():
	#--- Fft of cruithne semi major axis, not really used

	data = (np.loadtxt("cruithne_semi_major_axis.txt")[:500*365])
	data = 100*(data-np.mean(data))
	plt.plot(data)
	plt.show()
	ft = np.abs(np.fft.fft(data))
	plt.plot(ft)
	plt.show()


def mean_longitudinal_difference():
	#--- Plots difference between mean longitudes of earth and cruithne

	data = np.loadtxt("longitudes_ias15.txt")
	f = lambda x:np.ones(x)/float(x)
	ys = np.convolve(data[:,0]-data[:,1],f(1),"same")*360/(2*np.pi)
	plt.plot(ys)
	plt.show()



def semi_major_mean_long_diff_multi(n,m,s=20):
	#--- Plot of semi major axis against mean longitudinal difference, for all three methods

	assert n<m
	sma_ias = np.loadtxt("cruithne_semi_major_axis_ias15.txt")[n:m]
	sma_whf = np.loadtxt("cruithne_semi_major_axis_whfast.txt")[n:m]
	sma_jan = np.loadtxt("cruithne_semi_major_axis_janus.txt")[n:m]

	p = lambda x: ((x[:,1]-x[:,0]))*360/(2*np.pi)
	mld_ias = p(np.loadtxt("longitudes_ias15.txt")[n:m])
	mld_whf = p(np.loadtxt("longitudes_whfast.txt")[n:m])
	mld_jan = p(np.loadtxt("longitudes_janus.txt")[n:m])
	plt.scatter(mld_ias[::s],sma_ias[::s],label="IAS15",color="green",s=3,alpha=0.5)
	plt.scatter(mld_whf[::s],sma_whf[::s],label="WHFast",color="red",s=3,alpha=0.5)
	plt.scatter(mld_jan[::s],sma_jan[::s],label="JANUS",color="blue",s=3,alpha=0.5)
	plt.xlabel("Mean longitudinal difference (degrees)",size=30)
	plt.ylabel("Semi major axis (AU)",size=30)
	plt.legend()
	plt.show()



def eccentricity_mean_long_diff_multi(n,m,s=10):
	#--- Plot of eccentricity against mean longitudinal difference, for all three methods
	
	assert n<m
	ecc_ias = np.loadtxt("eccentricity_ias15.txt")[n:m]
	ecc_whf = np.loadtxt("eccentricity_whfast.txt")[n:m]
	ecc_jan = np.loadtxt("eccentricity_janus.txt")[n:m]

	p = lambda x: ((x[:,1]-x[:,0]))*360/(2*np.pi)
	mld_ias = p(np.loadtxt("longitudes_ias15.txt")[n:m])
	mld_whf = p(np.loadtxt("longitudes_whfast.txt")[n:m])
	mld_jan = p(np.loadtxt("longitudes_janus.txt")[n:m])
	plt.scatter(mld_ias[::s],ecc_ias[::s],label="IAS15",color="green",s=1,alpha=0.8)
	plt.scatter(mld_whf[::s],ecc_whf[::s],label="WHFast",color="red",s=1,alpha=0.8)
	plt.scatter(mld_jan[::s],ecc_jan[::s],label="JANUS",color="blue",s=1,alpha=0.8)
	plt.xlabel("Mean longitudinal difference (degrees)")
	plt.ylabel("Eccentricity")
	plt.legend()
	plt.show()

def semi_major_mean_long_diff():
	#--- Plot of semi major axis against mean longitudinal difference, for one method

	sma = np.loadtxt("cruithne_semi_major_axis_whfast.txt")[:1000]
	mld = np.loadtxt("longitudes_whfast.txt")[:1000]
	#b,a = signal.butter(3,0.002)
	print(sma.shape)
	mld = ((mld[:,1]-mld[:,0]))*360/(2*np.pi)

	plt.scatter(mld[::10],sma[::10])
	plt.xlabel("Mean longitudinal difference")
	plt.ylabel("Semi major axis")
	plt.show()


def venus_dist():
	#--- Distance between cruithne and venus, to attempt to check for close approach. Not really used much
	vd = np.loadtxt("ven_dist.txt")
	xs_vd = np.arange(vd.shape[0])*10/365+2019
	#b,a = signal.butter(3,0.02)
	#dist_smooth = signal.filtfilt(b,a,vd)
	plt.plot(xs_vd,vd)
	plt.xlabel("Time (years)")
	plt.ylabel("Distance between Venus and cruithne (AU)")
	#Sphere of influence of venus 0.004 au
	plt.hlines(0.004,2019,vd.shape[0]*10/365+2019)
	#h = 0.005
	#plt.scatter(xs_vd[vd<h],vd[vd<h],color="red")
	#plt.plot(xs_vd,dist_smooth)
	plt.show()

def main():
	"""
	file = sys.argv[1]
	data = np.loadtxt(file)
	
	#Split data by astronomical bodydef cruithne_semi_major_multiplot():

	com = data[:,0:3]
	sun = data[:,3:6]
	earth = data[:,6:9]
	moon = data[:,9:12]
	cruithne = data[:,9:12]
	#Shift data so com = origin
	sun+=com
	earth+=com
	cruithne+=com
	"""

	#file = sys.argv[1]
	#data = np.loadtxt(file)
	#plt.scatter(data[:,0],data[:,1])
	#plt.scatter(data[:,2],data[:,3])
	#plt.scatter(data[:,4],data[:,5])
	#plt.show()
	#cruithne_semi_major_multiplot()
	#JANUS_order_comparison_reverse()
	cruithne_reversed_sma()
	#single_reverse_comparison()
	#cruithne_semi_major()
	#cruithne_inc_multiplot()
	#cruithne_longitude_ascending_multiplot()
	#cruithne_argument_peri_multiplot()
	#cruithne_eccentricity_multiplot()
	#cruithne_reversed_eccentricity()
	#cruithne_mean_anomaly_multiplot()
	#mean_longitudinal_difference()
	#semi_major_mean_long_diff_multi(0,501*365)
	#semi_major_mean_long_diff_multi(500*365,1000*365)
	#eccentricity_mean_long_diff_multi(0,500*365)
	#eccentricity_mean_long_diff_multi(500*365,1000*365)

	#venus_dist()
	#semi_major_mean_long_diff()
	#fft(sun,earth,cruithne)
	#earth_cruithne_distance(earth,cruithne)
	#cruithne_freq_corot(earth,cruithne)
	#earth_cruithne_freq(earth,cruithne)
	#cruithne_rotating_frame(sun,earth,cruithne)

main()