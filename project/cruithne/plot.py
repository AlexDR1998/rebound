import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
#plt.rcParams['axes.facecolor'] = 'black'
import sys





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
	earth_fft = np.fft.fftn(earth)
	cruithne_fft = np.fft.fftn(cruithne)
	plt.plot((np.abs(earth_fft[:,0])))
	plt.plot((np.abs(earth_fft[:,1])))
	#plt.plot(cruithne_fft[:,0].imag)
	

	plt.show()
	print(cruithne_fft.shape)

def earth_cruithne_distance(earth,cruithne):
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
	data = np.loadtxt("cruithne_semi_major_axis_whfast.txt")
	xs = np.arange(data.shape[0])*10/365+2019

	b,a = signal.butter(3,0.002)
	data_smooth = signal.filtfilt(b,a,data)
	peaks,_=signal.find_peaks(data_smooth,height=1.0024)
	print(np.diff(peaks)*10/365)
	plt.plot(xs,data)
	#plt.plot(xs[peaks],data[peaks],"x",color="orange")
	plt.title("Fluctuation of semi major axis")
	plt.xlabel("Year")
	plt.ylabel("Semi major axis (AU)")
	plt.axvline(7000,color="red")
	plt.show()



def cruithne_semi_major_fft():
	data = (np.loadtxt("cruithne_semi_major_axis.txt")[:500*365])
	data = 100*(data-np.mean(data))
	plt.plot(data)
	plt.show()
	ft = np.abs(np.fft.fft(data))
	plt.plot(ft)
	plt.show()


def mean_longitudinal_difference():
	data = np.loadtxt("longitudes.txt")
	f = lambda x:np.ones(x)/float(x)
	ys = np.convolve(data[:,0]-data[:,1],f(1),"same")*360/(2*np.pi)
	plt.plot(ys)
	plt.show()

def semi_major_mean_long_diff():
	sma = np.loadtxt("cruithne_semi_major_axis.txt")[:100000]
	mld = np.loadtxt("longitudes.txt")[:100000]
	#b,a = signal.butter(3,0.002)
	print(sma.shape)
	mld = ((mld[:,1]-mld[:,0]))*360/(2*np.pi)

	plt.scatter(mld[::10],sma[::10])
	plt.xlabel("Mean longitudinal difference")
	plt.ylabel("Semi major axis")
	plt.show()

def main():

	file = sys.argv[1]
	data = np.loadtxt(file)
	
	#Split data by astronomical body
	com = data[:,0:3]
	sun = data[:,3:6]
	earth = data[:,6:9]
	moon = data[:,9:12]
	cruithne = data[:,9:12]
	#Shift data so com = origin
	sun+=com
	earth+=com
	cruithne+=com
	

	#file = sys.argv[1]
	#data = np.loadtxt(file)
	#plt.scatter(data[:,0],data[:,1])
	#plt.scatter(data[:,2],data[:,3])
	#plt.scatter(data[:,4],data[:,5])
	#plt.show()
	cruithne_semi_major()
	#mean_longitudinal_difference()
	#semi_major_mean_long_diff()
	#fft(sun,earth,cruithne)
	#earth_cruithne_distance(earth,cruithne)
	#cruithne_freq_corot(earth,cruithne)
	#earth_cruithne_freq(earth,cruithne)
	#cruithne_rotating_frame(sun,earth,cruithne)

main()