import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import sys


def megno_plot():
	#---Plots megno from .txt file
	
	data = np.loadtxt("megno_ias15.txt")
	plt.plot(np.arange(data.shape[0])*100/365,data)
	#plt.plot(np.arange(data.shape[0])*100/365,np.gradient(data,100/365.0))
	plt.show()

def megno_multi():
	#---Plots both IAS15 and WHFAST megno results together

	data_whf = np.loadtxt("megno_whfast_longer.txt")[:900000]
	data_ias = np.loadtxt("megno_ias15_longer.txt")[:900000]
	xs = np.arange(data_whf.shape[0])*10/365+2019
	plt.plot(xs,data_whf,label="WHFast",color="red")
	plt.plot(xs,data_ias,label="IAS15",color="green")
	plt.legend()
	plt.xlabel("Time (years)")
	plt.ylabel("MEGNO")
	plt.show()

def lyap_from_megno():
	#--- Plots Lyapunov timescales calculated from megno results

	data_whf = np.loadtxt("megno_whfast_longer.txt")[:900000]
	data_ias = np.loadtxt("megno_ias15_longer.txt")[:900000]
	xs = np.arange(data_whf.shape[0])*10/365+2019
	plt.plot(xs,xs/(2*data_whf),label="WHFast",color="red")
	plt.plot(xs,xs/(2*data_ias),label="IAS15",color="green")
	plt.legend()
	plt.xlabel("Time (years)")
	plt.ylabel("Lyapunov Time (years)")
	plt.show()

def lyap_comparison(s=20000*36):
	#--- Plots lyapunov timescales from different methods (MEGNO/LCE, IAS15/WHFAST) from day s onwards.
	#--- Calculates average lyapunov time from day s onwards

	l_data_whf = np.loadtxt("lyap_whfast_longer.txt")[:900000]
	l_data_ias = np.loadtxt("lyap_ias15_longer.txt")[:900000]
	m_data_whf = np.loadtxt("megno_whfast_longer.txt")[:900000]
	m_data_ias = np.loadtxt("megno_ias15_longer.txt")[:900000]
	xs = np.arange(m_data_whf.shape[0])*10/365+2019




	l_data_whf=l_data_whf[s:]
	l_data_ias=l_data_ias[s:]
	m_data_whf=m_data_whf[s:]
	m_data_ias=m_data_ias[s:]
	xs=xs[s:]

	lt_1=np.mean(1/(l_data_whf*365))
	lt_2=np.mean(1/(l_data_ias*365))
	lt_3=np.mean((xs/(2*m_data_whf)))
	lt_4=np.mean((xs/(2*m_data_ias)))

	d_lt_1=np.std(1/(l_data_whf*365))
	d_lt_2=np.std(1/(l_data_ias*365))
	d_lt_3=np.std((xs/(2*m_data_whf)))
	d_lt_4=np.std((xs/(2*m_data_ias)))

	print("LCE WHFast:   "+str(lt_1)+"   "+str(d_lt_1))
	print("LCE IAS15:    "+str(lt_2)+"   "+str(d_lt_2))
	print("MEGNO WHFast: "+str(lt_3)+"   "+str(d_lt_3))
	print("MEGNO IAS15:  "+str(lt_4)+"   "+str(d_lt_4))
	tot = np.array([lt_1,lt_2,lt_3,lt_4])
	print("Mean:         "+str(np.mean(tot))+"   "+str(np.sqrt(np.std(tot)**2 +d_lt_1**2 +d_lt_2**2 +d_lt_3**2 +d_lt_4**2)))

	plt.plot(xs,1/(l_data_whf*365),label="LCE WHFast",color="red")
	plt.plot(xs,xs/(2*m_data_whf),label="MEGNO WHFast",color="#990044")
	plt.plot(xs,1/(l_data_ias*365),label="LCE IAS15",color="#00FF00")
	plt.plot(xs,xs/(2*m_data_ias),label="MEGNO IAS15",color="#009999")
	plt.legend()
	plt.xlabel("Time (years)",size=20)
	plt.ylabel("Lyapunov Time (years)",size=20)
	plt.show()

def lyap_multi():
	#--- Plots lyapunov times from LCE only, for IAS15 and WHFAST. Calculates averages.
	#--- Obsolete, see lyap_comparison()

	data_whf = np.loadtxt("lyap_whfast_longer.txt")[:900000]
	data_ias = np.loadtxt("lyap_ias15_longer.txt")[:900000]
	xs = np.arange(data_whf.shape[0])*10/365+2019
	print(np.mean(1/(data_whf[100000:]*365)))
	print(np.mean(1/(data_ias[100000:]*365)))
	plt.plot(xs,1/(data_whf*365),label="WHFast",color="red")
	plt.plot(xs,1/(data_ias*365),label="IAS15",color="green")
	plt.xlabel("Time (years)")
	plt.ylabel("Lyapunov Time (years)")
	plt.legend()
	plt.show()


def lyap_plot():
	#--- Plots lyapunov time from LCE from .txt file
	data = np.loadtxt("lyap_ias15.txt")
	plt.plot(1/(data[:]*365))
	plt.show()

def main():
	#megno_plot()
	#lyap_plot()
	#megno_multi()
	#lyap_multi()
	lyap_comparison()
main()