import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import sys


def megno_plot():
	data = np.loadtxt("megno.txt")
	plt.plot(np.arange(data.shape[0])*100/365,data)
	plt.plot(np.arange(data.shape[0])*100/365,np.gradient(data,100/365.0))
	plt.show()

def lyap_plot():
	data = np.loadtxt("lyap.txt")
	plt.plot(1/(data[5000:]*365))
	plt.show()

def main():
	megno_plot()
	lyap_plot()

main()