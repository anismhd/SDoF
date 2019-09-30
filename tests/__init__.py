import numpy as np
import matplotlib.pyplot as plt
from SDoF.sdof import sdof
def test1():
	"""
	Testing the module for cosine loading with frequency of range of 0.1 to 100 Hz.
	The single SDoF of natural frequency of 2 Hz in damping ratio of 0.05. 
	"""
	damping_ratio = [0.01,0.05,0.1,0.5]
	natural_freq = 2.0
	frequency_sets = 10**np.linspace(-1,2,31)
	mass = 1.0
	damping = 2.0 * (2.0*np.pi*freq) * np.array(damping_ratio)
	stiffness = (2.0*np.pi*freq)**2
	sdof1 = sdof(mass,damping[0],stiffness)
	sdof2 = sdof(mass,damping[1],stiffness)
	sdof3 = sdof(mass,damping[2],stiffness)
	sdof4 = sdof(mass,damping[3],stiffness)
	plt.plot(2.0*np.pi*frequency_sets,sdof1.harmonic_transfer_func(), label='TF for wn=2Hz, damping_ratio=0.01')
	plt.plot(2.0*np.pi*frequency_sets,sdof2.harmonic_transfer_func(), label='TF for wn=2Hz, damping_ratio=0.05')
	plt.plot(2.0*np.pi*frequency_sets,sdof2.harmonic_transfer_func(), label='TF for wn=2Hz, damping_ratio=0.1')
	plt.plot(2.0*np.pi*frequency_sets,sdof2.harmonic_transfer_func(), label='TF for wn=2Hz, damping_ratio=0.5')
def test2():
	pass