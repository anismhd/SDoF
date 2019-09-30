from hyperbolicEqn import hyperbolicEqn
from freqSDoF import freqSDoF
import matplotlib.pyplot as plt
import numpy as np
if __name__ == '__main__':
	print "\nThis is test for SDoF solver"
	print "\nTwo types of solvers \n\t Time integration\n\t Frequency domain solution"
	print "\nTesting of module against analytical solutions"
	print "\tMass = 10 kg\n\tStiffnes = 40 N/m\n\tDamping=1 N/m/s"
	print "\tForce = sin(3t)"
	Wn = 5
	eta = 0.05
	M = 1
	K = Wn**2
	C = 2.0*Wn*eta
	t = np.linspace(0,20,10001)
	f = np.sin(2.5*t)
	disp, velc, accl = hyperbolicEqn(M, C, K, f, t, n=0)
	plt.plot(t,disp)
	disp, velc, accl = freqSDoF(M, C, K, f, t[1]-t[0])
	plt.plot(t,disp,label='Frequency',linestyle='--')
	plt.legend()
	plt.show()