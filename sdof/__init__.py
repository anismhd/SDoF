from numpy import zeros, sqrt, fft, log2, ceil

class sdof(object):
	"""docstring for sdof"""
	def __init__(self, mass, damping, stiffness):
		self.mass = mass
		self.stiffness = stiffness
		self.damping = damping

	def get_natural_freq(self):
		return sqrt(stiffness/mass)

	def get_damping_ratio(self):
		return (damping/(2*mass*self._get_natural_freq))

	def harmonic_transfer_func(self, freq):
		return 1.0/(-self.mass*freq**2+self.damping*freq*1j+self.stiffness)

	def __call__(self,force,dt,x0=0.0,v0=0.0, method='fft', n=0):
		if method == 'fft':
			# zeros padding of forcing function of fft
			length_force = len(force)
			nearest_fft_pnt = 2**int( ceil(log2(length_force)) )
			force_DFT = fft(F,nearest_fft_pnt)
			freq = 2.0*pi*fftfreq(nearest_fft_pnt,dt)
			response_DFT = self.harmonic_transfer_func(freq)*force_DFT
			return real(ifft(response_DFT)),\
				real(ifft((0+freq*1j)*response_DFT)),\
				real(ifft((0+freq*1j)**2*response_DFT))
		if method == 'time-integration':
			if n == 0:
				a = 0.5
				b = 0.5
			elif n == 1:
				a = 0.5
				b = 1.0/3.0
			elif n == 2:
				a = 0.5
				b = 0.0
			elif n == 3:
				a = 3.0/2.0
				b = 8.0/5.0
			elif n == 4:
				a = 3.0/2.0
				b = 2.0
			invM = 1.0/self.mass
			disp = zeros(len(force))
			vel = zeros(len(force))
			accl = zeros(len(force))
			disp[0] = x0
			vel[0] = v0
			accl[0] = invM * ( force[0] - self.stiffness * x0 - self.damping * v0);
			for i in range(1,len(force)):
				if i == len(force):
					dF = -force[-1]
				elif i > len(F):
					dF = 0.0
				else:
					dF = force[i] - force[i-1]
				KK = (2.0/(b * dt**2))*self.mass + \
					(2.0*a/(b*dt))*self.damping + self.stiffness
				FF = dF + ((2.0/(b * dt))*self.mass + \
					(2.0*a/b)*self.damping)*vel[i-1] + \
					((1.0/b)*self.mass  + dt*(1.0-a/b)* self.damping)*accl[i-1]
				dU = FF/KK
				disp[i] = disp[i-1] + dU
				vel[i] = vel[i-1] + dt*(1.0-a/b)*accl[i-1] + (2.0*a/(b*dt))*dU  - (2.0*a/b)*vel[i-1]
				accl[i] = accl[i-1] + (2/(b*dt**2))*dU - (2/(b*dt))*vel[i-1] - (1.0/b)*accl[i-1]
			return disp, vel, accl