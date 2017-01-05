'''
Anis Mohammed vengasseri
anis.mhd@gmail.com
https://github.com/anismhd

This a module to evaluate response of SDoF system in frequency domain
'''
from numpy import zeros
from numpy.fft import fft
def freqSDoF(M, C, K, F, dT):
	force_DFT = fft(F)
	freq = 2.0*pi*fftfreq(len(F),dT)
	Hw = 1/(-M*freq**2+C*freq*1j+K)
	response_DFT = Hw*force_DFT
	return real(ifft(response_DFT)),real(ifft((0+freq*1j)*response_DFT)),real(ifft((0+freq*1j)**2*response_DFT))
if __name__ == "__main__":
    def modifide_input(string,error_msg,dtype='str'):
        if not(dtype in ['int','float','float64','str']):
            print "Unknown data type... returning Nothing"
            return None
        while True:
            inp = raw_input(string)
            if dtype == 'int':
                try:
                    return int(inp)
                except ValueError:
                    print error_msg
            if dtype == 'float':
                try:
                    return float(inp)
                except ValueError:
                    print error_msg
            if dtype == 'float64':
                try:
                    return float64(inp)
                except ValueError:
                    print error_msg
            if dtype == 'str':
                try:
                    return str(inp)
                except ValueError:
                    print error_msg

    import matplotlib.pyplot as plt
    import numpy as np
    import sys
    import os
    print "\n\n\tThis is a python module to analyse SDoF system subjected to arbitrary forcing in frequency domain"
    print "\tby\t ANIS MOHAMMED VENGASSERI"
    print "\t\t anis.mhd@gmail.com"
    print "\t\t https://github.com/anismhd"
    print "\t The SDoF system system are solved in frequency domain..."