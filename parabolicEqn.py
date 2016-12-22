# 0 - the constant-average accelaration method (stable)
# 1 - the linear accelaration method (conditionally stable)
# 2 - the central difference method (conditionally stable) --- This need more explanation
# 3 - the Galerkin method (stable)
# 4 - the backward difference method (stable)
from numpy import zeros
def parabolicEqn(M, C, K, F, t, n=0, x0=0.0, v0=0.0):
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
    invM = 1.0/M
    disp = zeros(len(t))
    vel = zeros(len(t))
    accl = zeros(len(t))
    disp[0] = x0
    vel[0] = v0
    accl[0] = invM * ( F[0] - K * x0 - C * v0);
    for i in range(1,len(t)):
        if i == len(F):
            dF = -F[-1]
        elif i > len(F):
            dF = 0.0
        else:
            dF = F[i] - F[i-1]
        dT = t[i] - t[i-1]
        KK = (2.0/(b * dT**2))*M + (2.0*a/(b*dT))*C + K
        FF = dF + ((2.0/(b * dT))*M + (2.0*a/b)*C)*vel[i-1] + ((1.0/b)*M  + dT*(1.0-a/b)* C)*accl[i-1]
        dU = FF/KK
        disp[i] = disp[i-1] + dU
        vel[i] = vel[i-1] + dT*(1.0-a/b)*accl[i-1] + (2.0*a/(b*dT))*dU  - (2.0*a/b)*vel[i-1]
        accl[i] = accl[i-1] + (2/(b*dT**2))*dU - (2/(b*dT))*vel[i-1] - (1.0/b)*accl[i-1]
    return disp, vel, accl
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    import sys
    import os
    print "\n\n\tThis is a python module to analyse SDoF system subjected to arbitrary forcing"
    print "\tby\t ANIS MOHAMMED VENGASSERI"
    print "\t\t anis.mhd@gmail.com"
    print "\t\t https://github.com/anismhd"
    print "\t The time integration can be in of following types,"
    print "\t\t 0 - The constant-average accelaration method (stable)"
    print "\t\t 1 - The linear accelaration method (conditionally stable)"
    print "\t\t 3 - The Galerkin method (stable)"
    print "\t\t 4 - The backward difference method (stable) \n\n"
    analysis_type = {}
    analysis_type[0] = 'the constant-average accelaration method (stable)'
    analysis_type[1] = 'the linear accelaration method (conditionally stable)'
#    analysis_type[2] = 'the central difference method (conditionally stable) --- This need more explanation'
    analysis_type[3] = 'The Galerkin method (stable)'
    analysis_type[4] = 'The backward difference method (stable)'
# 3 - the Galerkin method (stable)
# 4 - the backward difference method (stable)
    if len(sys.argv) < 2:
        print "This is a demo output of module.."
        print "\nDEMO 1 - Comparison of analytical and module output - Free vibration"
        Wn = 5
        eta = 0.05
        x0 = 1.0
        v0 = 0.0
        print "\t\t{0:40s} = {1:.4f} ras/s".format('Natiral frequency of system ',Wn)
        print "\t\t{0:40s} = {1:.4f} ".format('Damping ratio frequency of system ',eta)
        print "\t\t{0:40s} = {1:.4f} m".format('Initial displacement of system ',x0)
        print "\t\t{0:40s} = {1:.4f} m/s".format('Initial velocity of system ',v0)
        plt.figure()
        t = np.linspace(0,20,1001)
        f = np.zeros(1001)
        Wd = np.sqrt(1-eta**2)*Wn
        disp = np.exp(-eta*Wn*t)*(x0*np.cos(Wd*t)+((v0+eta*Wn*x0)/Wd)*np.sin(Wd*t))
        disp0, vel0, accl0 = parabolicEqn(1.0, 2.0*Wn*eta, Wn**2, f, t, n=0, x0=x0, v0=v0)
        disp1, vel1, accl1 = parabolicEqn(1.0, 2.0*Wn*eta, Wn**2, f, t, n=0, x0=x0, v0=v0)
        disp3, vel3, accl3 = parabolicEqn(1.0, 2.0*Wn*eta, Wn**2, f, t, n=0, x0=x0, v0=v0)
        disp4, vel4, accl4 = parabolicEqn(1.0, 2.0*Wn*eta, Wn**2, f, t, n=0, x0=x0, v0=v0)
        plt.plot(t,disp,'--',label='Analytical')
        plt.plot(t,disp0, c='k', label='constant-average')
        plt.plot(t,disp1, c='c', label='linear')
        plt.plot(t,disp3, c='b', label='Galerkin method')
        plt.plot(t,disp4, c='r', label='backward difference')
        plt.title('DEMO 1 Free vibration')
        plt.ylabel('Displacement (m)')
        plt.xlabel('Time (sec)')
        plt.legend()
        plt.show()
    else:
        print "\t{0:40s} :: {1:40s}".format('Input file name',sys.argv[1])
        if ~os.path.isfile(sys.argv[1]):
            print "\t\tFile {0:s} does not exist....".format(sys.argv[1])
        entering = True
        while entering:
            M = input("\t{0:40s}= ".format('Enter the mass value'))
            C = input("\t{0:40s}= ".format('Enter the damping value'))
            K = input("\t{0:40s}= ".format('Enter the stiffness value'))
            x0 = input("\t{0:40s}= ".format('Enter the initial displacement'))
            v0 = input("\t{0:40s}= ".format('Enter the initial velocity'))
            print "\t0-{0:s}\n\t1-{1:s}\n\t3-{2:s}\n\t4-{3:s}".format(analysis_type[0],analysis_type[1],analysis_type[3],analysis_type[4])
            while True:
                n = input("\t{0:40s}= ".format('Enter the type of analysis from above list'))
                if n in [0,2,3,4]:
                    break
                else:
                    print "\t Invalid analysis type.. please try again.."

            M = float(M)
            C = float(C)
            K = float(K)
            x0 = float(x0)
            v0 = float(v0)
            print "\t\t{0:40s} = {1:s}".format('Input file name',sys.argv[1])
            print "\t\t{0:40s} = {1:.4f}".format('Mass value',M)
            print "\t\t{0:40s} = {1:.4f}".format('Damping value',C)
            print "\t\t{0:40s} = {1:.4f}".format('Stiffness value',K)
            print "\t\t{0:40s} = {1:.4f} ras/s".format('Natiral frequency of system ',np.sqrt(K/M))
            print "\t\t{0:40s} = {1:.4f} ".format('Damping ratio frequency of system ',C/(2.*M*np.sqrt(M/K)))
            print "\t\t{0:40s} = {1:.4f} m".format('Initial displacement of system ',x0)
            print "\t\t{0:40s} = {1:.4f} m/s".format('Initial velocity of system ',v0)
            print "\t\t{0:40s} = {1:s} m/s".format('Type of analysis',analysis_type[n])
            while True:
                key = raw_input("\t\tType yes Proceed to analysis or no to change the values (yes/no) ")
                if (key=='y') or (key=='yes'):
                    print "\t Please re-enter the values again"
                    break
                elif (key=='n') or (key=='no'):
                    entering = False
                    break
                else:
                    print "Invalid input please try again..."
        data = np.loadtxt(sys.argv[1])
        disp, vel, accl = parabolicEqn(M, C, K, data[:,1], data[:,0], n=n, x0=x0, v0=v0)

'''
        kobe = np.loadtxt('Kobe.txt')
        disp0, vel0, accl0 = parabolicEqn(1, 2.0*6.0*0.05, 36.0, kobe[:,1]*9.81, kobe[:,0], n=0, x0=0.0, v0=0.0)
        disp1, vel1, accl1 = parabolicEqn(1, 2.0*6.0*0.05, 36.0, kobe[:,1]*9.81, kobe[:,0], n=1, x0=0.0, v0=0.0)
        disp3, vel3, accl3 = parabolicEqn(1, 2.0*6.0*0.05, 36.0, kobe[:,1]*9.81, kobe[:,0], n=3, x0=0.0, v0=0.0)
        disp4, vel4, accl4 = parabolicEqn(1, 2.0*6.0*0.05, 36.0, kobe[:,1]*9.81, kobe[:,0], n=4, x0=0.0, v0=0.0)
        plt.subplot(2, 2, 1)
        plt.plot( kobe[:,0], kobe[:,1], '--')
        plt.title('Input Accelarogram')
        plt.ylabel('Accelaration (g)')
        plt.xlabel('Time (sec)')

        plt.subplot(2, 2, 2)
        plt.plot( kobe[:,0], disp0, c='k', label='constant-average')
        plt.plot( kobe[:,0], disp1, c='c', label='linear')
        plt.plot( kobe[:,0], disp3, c='b', label='Galerkin method')
        plt.plot( kobe[:,0], disp4, c='r', label='backward difference')
        plt.title('Output Displacement')
        plt.ylabel('Displacement (m)')
        plt.xlabel('Time (sec)')
        plt.legend()

        plt.subplot(2, 2, 3)
        plt.plot( kobe[:,0], vel0, c='k', label='constant-average')
        plt.plot( kobe[:,0], vel1, c='c', label='linear')
        plt.plot( kobe[:,0], vel3, c='b', label='Galerkin method')
        plt.plot( kobe[:,0], vel4, c='r', label='backward difference')
        plt.title('Output Velocity')
        plt.ylabel('Velocity (m/s)')
        plt.xlabel('Time (sec)')
        plt.legend()

        plt.subplot(2, 2, 4)
        plt.plot( kobe[:,0], accl0/9.81, c='k', label='constant-average')
        plt.plot( kobe[:,0], accl1/9.81, c='c', label='linear')
        plt.plot( kobe[:,0], accl3/9.81, c='b', label='Galerkin method')
        plt.plot( kobe[:,0], accl4/9.81, c='r', label='backward difference')
        plt.title('Output Accelarogram')
        plt.ylabel('Accelaration (g)')
        plt.xlabel('Time (sec)')
        plt.legend()

        plt.tight_layout()
        plt.show()
'''