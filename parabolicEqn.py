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
    kobe = np.loadtxt('Kobe.txt')
    disp0, vel0, accl0 = parabolicEqn(1, 2.0*6.0*0.05, 36.0, kobe[:,1]*9.81, kobe[:,0], n=0, x0=0.0, v0=0.0)
    disp1, vel1, accl1 = parabolicEqn(1, 2.0*6.0*0.05, 36.0, kobe[:,1]*9.81, kobe[:,0], n=1, x0=0.0, v0=0.0)
#    disp2, vel2, accl2 = parabolicEqn(1, 2.0*6.0*0.05, 36.0, kobe[:,1]*9.81, kobe[:,0], n=2, x0=0.0, v0=0.0)
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
#    plt.plot( kobe[:,0], disp2, c='g')
    plt.plot( kobe[:,0], disp3, c='b', label='Galerkin method')
    plt.plot( kobe[:,0], disp4, c='r', label='backward difference')
    plt.title('Output Displacement')
    plt.ylabel('Displacement (m)')
    plt.xlabel('Time (sec)')
    plt.legend()

    plt.subplot(2, 2, 3)
    plt.plot( kobe[:,0], vel0, c='k', label='constant-average')
    plt.plot( kobe[:,0], vel1, c='c', label='linear')
#    plt.plot( kobe[:,0], vel2, c='g')
    plt.plot( kobe[:,0], vel3, c='b', label='Galerkin method')
    plt.plot( kobe[:,0], vel4, c='r', label='backward difference')
    plt.title('Output Velocity')
    plt.ylabel('Velocity (m/s)')
    plt.xlabel('Time (sec)')
    plt.legend()

    plt.subplot(2, 2, 4)
    plt.plot( kobe[:,0], accl0/9.81, c='k', label='constant-average')
    plt.plot( kobe[:,0], accl1/9.81, c='c', label='linear')
#    plt.plot( kobe[:,0], accl2, c='g')
    plt.plot( kobe[:,0], accl3/9.81, c='b', label='Galerkin method')
    plt.plot( kobe[:,0], accl4/9.81, c='r', label='backward difference')
    plt.title('Output Accelarogram')
    plt.ylabel('Accelaration (g)')
    plt.xlabel('Time (sec)')
    plt.legend()

    plt.tight_layout()
    plt.show()