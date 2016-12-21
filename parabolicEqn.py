# 0 - the constant-average accelaration method (stable)
# 1 - the linear accelaration method (conditionally stable)
# 2 - the central difference method (conditionally stable)
# 3 - the Galerkin method (stable)
# 4 - the backward difference method (stable)
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
    disp = np.zeros(len(t))
    vel = np.zeros(len(t))
    accl = np.zeros(len(t))
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
