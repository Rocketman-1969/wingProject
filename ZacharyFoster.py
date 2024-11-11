import numpy as np

def interpolate(x1,x2,x,y1,y2):
    return (y2-y1)/(x2-x1)*(x-x1)+y1

def theta2s(theta):
    return -1.*np.cos(theta)

def s2theta(s):
    return np.arccos(-s)

def omega(th):
    return abs(np.cos(th))

def chi(th):
    s = theta2s(th)
    if s < -sat or abs(s) <= sar or sat < s:
        return 0.
    elif -sat <= s and s < -sar:
        return flapEfficiency(th)
    else:
        return -flapEfficiency(th)

def flapEfficiency(th):
    s = abs(theta2s(th))
    hs = interpolate(sar,sat,s,har,hat)
    c = chord(th)
    te = 0.75*c
    fcf = (te - hs) / c
    thf = np.arccos(2*fcf-1)
    ei = 1 - (thf - np.sin(thf)) / np.pi
    return de*he*ei

def chord(th):
    # ~ s = theta2s(th)
    # ~ cw = 1.
    # ~ cr, ct = 2.*cw/(1.+Rt), 2.*Rt*cw/(1.+Rt)
    # ~ return cr + abs(s) * (ct - cr)
    return 4/np.pi*np.sin(th)


def get_theta(N):
    theta = np.zeros(N)
    for i in range(N):
        theta[i] = i*np.pi/(N-1)

    return theta

def get_C(N, CLa, theta, span):
    #alocate memory
    C = np.zeros([N,N])
    #Calculate C
    C=np.zeros([N,N])
    #calcualte elements of c matrix
    for j in range(N):
        #calculate first and last row
        C[0,j]=(j+1)**2
        C[-1,j]=((-1)**j)*(j+1)**2
        for i in range(1,N-1):
            #calculate everything else
            C[i,j]=(4*span/(chord(theta[i])*CLa)+((j+1)/np.sin(theta[i])))*np.sin((j+1)*theta[i])
    print(C)
    return C

def get_coef(C, N, theta,de,he):
    #alocate memory
    aj = np.ones([N,1])
    bj = omega(theta).reshape(-1, 1)
    cj = np.zeros([N,1])
    for i in range(len(theta)-1):
        cj[i] = chi(theta[i])
    dj = np.cos(theta).reshape(-1, 1)

    #solve Cinv
    eye = np.eye(N)
    Cinv = np.linalg.solve(C,eye)

    #solve aj
    aj = np.matmul(Cinv,aj)
    #solve bj
    bj = np.matmul(Cinv,bj)
    #solve cj
    cj = np.matmul(Cinv,cj)
    #solve dj
    dj = np.matmul(Cinv,dj)

    coef = {'aj':aj, 'bj':bj, 'cj':cj, 'dj':dj}

    return coef

def get_A(coef, aoa, Omega, das, pbar):
    aj, bj, cj, dj = coef['aj'], coef['bj'], coef['cj'], coef['dj']
    A = np.zeros([N])
    for j in range(N):
        A[j] = aj[j, 0] * np.deg2rad(aoa) - bj[j, 0] * np.deg2rad(Omega) + cj[j, 0] * np.deg2rad(das) + dj[j, 0] * pbar
    return A

def get_wing_coef(A,Ra, pbar, N):
    # Calculate aerodynamic coefficients using equations (3) to (6)
    CL = np.pi * Ra * A[0]  # Equation (3)

    # Equation (4) for induced drag coefficient C_Di
    CDi = np.pi * Ra * (-0.5 * A[1] * pbar + sum(j * A[j-1]**2 for j in range(1, N)))

    # Equation (5) for rolling moment coefficient C_l
    Cl = -np.pi * Ra * A[1] / 4

    # Equation (6) for yawing moment coefficient C_n
    Cn = (np.pi * Ra / 4) * (-0.5 * (A[0] + A[2]) * pbar + sum((2 * j - 1) * A[j - 2] * A[j-1] for j in range(2, N)))

    return CL, CDi, Cl, Cn

isElliptic = True       ## flag that can be used to know if the planform is elliptic (useful when creating the C matrix)

Ra = 8
bw = Ra         ## implies cw = 1, which means chord(th) should give a cw of 1
cw = 1
Rt = 0.4
CLa = 2*np.pi

aoa_deg     = 5
Omega_deg   = -1
das_deg      = 7
pbar        = -0.1

sar     = 0.5
sat     = 0.9
fcf_ar  = 0.25
fcf_at  = 0.35
de      = 1.
he      = 0.85

car = chord(s2theta(sar))
cat = chord(s2theta(sat))
har = (0.75-fcf_ar)*car
hat = (0.75-fcf_at)*cat

N = 99

theta = get_theta(N)
ctheta = chord(theta)


C = get_C(N, CLa, theta, bw)

coef = get_coef(C, N, theta,de,he)
print("aj: ", coef['aj'])
print("bj: ", coef['bj'])
print("cj: ", coef['cj'])


A = get_A(coef, aoa_deg, Omega_deg, das_deg, pbar)
CL, CDi, Cl, Cn = get_wing_coef(A,Ra, pbar, N)
print("CL: ", CL)
