import numpy as np
import matplotlib.pyplot as plt
from plot_config import apply_plot_settings  # Import the plot settings

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

def get_C(N, CLa, theta, span, isElliptic):
    #alocate memory
    C = np.zeros([N,N])
    #Calculate C
    C=np.zeros([N,N])
    #calcualte elements of c matrix
    for j in range(N):
        #calculate first and last row
        if isElliptic:
            C[0,j]=4*(j+1)**2
            C[-1,j]=((-1)**j)*4*(j+1)**2
        else:
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
    Aalph = np.zeros([N])
    Aomeg = np.zeros([N])
    Adelt = np.zeros([N])
    Apbar = np.zeros([N])
    Atotl = np.zeros([N])
    for j in range(N):
        Aalph[j] = aj[j, 0] * np.deg2rad(aoa)
        Aomeg[j] = -1*bj[j, 0] * np.deg2rad(Omega)
        Adelt[j] = cj[j, 0] * np.deg2rad(das)
        Apbar[j] = dj[j, 0] * pbar

        Atotl[j] = aj[j, 0] * np.deg2rad(aoa) - bj[j, 0] * np.deg2rad(Omega) + cj[j, 0] * np.deg2rad(das) + dj[j, 0] * pbar

    return Atotl, Aalph, Aomeg, Adelt, Apbar


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

def get_section_lift(A, Ra, theta, n):
    # Calculate the lift distribution for each section
    CL = np.zeros(n)
    
    for i in range(n):
        CL[i] = 4*Ra*sum(A[j]*np.sin((j+1)*theta[i]) for j in range(n))

    return CL

def generate_plot(chord, theta, sar, sat, har, hat, RA, ax):
    leading_edge = (chord * .25)/RA
    trailing_edge = (chord * -.75)/RA
    z = theta2s(theta)/2

    ax.plot([z, z], [leading_edge, trailing_edge], 'b:', linewidth=0.5)
    ax.plot(z, leading_edge, 'k-', linewidth=0.5)
    ax.plot(z, trailing_edge, 'k-', linewidth=0.5)

    ail_root = sar/2
    ail_tip = sat/2
    ail_x = [ail_root, ail_tip]
    ail_y = [-har/RA, -hat/RA]
    ax.plot(ail_x, ail_y, 'r:', linewidth=0.5)
    ax.plot([-x for x in ail_x], ail_y, 'r:', linewidth=0.5)  # Reflect across y-axis

    # Plot lines from aileron coordinates to trailing edge
    for x, y in zip(ail_x, ail_y):
        ax.plot([x, x], [y, trailing_edge[np.argmin(abs(z - x))]], 'r:', linewidth=0.5)
    for x, y in zip([-x for x in ail_x], ail_y):
        ax.plot([x, x], [y, trailing_edge[np.argmin(abs(z - x))]], 'r:', linewidth=0.5)

    # Set axis limits to match the reference style
    x_max = max(z) + 0.02
    x_min = min(z) - 0.02
    y_max = max(leading_edge) + 0.02
    y_min = min(trailing_edge) - 0.02
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # Adding labels
    ax.set_xlabel(r'$z/b$')
    ax.set_ylabel(r'$x/b$')

    ax.invert_xaxis()  # Invert x-axis to match the reference style

    # Apply plot settings
    apply_plot_settings(ax, aspect_equal=True)



def generate_Atot_plot(CL, theta, y_label, ax):
    z = theta2s(theta)/2
    # Set axis limits
    x_max = max(z) + 0.02
    x_min = min(z) - 0.02
    y_max = max(CL) + 0.01
    y_min = min(CL) - 0.01
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    # Determine tick interval based on CL range
    range_CL = y_max
    if range_CL > 1:
        My = 2
        miny = 0.5
    elif range_CL > .2:
        My = 0.2
        miny = 0.05
    elif range_CL > 0.1:
        My = 0.1
        miny = 0.02
    else:
        My = 0.02
        miny = 0.005

    apply_plot_settings(ax, False, 0.2, My, 0.05, miny)
    
    
    ax.plot(z, CL, 'b-', linewidth=0.5)
    ax.set_xlabel('z/b')  # Set x-axis label
    ax.set_ylabel(y_label)  # Set y-axis label

    

    ax.invert_xaxis()  # Invert x-axis to match the reference style

    



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


C = get_C(N, CLa, theta, bw, isElliptic)

coef = get_coef(C, N, theta, de, he)
print("aj: ", coef['aj'])
print("bj: ", coef['bj'])
print("cj: ", coef['cj'])

Atotl, Aalph, Aomeg, Adelt, Apbar = get_A(coef, aoa_deg, Omega_deg, das_deg, pbar)

CL_totl, CDi_totl, Cl_totl, Cn_totl = get_wing_coef(Atotl, Ra, pbar, N)
CL_alph, CDi_alph, Cl_alph, Cn_alph = get_wing_coef(Aalph, Ra, pbar, N)
CL_omeg, CDi_omeg, Cl_omeg, Cn_omeg = get_wing_coef(Aomeg, Ra, pbar, N)
CL_delt, CDi_delt, Cl_delt, Cn_delt = get_wing_coef(Adelt, Ra, pbar, N)
CL_pbar, CDi_pbar, Cl_pbar, Cn_pbar = get_wing_coef(Apbar, Ra, pbar, N)
print("Total Coefficients")
print("CL_totl: ", CL_totl)
print("CDi_totl: ", CDi_totl)
print("Cl_totl: ", Cl_totl)
print("Cn_totl: ", Cn_totl)

print("Alpha Components")
print("CL_alph: ", CL_alph)
print("CDi_alph: ", CDi_alph)
print("Cl_alph: ", Cl_alph)
print("Cn_alph: ", Cn_alph)

print("Omega Components")
print("CL_omeg: ", CL_omeg)
print("CDi_omeg: ", CDi_omeg)
print("Cl_omeg: ", Cl_omeg)
print("Cn_omeg: ", Cn_omeg)

print("Delta Components")
print("CL_delt: ", CL_delt)
print("CDi_delt: ", CDi_delt)
print("Cl_delt: ", Cl_delt)
print("Cn_delt: ", Cn_delt)

print("pbar Components")
print("CL_delt: ", CL_delt)
print("CDi_delt: ", CDi_delt)
print("Cl_delt: ", Cl_delt)
print("Cn_delt: ", Cn_delt)

print("Symmetric Components")
print("CL: ", CL_alph+CL_omeg)
print("CDi: ", CDi_alph+CDi_omeg)
print("Cl: ", Cl_alph+Cl_omeg)
print("Cn: ", Cn_alph+Cn_omeg)

print("Antisymmetric Components")
print("CL: ", CL_delt+CL_pbar)
print("CDi: ", CDi_delt+CDi_pbar)
print("Cl: ", Cl_delt+Cl_pbar)
print("Cn: ", Cn_delt+Cn_pbar)

#calculate symmetric and antisymmetric lift distributions
CL_symt = CL_alph+CL_omeg
CDi_symt = CDi_alph+CDi_omeg
Cl_symt = Cl_alph+Cl_omeg
Cn_symt = Cn_alph+Cn_omeg

CL_asymt = CL_delt+CL_pbar
CDi_asymt = CDi_delt+CDi_pbar
Cl_asymt = Cl_delt+Cl_pbar
Cn_asymt = Cn_delt+Cn_pbar


CL_totl_dist = get_section_lift(Atotl, Ra, theta, N)
CL_alph_dist = get_section_lift(Aalph, Ra, theta, N)
CL_omeg_dist = get_section_lift(Aomeg, Ra, theta, N)
CL_delt_dist = get_section_lift(Adelt, Ra, theta, N)
CL_pbar_dist = get_section_lift(Apbar, Ra, theta, N)

# Define the number of rows and columns explicitly for a 4x2 grid layout
fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(15, 12))

fig.suptitle("Lift Distributions", fontsize=16)

# Define titles for each subplot
titles = [
    "Wing planform",
    "CL={:.14f} CD={:.14f} Cl={:.14f} Cn={:.14f}".format(CL_totl, CDi_totl, Cl_totl, Cn_totl),
    "CL={:.14f} CD={:.14f} Cl={:.14f} Cn={:.14f}".format(CL_alph, CDi_alph, Cl_alph, Cn_alph),
    "CL={:.14f} CD={:.14f} Cl={:.14f} Cn={:.14f}".format(CL_delt, CDi_delt, Cl_delt, Cn_delt),
    "CL={:.14f} CD={:.14f} Cl={:.14f} Cn={:.14f}".format(CL_omeg, CDi_omeg, Cl_omeg, Cn_omeg),
    "CL={:.14f} CD={:.14f} Cl={:.14f} Cn={:.14f}".format(CL_pbar, CDi_pbar, Cl_pbar, Cn_pbar),
    "CL={:.14f} CD={:.14f} Cl={:.14f} Cn={:.14f}".format(CL_symt, CDi_symt, Cl_symt, Cn_symt),
    "CL={:.14f} CD={:.14f} Cl={:.14f} Cn={:.14f}".format(CL_asymt, CDi_asymt, Cl_asymt, Cn_asymt),
]

# Generate each subplot with titles
generate_plot(ctheta, theta, sar, sat, har, hat, Ra, axs[0, 0])
axs[0, 0].text(0.5, 1.05, titles[0], transform=axs[0, 0].transAxes, ha="center", va="bottom")

generate_Atot_plot(CL_totl_dist, theta, "Total", axs[0, 1])
axs[0, 1].text(0.5, 1.05, titles[1], transform=axs[0, 1].transAxes, ha="center", va="bottom")

generate_Atot_plot(CL_alph_dist, theta, r"$\alpha$ component", axs[1, 0])
axs[1, 0].text(0.5, 1.05, titles[2], transform=axs[1, 0].transAxes, ha="center", va="bottom")

generate_Atot_plot(CL_delt_dist, theta, r"$\delta_{as}$ component", axs[1, 1])
axs[1, 1].text(0.5, 1.05, titles[3], transform=axs[1, 1].transAxes, ha="center", va="bottom")

generate_Atot_plot(CL_omeg_dist, theta, r"$\Omega$ component", axs[2, 0])
axs[2, 0].text(0.5, 1.05, titles[4], transform=axs[2, 0].transAxes, ha="center", va="bottom")

generate_Atot_plot(CL_pbar_dist, theta, r"$\bar{p}$ component", axs[2, 1])
axs[2, 1].text(0.5, 1.05, titles[5], transform=axs[2, 1].transAxes, ha="center", va="bottom")

generate_Atot_plot(CL_alph_dist + CL_omeg_dist, theta, "Symmetric component", axs[3, 0])
axs[3, 0].text(0.5, 1.05, titles[6], transform=axs[3, 0].transAxes, ha="center", va="bottom")

generate_Atot_plot(CL_delt_dist + CL_pbar_dist, theta, "Antisymmetric component", axs[3, 1])
axs[3, 1].text(0.5, 1.05, titles[7], transform=axs[3, 1].transAxes, ha="center", va="bottom")

# Adjust layout for readability
plt.tight_layout()
plt.show()
