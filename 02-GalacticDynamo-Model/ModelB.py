#working on model B
#m=0seed,S0=1andS1=0
#constant outflow is included
import numpy as np
import math
import matplotlib.pyplot as plt

#functions
def first_derivative(f):
    dfdr = np.zeros_like(f)
    for i in range(3, len(f)-3):
        dfdr[i] = ( -f[i-3] + 9*f[i-2] -45*f[i-1] +45*f[i+1] -9*f[i+2] +f[i+3] ) / (60*dr)

    return dfdr

def second_derivative(f):
    d2fdr2 = np.zeros_like(f)
    for i in range(3, len(f)-3):
        d2fdr2[i] = ( 2*f[i-3] -27*f[i-2] +270*f[i-1] -490*f[i] +270*f[i+1] -27*f[i+2] +2*f[i+3] ) / (180*dr**2)

    return d2fdr2

# periodic: index wraps with % nphi and no ghost zones are needed
def first_derivative_phi(f):
    dfdp = np.zeros_like(f)
    for j in range(nphi):
        dfdp[:, j] = ( -f[:, (j-3)%nphi] + 9*f[:, (j-2)%nphi] -45*f[:, (j-1)%nphi] +45*f[:, (j+1)%nphi] -9*f[:, (j+2)%nphi] +f[:, (j+3)%nphi] ) / (60*dphi)

    return dfdp

def second_derivative_phi(f):
    d2fdp2 = np.zeros_like(f)
    for j in range(nphi):
        d2fdp2[:, j] = ( 2*f[:, (j-3)%nphi] -27*f[:, (j-2)%nphi] +270*f[:, (j-1)%nphi] -490*f[:, j] +270*f[:, (j+1)%nphi] -27*f[:, (j+2)%nphi] +2*f[:, (j+3)%nphi] ) / (180*dphi**2)

    return d2fdp2

a=0

#constants

R = 20 #kpc radial disc boundary
h0 = 0.5 #kpc halfthickness at r=R/2
t0 = 0.73 #Gyr vertical diffusion time
ctau = 1 #tau/tauc
l = 1/5 #kpc turbulence scale length                                                   #h0/5
u = 15 #kmps rms velocity                                                        #15*h0/t0
Tfinal = 10/t0 #corresponds roughly to 10 Gyr when multiplied with t0
deltat = 0.0002 # can use upto 0.006 if its too much
rD = 0.5 #kpc disc flaring radius                                                     #R/2
hD = 1/math.sqrt(2) #kpc disc half thickness at r=0 # some dimensionless value probably  #h0/sqrt2
w0 = 96 #units                                                                          #96/t0
rw = 0.1 #kpc Brandt radius                                                          #R/10
etat = 1
Rm = 1.0e6
NG = 3                                                                   #h0**2/t0
RU = 0.45 #coefficient for vertical velocity dimensionless Uz = RU in this model

tst = []

#creating grid

r = np.linspace(0, 1, 100)

dr = r[1] - r[0]
r_kpc = R*r
#ghost
rghost = np.zeros(106)
rghost[3:103] = r
for i in range(3):
    rghost[2-i] = r[0] - (i+1)*dr
    rghost[103+i] = r[-1] + (i+1)*dr

#NEW (phi): phi grid, 60 points as in the paper (100 x 60). endpoint=False because phi=0 and 2pi are the same point
nphi = 60
phi = np.linspace(0, 2*np.pi, nphi, endpoint=False)
dphi = phi[1] - phi[0]

#CHANGED: r is made a column (100,1) so that every profile below (h, w, G, alpha_k, D, Beq) is also (100,1)
#and broadcasts against the (100,60) fields. r_kpc and rghost above were made before this so they are untouched
r = r[:, None]

#not the constants
h = hD * np.sqrt((1 + (r/rD)**2))
w = w0/np.sqrt((1 + (r/rw)**2))
alpha_k = l**2 * w / h
G = -w * ((r/rw)**2/(1 + (r/rw)**2) )
lamda = h0/R
K = 2 * ((1/l)**2)                                                              #2*(h0/l)**2
D = alpha_k*G*h**3/etat**2
Beq = np.exp(-r) #using Beq as 1 since is arbitrary

j = np.argmin(np.abs(r - 0.5))

# seed amplitudes, m=0 seed means S0=1, S1=0
S0 = 1.0
S1 = 0.0

# the actual formula Bphi = (r / R) * (1 - r / R)**2 * np.exp(-r / R) * (S0 + S1 * np.cos(phi)) last bracket is just 1 cos of S1,S0 values.

Bphi = 1e-10 * r * (1.0-r)**2 * np.exp(-r) * (S0 + S1*np.cos(phi))   #shape (100,60)
Br = -Bphi
alpham = np.zeros((100, nphi))                              

#ghost zone check indexing
Bpghost = np.zeros((106, nphi))                                       #(ghost cells are only in r, phi is periodic)
Brghost = np.zeros((106, nphi))
alphamghost = np.zeros((106, nphi))
Bpghost[3:103] = Bphi
Brghost[3:103] = Br
alphamghost[3:103] = alpham
for i in range(3):
    Brghost[2-i] = -Br[i+1]
    Bpghost[2-i] = -Bphi[i+1]
    alphamghost[2-i] = -alpham[i+1]

    Brghost[103+i] = -Br[98-i]
    Bpghost[103+i] = -Bphi[98-i]
    alphamghost[103+i] = -alpham[98-i]

#timewant
times_gyr = [3.1, 3.5, 4.2, 4.6, 5.5, 7.0, 10.0]
# Convert to dimensionless time
time = [tim/t0 for tim in times_gyr]

allB = []
alpha_m_history = []
dalpha_m_dt_history = []
Q1_history = []
Q2_history = []

#start the time loop
#implement CFL later
t=0
print("HERE!!! before timestep")
while t <= (Tfinal + 1e-9):
  print("HERE!!!  tfinal = ",Tfinal,'t = ', t)


  alpha = alpha_k + alpham

  #using finite difference to get all the spatial derivatives
  dBp_dr = first_derivative(Bpghost)[3:103]
  dBr_dr = first_derivative(Brghost)[3:103]
  d2Bp_dr2 = second_derivative(Bpghost)[3:103]
  d2Br_dr2 = second_derivative(Brghost)[3:103]

  #NEW (phi): phi derivatives (taken on the normal arrays, no ghosts needed)
  dBp_dphi = first_derivative_phi(Bphi)
  dBr_dphi = first_derivative_phi(Br)
  d2Bp_dphi2 = second_derivative_phi(Bphi)
  d2Br_dphi2 = second_derivative_phi(Br)
  dalpham_dphi = first_derivative_phi(alpham)

  #getting the equations
  k1rhsBr = np.zeros_like(Br)
  k1rhsBp = np.zeros_like(Bphi)
  k1rhsalpham = np.zeros_like(alpham)
  i = slice(1, -1)

  k1rhsBr[i] = ((-2 * ctau * alpha[i])/(math.pi * h[i]))*Bphi[i] - ((ctau * math.pi**2)/(4 * h[i]**2))*Br[i] + (ctau * lamda**2)*(-(1/r[i]**2)*Br[i] + (1/r[i])*dBr_dr[i] + d2Br_dr2[i]) - w[i]*dBr_dphi[i] + (ctau * lamda**2)*((1/r[i]**2)*d2Br_dphi2[i] - (2/r[i]**2)*dBp_dphi[i]) - RU*Br[i]/h[i]

  k1rhsBp[i] = G[i] * Br[i] - (2 * ctau / (np.pi * h[i])) * alpha[i] * Br[i] - (ctau * np.pi**2 / (4 * h[i]**2)) * Bphi[i] + (ctau * lamda**2) * (-(1/r[i]**2)*Bphi[i] + (1/r[i])*dBp_dr[i] + d2Bp_dr2[i]) - w[i]*dBp_dphi[i] + (ctau * lamda**2)*((1/r[i]**2)*d2Bp_dphi2[i] + (2/r[i]**2)*dBr_dphi[i]) - RU*Bphi[i]/h[i]

  k1rhsalpham[i] = -K * (ctau * alpha[i] * ((Br[i]**2 + Bphi[i]**2) / Beq[i]**2) + (( 3 * ctau * np.sqrt(-D[i]) / (8 * np.sqrt(np.pi) * h[i] * Beq[i]**2) ) * (Br[i] * Bphi[i]) )) - w[i]*dalpham_dphi[i] - RU*alpham[i]/h[i]



  #beginning RK stage

  #approximate step
  rkBphi = Bphi + deltat * k1rhsBp
  rkBr = Br + deltat * k1rhsBr
  rkalpham = alpham + deltat * k1rhsalpham

  rkalpha = alpha_k + rkalpham

  #creating ghosts for new calculation
  rkBpghost = np.zeros((106, nphi))                                   #CHANGED: 2D
  rkBrghost = np.zeros((106, nphi))
  rkalphamghost = np.zeros((106, nphi))
  rkBpghost[3:103] = rkBphi
  rkBrghost[3:103] = rkBr
  rkalphamghost[3:103] = rkalpham
  for i in range(3):
      rkBrghost[2-i] = -rkBr[i+1]
      rkBpghost[2-i] = -rkBphi[i+1]
      rkalphamghost[2-i] = -rkalpham[i+1]

      rkBrghost[103+i] = -rkBr[98-i]
      rkBpghost[103+i] = -rkBphi[98-i]
      rkalphamghost[103+i] = -rkalpham[98-i]


  #calculating new slopes
  rkdBp_dr = first_derivative(rkBpghost)[3:103]
  rkdBr_dr = first_derivative(rkBrghost)[3:103]
  rkd2Bp_dr2 = second_derivative(rkBpghost)[3:103]
  rkd2Br_dr2 = second_derivative(rkBrghost)[3:103]

  #NEW (phi): phi derivatives for the rk stage
  rkdBp_dphi = first_derivative_phi(rkBphi)
  rkdBr_dphi = first_derivative_phi(rkBr)
  rkd2Bp_dphi2 = second_derivative_phi(rkBphi)
  rkd2Br_dphi2 = second_derivative_phi(rkBr)
  rkdalpham_dphi = first_derivative_phi(rkalpham)

  #new equations
  k2rhsBr = np.zeros_like(Br)
  k2rhsBp = np.zeros_like(Bphi)
  k2rhsalpham = np.zeros_like(alpham)
  i = slice(1, -1)

  k2rhsBr[i] = ((-2 * ctau * rkalpha[i])/(math.pi * h[i]))*rkBphi[i] - ((ctau * math.pi**2)/(4 * h[i]**2))*rkBr[i] + (ctau * lamda**2)*(-(1/r[i]**2)*rkBr[i] + (1/r[i])*rkdBr_dr[i] + rkd2Br_dr2[i]) - w[i]*rkdBr_dphi[i] + (ctau * lamda**2)*((1/r[i]**2)*rkd2Br_dphi2[i] - (2/r[i]**2)*rkdBp_dphi[i]) - RU*rkBr[i]/h[i]

  k2rhsBp[i] = G[i] * rkBr[i] - (2 * ctau / (np.pi * h[i])) * rkalpha[i] * rkBr[i] - (ctau * np.pi**2 / (4 * h[i]**2)) * rkBphi[i] + (ctau * lamda**2) * (-(1/r[i]**2)*rkBphi[i] + (1/r[i])*rkdBp_dr[i] + rkd2Bp_dr2[i]) - w[i]*rkdBp_dphi[i] + (ctau * lamda**2)*((1/r[i]**2)*rkd2Bp_dphi2[i] + (2/r[i]**2)*rkdBr_dphi[i]) - RU*rkBphi[i]/h[i]

  k2rhsalpham[i] = -K * (ctau * rkalpha[i] * ((rkBr[i]**2 + rkBphi[i]**2) / Beq[i]**2)  + (( 3 * ctau * np.sqrt(-D[i]) / (8 * np.sqrt(np.pi) * h[i] * Beq[i]**2) ) * (rkBr[i] * rkBphi[i]))) - w[i]*rkdalpham_dphi[i] - RU*rkalpham[i]/h[i]

  #getting average slope
  rhsBr = (k1rhsBr + k2rhsBr)/2
  rhsBphi = (k1rhsBp + k2rhsBp)/2
  rhsalpham = (k1rhsalpham + k2rhsalpham)/2

  #final change
  Br += deltat * rhsBr
  Bphi += deltat * rhsBphi
  alpham += deltat * rhsalpham
  alpha = alpha_k + alpham

  Bpghost[3:103] = Bphi
  Brghost[3:103] = Br
  alphamghost[3:103] = alpham
  for i in range(3):
      Brghost[2-i] = -Br[i+1]
      Bpghost[2-i] = -Bphi[i+1]
      alphamghost[2-i] = -alpham[i+1]

      Brghost[103+i] = -Br[98-i]
      Bpghost[103+i] = -Bphi[98-i]
      alphamghost[103+i] = -alpham[98-i]

  #reinforce boundaries
  Br[0] = Br[-1] = 0.0
  Bphi[0] = Bphi[-1] = 0.0
  alpham[0] = alpham[-1] = 0.0
  B = np.sqrt(Br**2 + Bphi**2)

#new section
  alpha = alpha_k + alpham
  Q1 = ctau * alpha * ((Br**2 + Bphi**2) / Beq**2)
  Q2 = (
      3 * ctau * np.sqrt(-D)
      / (8 * np.sqrt(np.pi) * h * Beq**2)
  ) * (Br * Bphi)

  dalpha_m_dt = -K * (Q1 + Q2)

  #CHANGED: store only the azimuthal mean (axis=1 is phi). storing the full (100,60) array every step would need ~13 GB
  alpha_m_history.append(alpham.mean(axis=1))
  Q1_history.append(Q1.mean(axis=1))
  Q2_history.append(Q2.mean(axis=1))
  dalpha_m_dt_history.append(dalpha_m_dt.mean(axis=1))


  #snapshot
  for tim in time:
    if np.abs(t-tim) < deltat/2:
      B = np.sqrt(Br**2 + Bphi**2)
      allB.append(B.mean(axis=1))      #CHANGED: azimuthal mean so we get a 1D profile in r to plot
  a = a+1
  t = t+deltat

#NEW (phi): check, how much Bphi still depends on phi at the end. for the m=0 seed this should be ~0 (round off)

allB = np.array(allB)

diagnostic_time = np.arange(len(Q1_history)) * deltat * t0

#now plot
print("HERE!!! endo of loop")
plt.figure(figsize=(7,6))

styles = [
    ('k', '-', 2.5),
    ('#6a3d9a', ':', 2.2),
    ('#4f83cc', '--', 2.5),
    ('#66bb44', '-.', 2.5),
    ('gold', '-', 2.5),
    ('darkorange', '-.', 2.5),
    ('red', '--', 2.5),
]

for i in range(len(times_gyr)):
    color, ls, lw = styles[i]
    plt.plot(
        r_kpc,
        allB[i],           # divide by B0 here if B0 != 1
        color=color,
        linestyle=ls,
        linewidth=lw,
        label=f"$t={times_gyr[i]:.1f}$ Gyr"
    )

# 0.4 Beq/B0 curve
plt.plot(
    r_kpc,
    0.4*Beq[:, 0],         #CHANGED: Beq is (100,1) now so take the column
    'k--',
    linewidth=1.5,
    dashes=(5,5),
    label=r"$0.4B_{\rm eq}/B_0$"
)

plt.xlim(0,12)
plt.ylim(0,0.6)

plt.xlabel(r"$r\ ({\rm kpc})$", fontsize=15)
plt.ylabel(r"$B/B_0$", fontsize=15)

plt.legend(frameon=False, fontsize=12, loc="upper right")

plt.tick_params(direction="in", which="both", top=True, right=True)

plt.tight_layout()
plt.show()