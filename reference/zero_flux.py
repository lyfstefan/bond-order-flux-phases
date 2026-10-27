"""Author-source functions: import-safe extraction, not a confirmed final figure version.
Only trailing whitespace was removed. See docs/source_audit.md.
"""
import numpy as np
from numpy import pi,sin,cos,exp,sqrt,tan
from scipy.optimize import shgo
t=1.0
tp=0.2
Nx=100
Ny=100
siteN=Nx*Ny
error=3

def E1(kx,ky,J1,J2,J3,J4,V): #dispersion of the lower band in the translation-broken phase
    A=2*(t**2+J1*J3*V**2)
    B=2*(4*t**2+(J1+J3)*(J2+J4)*V**2)
    C=2*(t**2+J2*J4*V**2)
    D=2*(J3-J1)*t*V
    E=2*(J4-J2)*t*V
    F=2*(J3-J1)*(J4-J2)*V**2
    return -4*tp*cos(kx)*cos(ky)\
              -sqrt( np.abs(\
                  4*t**2+(J1**2+J2**2+J3**2+J4**2)*V**2+A*cos(2*kx)+B*cos(kx)*cos(ky)+C*cos(2*ky) \
                      +2*E*cos(kx)*sin(ky)+D*sin(2*kx)+2*D*cos(ky)*sin(kx)+F*sin(kx)*sin(ky)+E*sin(2*ky)
                        )
                          )

def E2(kx,ky,J1,J2,J3,J4,V): #dispersion of the lower band in the translation-broken phase
    A=2*(t**2+J1*J3*V**2)
    B=2*(4*t**2+(J1+J3)*(J2+J4)*V**2)
    C=2*(t**2+J2*J4*V**2)
    D=2*(J3-J1)*t*V
    E=2*(J4-J2)*t*V
    F=2*(J3-J1)*(J4-J2)*V**2
    return -4*tp*cos(kx)*cos(ky)\
              +sqrt( np.abs(\
                  4*t**2+(J1**2+J2**2+J3**2+J4**2)*V**2+A*cos(2*kx)+B*cos(kx)*cos(ky)+C*cos(2*ky) \
                      +2*E*cos(kx)*sin(ky)+D*sin(2*kx)+2*D*cos(ky)*sin(kx)+F*sin(kx)*sin(ky)+E*sin(2*ky)
                        )
                          )

def Nparticle12(J1,J2,J3,J4,V,mu):   # the number of particles in the translation-broken phase for a given phi, chemical potential mu
        kx=np.arange(0,2*pi,2*pi/Nx)-pi
        ky=np.arange(0,2*pi,2*pi/Ny)-pi
        Kx,Ky=np.meshgrid(kx,ky)
        E1_matrix=E1(Kx,Ky,J1,J2,J3,J4,V)
        E2_matrix=E2(Kx,Ky,J1,J2,J3,J4,V)

        N1_matrix=(E1_matrix<=mu).astype(int)
        N2_matrix=(E2_matrix<=mu).astype(int)

        N1=N1_matrix.sum()
        N2=N2_matrix.sum()

        Nparticle12=(N1+N2)/2

        return Nparticle12

def Mu12(J1,J2,J3,J4,N,V):    # Chemical potential
    start=-1000
    end=1000
    n=0
    while Nparticle12(J1,J2,J3,J4,V,start)<=Nparticle12(J1,J2,J3,J4,V,end):
        n=n+1
        if n<100:
            mid=(end+start)/2
            if Nparticle12(J1,J2,J3,J4,V,mid)>N+error/2:
                end=mid
            elif Nparticle12(J1,J2,J3,J4,V,mid)<N-error/2:
                start=mid
            else:
                return mid
        else:
            return mid

def Etot(x,N,V):  # total energy
    J1=x[0]
    J2=x[1]
    J3=x[2]
    J4=x[3]
    mu=Mu12(J1,J2,J3,J4,N,V)

    kx=np.arange(0,2*pi,2*pi/Nx)-pi
    ky=np.arange(0,2*pi,2*pi/Ny)-pi
    Kx,Ky=np.meshgrid(kx,ky)
    E1_matrix=E1(Kx,Ky,J1,J2,J3,J4,V)
    E2_matrix=E2(Kx,Ky,J1,J2,J3,J4,V)

    E1_matrix[E1_matrix>mu]=0
    E2_matrix[E2_matrix>mu]=0

    EE1=E1_matrix.sum()
    EE2=E2_matrix.sum()

    E=EE1/2+EE2/2+V*siteN/4*(J1**2+J2**2+J3**2+J4**2)
    return E/siteN

def constraint0(x):
    return x[3]-x[1]+1e-5

def constraint1(x):
    return x[1]-x[2]+1e-5

def constraint2(x):
    return x[2]-x[0]+1e-5

def constraint3(x):
    return x[2]-x[1]+1e-5

def constraint4(x):
    return x[1]-x[0]+1e-5

def constraint5(x):
    return x[3]-x[1]+1e-5

def constraint6(x):
    return x[3]-x[2]+1e-5

def constraint7(x):
    return x[2]-x[1]+1e-5

def constraint8(x):
    return x[1]-x[0]+1e-5

def constraint9(x):
    return -x[0]-x[1]+1e-5

def constraint10(x):
    return -x[0]-x[2]+1e-5

def constraint11(x):
    return -x[0]-x[3]+1e-5

def constraintZeroflux(x,V):
    phi1=np.arctan(V*x[0]/t)
    phi2=np.arctan(V*x[1]/t)
    phi3=np.arctan(V*x[2]/t)
    phi4=np.arctan(V*x[3]/t)
    phi_square=phi4-phi1+phi2-phi3
    return phi_square

def Order(N,V):
    con0 = [{'type': 'ineq', 'fun': constraint0},
           {'type': 'ineq', 'fun': constraint1},
           {'type': 'ineq', 'fun': constraint2},
           {'type': 'ineq', 'fun': constraint9},
            {'type': 'ineq', 'fun': constraint10},
            {'type': 'ineq', 'fun': constraint11},
            {'type': 'eq', 'fun': constraintZeroflux, 'args': (V,)}]

    con1 = [{'type': 'ineq', 'fun': constraint3},
           {'type': 'ineq', 'fun': constraint4},
           {'type': 'ineq', 'fun': constraint5},
           {'type': 'ineq', 'fun': constraint9},
            {'type': 'ineq', 'fun': constraint10},
            {'type': 'ineq', 'fun': constraint11},
            {'type': 'eq', 'fun': constraintZeroflux, 'args': (V,)}]

    con2 = [{'type': 'ineq', 'fun': constraint6},
           {'type': 'ineq', 'fun': constraint7},
           {'type': 'eq', 'fun': constraint8},
           {'type': 'ineq', 'fun': constraint9},
            {'type': 'ineq', 'fun': constraint10},
            {'type': 'ineq', 'fun': constraint11},
           {'type': 'eq', 'fun': constraintZeroflux, 'args': (V,)}]

#    con0 = [{'type': 'ineq', 'fun': constraint0},
#           {'type': 'ineq', 'fun': constraint1},
#           {'type': 'ineq', 'fun': constraint2}]

#    con1 = [{'type': 'ineq', 'fun': constraint3},
#           {'type': 'ineq', 'fun': constraint4},
#           {'type': 'ineq', 'fun': constraint5}]

#    con2 = [{'type': 'ineq', 'fun': constraint6},
#           {'type': 'ineq', 'fun': constraint7},
#           {'type': 'eq', 'fun': constraint8}]

    # Define the i nitial guess for the optimization algorithm
#   bounds = [(-2, 2), (-2, 2),(-2, 2), (0, 2)]
    bounds = [(-1, 0), (-1, 1),(-1, 1), (-1, 1)]
    result  = shgo(Etot,args=(N,V),bounds=bounds,constraints=con0)
    result1 = shgo(Etot,args=(N,V),bounds=bounds,constraints=con1)
    result2 = shgo(Etot,args=(N,V),bounds=bounds,constraints=con2)
    result = result if result1.fun>result.fun else result1
    result = result if result2.fun>result.fun else result2
    return result.x[0],result.x[1],result.x[2],result.x[3],result.fun
