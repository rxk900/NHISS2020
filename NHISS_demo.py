import numpy as np

AMU = {"H": 1.00794, "He": 4.002602, "Li": 6.941, "Be": 9.012182, "B": 10.811, "C": 12.0107, "N": 14.0067, "O": 15.9994, "F": 18.9984032, "Ne": 20.1797}

def CoM(M,X,Y,Z):
    return np.sum(M*X)/np.sum(M),np.sum(M*Y)/np.sum(M), np.sum(M*Z)/np.sum(M)

def Inertia(M,X,Y,Z):
    Ix=np.sum(M*(Y*Y+Z*Z))
    Iy=np.sum(M*(X*X+Z*Z))
    Iz=np.sum(M*(X*X+Y*Y))
    return Ix,Iy,Iz

def RotCon(Ix,Iy,Iz):
    import numpy as np
    # 1 AMU= 1./6.022e+23 kg^-1
    # I in kg m^2 is I in au (1/6.022e+23)*1e-20
    Ix_SI = (Ix/6.022e+23)*1e-20
    Iy_SI = (Iy/6.022e+23)*1e-20
    Iz_SI = (Iz/6.022e+23)*1e-20
    print(Ix_SI,Iy_SI,Iz_SI)
    h=6.62607004e-34 # m^2 kg s^-1
    pi=np.pi
    Ix_Hz=h/(8*pi*pi*Ix_SI)
    Iy_Hz=h/(8*pi*pi*Iy_SI)
    Iz_Hz=h/(8*pi*pi*Iz_SI)
    return Ix_Hz/1e3,Iy_Hz/1e3,Iz_Hz/1e3

atom=np.zeros(10,dtype=str)
mass=np.zeros(10)
x=np.zeros(10)
y=np.zeros(10)
z=np.zeros(10)

# Read in structures for the 6 candidate molecules and print moments of inertia
for i in range (6):
    filename="FA-DFM-"+str(i+1)+".xyz"
    print("Structure from ",filename)
    molecule=np.genfromtxt(filename,skip_header=2,dtype=None)
    count=0
    for line in molecule:
        atom[count],x[count],y[count],z[count]=molecule[count]
        mass[count]=AMU[atom[count]]
        count += 1
    Ix,Iy,Iz=Inertia(mass,x,y,z)
    print("Moments of inertia is %.2f  %.2f  %.2f in amu A^-2" % (Ix,Iy,Iz))
    print("Moments of inertia is %.2f  %.2f  %.2f in MHz" % RotCon(Ix,Iy,Iz))
