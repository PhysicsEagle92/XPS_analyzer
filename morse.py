
import numpy as np
from scipy.constants import h, hbar, c, u
from scipy.special import factorial
from scipy.special import genlaguerre, gamma
import csv as csv

# Factor for conversion from cm-1 to J
FAC = 100 * h * c

class Morse:
    """A class representing the Morse oscillator model of a diatomic."""

    def __init__(self, mA, mB, we, wexe, re, Te):
        """Initialize the Morse model for a diatomic molecule.

        mA, mB are the atom masses (atomic mass units).
        we, wexe are the Morse parameters (cm-1).
        re is the equilibrium bond length (m).
        Te is the electronic energy (minimum of the potential well; origin
            of the vibrational state energies).

        """

        self.mA, self.mB = mA, mB
        self.mu = mA*mB/(mA+mB) * u
        self.we, self.wexe = we, wexe
        self.re = re
        self.Te = Te

        self.De = we**2 / 4 / wexe * FAC
        self.ke = (2 * np.pi * c * 100 * we)**2 * self.mu
        #  Morse parameters, a and lambda.
        self.a = self.calc_a()
        self.lam = np.sqrt(2 * self.mu * self.De) / self.a / hbar
        # Maximum vibrational quantum number.
        self.vmax = int(np.floor(self.lam - 0.5))

        self.make_rgrid()
        self.V = self.Vmorse(self.r)

    def make_rgrid(self, n=10000, rmin=None, rmax=None, retstep=False):
        """Make a suitable grid of internuclear separations."""

        self.rmin, self.rmax = rmin, rmax
        if rmin is None:
            # minimum r where V(r)=De on repulsive edge
            self.rmin = self.re - np.log(2) / self.a
        if rmax is None:
            # maximum r where V(r)=f.De
            f = 0.999
            self.rmax = self.re - np.log(1-f)/self.a
        self.r, self.dr = np.linspace(self.rmin, self.rmax, n,
                                      retstep=True)
        if retstep:
            return self.r, self.dr
        return self.r

    def calc_a(self):
        """Calculate the Morse parameter, a.

        Returns the Morse parameter, a, from the equilibrium
        vibrational wavenumber, we in cm-1, and the dissociation
        energy, De in J.

        """

        return (self.we * np.sqrt(2 * self.mu/self.De) * np.pi *
                c * 100)

    def Vmorse(self, r):
        """Calculate the Morse potential, V(r).

        Returns the Morse potential at r (in m) for parameters De
        (in J), a (in m-1) and re (in m).

        """

        return self.De * (1 - np.exp(-self.a*(r - self.re)))**2

    def Emorse(self, v):
        """Calculate the energy of a Morse oscillator in state v.

        Returns the energy of a Morse oscillator parameterized by
        equilibrium vibrational frequency we and anharmonicity
        constant, wexe (both in cm-1).

        """
        vphalf = v + 0.5
        return (self.we * vphalf - self.wexe * vphalf**2) * FAC

    def calc_turning_pts(self, E):
        """Calculate the classical turning points at energy E.

        Returns rm and rp, the classical turning points of the Morse
        oscillator at energy E (provided in J). rm < rp.

        """

        b = np.sqrt(E / self.De)
        return (self.re - np.log(1+b) / self.a,
                self.re - np.log(1-b) / self.a)

    def calc_psi(self, v, r=None, normed=True, psi_max=1):
        """Calculates the Morse oscillator wavefunction, psi_v.

        Returns the Morse oscillator wavefunction at vibrational
        quantum number v. The returned function is "normalized" to
        give peak value psi_max.

        """

        if r is None:
            r = self.r
        z = 2 * self.lam * np.exp(-self.a*(r - self.re))
        alpha = 2*(self.lam - v) - 1
        psi = (z**(self.lam-v-0.5) * np.exp(-z/2) *
               genlaguerre(v, alpha)(z))
        psi *= psi_max / np.max(psi)
        return psi * 0.37

    def calc_psi_z(self, v, z):
        alpha = 2*(self.lam - v) - 1
        psi = (z**(self.lam-v-0.5) * np.exp(-z/2) *
               genlaguerre(v, alpha)(z))
        Nv = np.sqrt(factorial(v) * (2*self.lam - 2*v - 1) /
                     gamma(2*self.lam - v))
        return Nv * psi

    def plot_V(self, ax, **kwargs):
        """Plot the Morse potential on Axes ax."""

        ax.plot(self.r*1.e10, self.V / FAC + self.Te, **kwargs)

    def get_vmax(self):
        """Return the maximum vibrational quantum number."""

        return int(self.we / 2 / self.wexe - 0.5)

    def draw_Elines(self, vlist, ax, **kwargs):
        """Draw lines on Axes ax representing the energy level(s) in vlist."""

        if isinstance(vlist, int):
            vlist = [vlist]
        for v in vlist:
            E = self.Emorse(v)
            rm, rp = self.calc_turning_pts(E)
            ax.hlines(E / FAC + self.Te, rm*1.e10, rp*1e10, **kwargs)

    def label_levels(self, vlist, ax):
        if isinstance(vlist, int):
            vlist = [vlist]

        for v in vlist:
            E = self.Emorse(v)
            rm, rp = self.calc_turning_pts(E)
            ax.text(s=r'$v={}$'.format(v), x=rp*1e10 + 0.6,
                    y=E / FAC + self.Te, va='center')
    
    def plot_psi(self, vlist, ax, r_plot=None, scaling=1, **kwargs):
        """Plot the Morse wavefunction(s) in vlist on Axes ax."""
        if isinstance(vlist, int):
            vlist = [vlist]
        for v in vlist:
            E = self.Emorse(v)
            if r_plot is None:
                rm, rp = self.calc_turning_pts(E)
                x = self.r[self.r<rp*1.2]
            else:
                x = r_plot
            psi = self.calc_psi(v, r=x, psi_max=self.we/2)
            psi_plot = psi*scaling + self.Emorse(v)/FAC + self.Te
            ax.plot(x*1.e10, psi_plot, **kwargs)
            
    
    def guassian(self, x, mu, sigma, A):
        return A*np.exp(-(x-mu)**2/(2*sigma**2))
    
    def parabola(self, x, a, b, c,p):
        return a*(x-p)**2 + b*(x-p) + c
    
    
    def plot_psi_map(self, vlist, ax,h,r_plot=None, scaling=1, **kwargs):
        """function to plot psi as a heat map around each energy level"""
        #define the x and y values
        w = 3*self.re
        xl = 8000
        yl = 5000
        xm = np.linspace(0, w, xl)
        ym = np.linspace(self.Te, h, yl)
        #define the meshgrid
        mesh = np.zeros((xl,yl))
        #define the energy levels
        for v in vlist:
            E = self.Emorse(v)
            if r_plot is None:
                rm, rp = self.calc_turning_pts(E)
                x = self.r[self.r<rp*1.2]
                print(len(x))
            else:
                x = r_plot
            psi = self.calc_psi(v, r=x, psi_max=self.we/2)
            psi_plot = psi*scaling + self.Emorse(v)/FAC + self.Te
            amplitudes = ((psi_plot - np.min(psi_plot))/(np.max(psi_plot) - np.min(psi_plot)))
            sigma = np.ones(len(psi_plot))*0.5
            for i in range(len(psi_plot)):
                #find disstance from self.Te to self.Emorse in Y
                y_pixel = (h - self.Te)/yl
                y_start = int(np.round((self.Emorse(v)/FAC + self.Te)/y_pixel,0))

                y_end = int(np.round(psi_plot[i]/y_pixel,0))
                y_center = int(np.round((y_start+y_end)/2))            

                max_point = np.where(psi_plot == np.max(psi_plot))[0][0]
                
                amplitudes = (abs(psi_plot) - np.min(abs(psi_plot)))/(np.max(abs(psi_plot)) - np.min(abs(psi_plot)))
                #find distance from 0 to w in X
                x_pixel = xm[-1]/xl
                x_point = int(np.round(x[i]/x_pixel,0))
                #print(x_point,xm[x_point])
                #print(self.Emorse(v)/FAC + self.Te,psi_plot[i],ym[y_center],x_point,y_center)
                if y_start < y_end:
                    for j in range(y_start,y_end):
                        p_value = -self.parabola(ym[j],1,1,-5,ym[y_center])*(psi_plot[i]**2)*0.00000000001 + 1 + 0.03*int(v)
                        #mesh[x_point,j] = p_value
                        if p_value > 0:
                            mesh[x_point,j] = p_value
                        else:
                            mesh[x_point,j] = 0
                else:
                    for j in range(y_end,y_start):
                        p_value = -self.parabola(ym[j],1,1,-5,ym[y_center])*(psi_plot[i]**2)*0.00000000001 + 1 + 0.03*int(v)
                        mesh[x_point,j] = p_value
                        if p_value > 0:
                            mesh[x_point,j] = p_value
                        else:
                            mesh[x_point,j] = 0
        
        mesh = np.rot90(mesh)
        file_name = 'heatmap.csv'
        np.savetxt(file_name, mesh, delimiter=',')
                  
        ax.imshow(mesh, cmap='hot', interpolation='nearest')
        #ax.set_aspect('auto')
        ax.set_xlabel('r (Angstrom)')
        ax.set_ylabel('Energy (cm$^{-1}$)')
        ax.set_title('Morse Potential')
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.tick_params(axis='both', which='minor', labelsize=8)


 

