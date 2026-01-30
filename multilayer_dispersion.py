# class for defining magnetic film geometrical parameters  
import numpy as np
from numpy import pi
μ0 = 4*pi *1e-7  # Permeability of free space (H/m or T·m/A)

class MagneticFilm:
    def __init__(self, number_of_layers, total_thickness):
        """
        Initialize a MagneticFilm object.

        Args:
            number_of_layers (int): Number of layers in the film.
            thickness (float): Total thickness of the film (nm, μm, etc.).
        """
        self.number_of_layers = number_of_layers
        self.layer_thickness = total_thickness / number_of_layers  # thickness of each layer  

        self.magnetic_parameters = self.Magnetic_parameters(self)



    class Magnetic_parameters:
        def __init__(self,outer):
            self.outer = outer
            N = self.outer.number_of_layers
            d = self.outer.layer_thickness
            self.define_funct("th_u", lambda eta: 0) # anisotropy angle theta
            self.define_funct("phi_u", lambda eta: 0) # anisotropy angle phi
            self.define_funct("Ku", lambda eta: 0) # anisotropy constant

            self.define_funct("Ms", lambda eta: 800e3) # saturation magnetization (A/m)
            self.define_funct("phi", lambda eta: 0) # magnetization angle phi
            self.define_funct("th", lambda eta: 0) # magnetization angle theta

            self.define_funct("Aex", lambda eta: 11e-12) # exchange constant (J/m)

            self.define_funct("Db", lambda eta: 0) # DMI

            self.define_funct("gamma", lambda eta: 2*pi*29.58e9) # gyromagnetic ratio (m/As)

            self.define_funct("xi", lambda eta: d*(eta+1/2)) # position of the center of each layer

            self.define_interlayer_interaction("J", (2/d)*self.Aex_array) # interlayer exchange coupling (J/m^2)
            self.define_interlayer_interaction("Db_inter", (1/d)*self.Db_array) # interlayer DMI coupling (J/m^2)
            self.define_interlayer_interaction("J2", np.array([0 for _ in range(N-1)])) # interlayer biquadratic exchange coupling (J/m^2)

        # defines a function named "name" and calculates its vectorized for, so applies a value on each of the magnetic layers
        def define_funct(self, name, func):
            N = self.outer.number_of_layers
            setattr(self, name, func)
            setattr(self, f"{name}_array", np.array([func(eta) for eta in range(N)]))

        # interlayer function ocouples with adjacent layers only
        def define_interlayer_interaction(self, name, array_fun):
            N = self.outer.number_of_layers
            ic = np.zeros((N, N))

            for i in range(N-1):
                ic[i, i+1] = array_fun[i]
                ic[i+1, i] = array_fun[i]

            setattr(self, f"{name}_array", ic)

        def magnetic_parameters_summary(self):
            params = vars(self)
            summary = {key: value for key, value in params.items() if key.endswith('_array')}
            return summary
        
        


# defines the magnetic tensor and calculates the dynamic matrix considering all magnetic interactions
from numpy import abs, sin, cos, sinh, exp, sign

class MagneticTensor: # calcualte all tensors such that hi = -Nij mj. 
    def __init__(self, film): 

        self.d = film.layer_thickness
        self.N = film.number_of_layers
        self.mag_parmts = film.magnetic_parameters

        self.kr = np.eye(self.N+1)
        self.comps = {0: 'x', 1 : 'y'}
    
    # nu denotes the layer index in which the field is calculated, eta the layer index of the magnetization source 
    
    def matrix(self, N_comp):
        N = self.N
        comps = self.comps
        layers = range(N)

        D = np.zeros((2*N, 2*N), dtype=complex)

    # fills the matrix D with the magnetic interactions and external field.
        for i in layers:
            for comp_i in comps:
                comp_i_ = (comp_i+1)%2 # gives the opposite component of comp_i
                sgn_i = (2*comp_i-1)  # gives -1 if comp is x and +1 if comp is y
                row = 2*i + comp_i
                for j in layers:
                    for comp_j in comps:
                        col = 2*j + comp_j
                        D[row, col] = sgn_i * N_comp(i, j, comp_i_, comp_j)
        
        return D
        
    def D_dip(self, k): # Calculate dipolar interaction tensor Nij
        d = self.d
        kr = self.kr

        th = self.mag_parmts.th_array
        phi = self.mag_parmts.phi_array
        xi = self.mag_parmts.xi_array

        absk = abs(k)
        zeta = 2 * sinh(absk*d/2) * exp(-absk*d/2) / ( d*absk )
        P = 2 * sinh(absk*d/2) **2  / ( d*absk )
        sgn_k = sign(k)
        
        def N_dip(nu,eta,i,j):
            if i==0 and j == 0:
                N_dipolar_V = sin(phi[eta]) * (kr[eta,nu]*sin(phi[eta])*(1-zeta) + (1-kr[eta,nu])*sin(phi[nu])*P* exp(-absk*abs(xi[eta]-xi[nu])))
                N_dipolar_S = 0
            elif i==0 and j == 1:
                N_dipolar_V = cos(phi[eta])*cos(th[eta]) * (kr[eta,nu]*sin(phi[eta])*(1-zeta) + (1-kr[eta,nu])*sin(phi[nu])*P* exp(-absk*abs(xi[eta]-xi[nu])))
                N_dipolar_S = -1j*(1-kr[eta,nu])* sin(th[eta])*sin(phi[nu])*sgn_k*sign(nu-eta)*P* exp(-absk*abs(xi[eta]-xi[nu]))
            elif i==1 and j == 0:
                N_dipolar_V = sin(phi[eta]) * (kr[eta,nu]*cos(phi[eta])*cos(th[eta])*(1-zeta) - (1-kr[eta,nu])*P* exp(-absk*abs(xi[eta]-xi[nu]))*(1j*sgn_k*sign(nu-eta)*sin(th[nu])-cos(phi[nu])*cos(th[nu])))
                N_dipolar_S = 0
            elif i==1 and j == 1:
                N_dipolar_V = cos(phi[eta])*cos(th[eta]) * (kr[eta,nu]*cos(phi[eta])*cos(th[eta])*(1-zeta) - (1-kr[eta,nu])*P* exp(-absk*abs(xi[eta]-xi[nu]))*(1j*sgn_k*sign(nu-eta)*sin(th[nu])-cos(phi[nu])*cos(th[nu])))
                N_dipolar_S = kr[eta,nu]*sin(th[eta])**2*zeta - (1-kr[eta,nu])*sin(th[eta])*P* exp(-absk*abs(xi[eta]-xi[nu]))*(sin(th[nu])+1j*sgn_k*sign(nu-eta)*cos(th[nu])*cos(phi[nu]))
            return self.mag_parmts.Ms_array[nu]*(N_dipolar_V + N_dipolar_S)
        
        # try_f("dip", N_dip, 0)
        
        D_dip = self.matrix(N_dip)
        return D_dip

        

    # exchange interaction considers the taylor expansion of the exchange energy
    def D_ex(self, k):
        d=self.d
        kr = self.kr

        Ms = self.mag_parmts.Ms_array
        Aex = self.mag_parmts.Aex_array
        J = self.mag_parmts.J_array

        th = self.mag_parmts.th_array
        phi = self.mag_parmts.phi_array

        k2 = k**2

        def N_ex(nu, eta,  i, j):
            if i==0 and j == 0:
                N_ex = 2*Aex[nu]/(μ0*Ms[eta]**2)*k2*kr[nu,eta] - (J[nu,eta])/(μ0*Ms[eta]*Ms[nu]*d) *(kr[nu+1,eta]+kr[nu,eta+1])*cos(phi[nu]-phi[eta]) 
            elif i==0 and j == 1:
                N_ex = -(J[nu,eta])/(μ0*Ms[eta]*Ms[nu]*d) *(kr[nu+1,eta]+kr[nu,eta+1])*cos(th[eta])*sin(phi[nu]-phi[eta])
            elif i==1 and j == 0:
                N_ex = (J[nu,eta])/(μ0*Ms[eta]*Ms[nu]*d) *(kr[nu+1,eta]+kr[nu,eta+1])*cos(th[eta])*sin(phi[nu]-phi[eta])
            elif i==1 and j == 1:
                N_ex =  2*Aex[nu]/(μ0*Ms[eta]**2)*k2*kr[nu,eta] - (J[nu,eta])/(μ0*Ms[eta]*Ms[nu]*d) *(kr[nu+1,eta]+kr[nu,eta+1])* (cos(th[nu])*cos(th[eta])*cos(phi[nu]-phi[eta])+sin(th[nu])*sin(th[eta]))
            return Ms[nu]*N_ex 
        
        D_ex = self.matrix(N_ex)
        return D_ex
        
    def D_bqex(self, k):
        d=self.d
        kr = self.kr

        Ms = self.mag_parmts.Ms_array
        Aex = self.mag_parmts.Aex_array
        J = self.mag_parmts.J_array
        J2 = self.mag_parmts.J2_array

        th = self.mag_parmts.th_array
        phi = self.mag_parmts.phi_array

        k2 = k**2

        def N_bqex(nu, eta,  i, j):
            L = 2*J2[nu,eta]/(μ0*Ms[eta]*Ms[nu]*d)
            if i==0 and j == 0:
                k1 = - (cos(th[eta])*cos(th[nu])*cos(phi[eta]-phi[nu]) + cos(2*(phi[eta]-phi[nu]))*sin(th[eta])*sin(th[nu]))
                k2 = - sum((kr[mu+1,nu]+kr[mu,nu+1])*sin(th[mu])**2*sin(phi[mu]-phi[eta])**2 for mu in range(self.N))
                N_bqex = L*(kr[nu+1,eta]+kr[nu,eta+1])*k1 + L*(kr[nu,eta])*k2 
            elif i==0 and j == 1:
                k1 = sin(phi[eta]-phi[nu])* ( cos(2*th[eta])*cos(th[nu]) + cos(phi[eta]-phi[nu])*sin(2*th[eta])*sin(th[nu]) )
                k2 = sum((kr[mu+1,nu]+kr[mu,nu+1])*sin(th[mu])*sin(phi[eta]-phi[mu])*( cos(th[mu])*sin(th[eta]) - cos(th[eta])*sin(th[mu])*cos(phi[eta]-phi[mu]) ) for mu in range(self.N))
                N_bqex = L*(kr[nu+1,eta]+kr[nu,eta+1])*k1 + L*(kr[nu,eta])*k2 
            elif i==1 and j == 0:
                k1 = - sin(phi[eta]-phi[nu])*(cos(th[eta])*cos(2*th[nu]) + cos(phi[eta]-phi[nu])*sin(th[eta])*sin(2*th[nu]))
                k2 = - sum((kr[mu+1,nu]+kr[mu,nu+1])*sin(th[mu])*( -2*cos(th[mu])*sin(th[eta])*sin(phi[eta]-phi[mu]) + cos(th[eta])*sin(th[mu])*sin(2*(phi[eta]-phi[mu]))) for mu in range(self.N))
                N_bqex = L*(kr[nu+1,eta]+kr[nu,eta+1])*k1 + L*(kr[nu,eta])*k2 
            elif i==1 and j == 1:
                k1 = - ( 4*cos(2*th[eta])*cos(2*th[nu])*cos(phi[eta]-phi[nu])+(3+cos(2*(phi[eta]-phi[nu])))*sin(2*th[eta])*sin(2*th[nu]) )
                k2 = - sum((kr[mu+1,nu]+kr[mu,nu+1])*(cos(th[mu])*cos(th[eta])-cos(th[eta])*sin(th[mu])*cos(phi[eta]-phi[mu]))**2 for mu in range(self.N))
                N_bqex = L*(kr[nu+1,eta]+kr[nu,eta+1])*k1 + L*(kr[nu,eta])*k2 
            return Ms[nu]*N_bqex 
        
        # try_f("ex", N_ex, 0)
        
        D_bqex = self.matrix(N_bqex)
        return D_bqex
    
    def D_u(self, k):
        d=self.d
        kr = self.kr

        Ms = self.mag_parmts.Ms_array
        Aex = self.mag_parmts.Aex_array

        th = self.mag_parmts.th_array
        phi = self.mag_parmts.phi_array

        k2 = k**2

        def N_u(nu, eta,  i, j):
            if i==0 and j == 0:
                Nu = 0 
            elif i==0 and j == 1:
                Nu = 0
            elif i==1 and j == 0:
                Nu = 0
            elif i==1 and j == 1:
                Nu = 0
            return Ms[nu]*Nu 
        
        # try_f("ex", N_ex, 0)
        
        D_ex = self.matrix(N_u)
        return D_ex
    
    
    def D_bdmi(self, k):
        d=self.d
        kr = self.kr

        Db = self.mag_parmts.Db_array
        Db_inter = self.mag_parmts.Db_inter_array
        Ms = self.mag_parmts.Ms_array

        th = self.mag_parmts.th_array
        phi = self.mag_parmts.phi_array
        
        def N_bdmi(nu, eta, i, j):
            D_inter = (Db_inter[nu,eta]/(μ0*Ms[eta]*Ms[nu]) )* (kr[nu+1,eta]-kr[nu,eta+1]) 
            D_intra = (2*Db[nu] *1j * k/(μ0*Ms[eta]*Ms[nu]) ) * kr[nu,eta]
            if i==0 and j == 0:
                N_dmi = -D_inter * sin(phi[eta]-phi[nu])
            elif i==0 and j == 1:
                N_dmi = -D_intra*sin(th[eta]) - D_inter*cos(th[eta])*cos(phi[eta]-phi[nu])
            elif i==1 and j == 0:
                N_dmi = D_intra*cos(phi[eta]) + D_inter*cos(th[nu])*cos(phi[eta]-phi[nu])
            elif i==1 and j == 1:
                N_dmi = -D_intra*cos(th[eta])*sin(phi[eta]) - D_inter*cos(th[nu])*cos(th[eta])*sin(phi[eta]-phi[nu])
            return Ms[nu]*N_dmi
        
        # try_f("ex", N_bdmi, 1)

        D_bdmi = self.matrix(N_bdmi)
        return D_bdmi
    
    def D_Heq(self, H0, th_H0, phi_H0):
        d=self.d
        kr = self.kr

        Ms = self.mag_parmts.Ms_array
        Aex = self.mag_parmts.Aex_array
        J = self.mag_parmts.J_array
        Ku = self.mag_parmts.Ku_array
        Db = self.mag_parmts.Db_array

        th_u = self.mag_parmts.th_u_array
        phi_u = self.mag_parmts.phi_u_array

        th = self.mag_parmts.th_array
        phi = self.mag_parmts.phi_array
        
        def N_Heq(nu, eta, i, j):
            Hext = H0*( sin(th[nu])*sin(th_H0)*cos(phi[nu]-phi_H0) + cos(th[nu])*cos(th_H0) )
            Hu = 2*Ku[nu]/(μ0*Ms[nu])* (sin(th[nu])*sin(th_u[nu])*cos(phi[nu]-phi_u[nu]) )**2
            Hdip = -Ms[nu] * cos(th[nu]) 
            Hex = sum(J[nu,p]/(μ0*Ms[nu]*d) * (kr[nu+1,p]+kr[nu,p+1]) * (cos(th[nu])*cos(th[p])+cos(phi[nu]-phi[p])*sin(th[nu])*sin(th[p])) for p in range(self.N))
            Hdmi = sum(Db[nu]/(μ0*Ms[nu]*d) * (kr[nu+1,p]-kr[nu,p+1]) * (sin(th[nu]) * sin(th[p]) * sin(phi[p]-phi[nu])) for p in range(self.N))
            if i==0 and j == 0:
                H = kr[nu,eta]*(Hext + Hu + Hdip + Hex + Hdmi)
            elif i==0 and j == 1:
                H = 0
            elif i==1 and j == 0:
                H = 0
            elif i==1 and j == 1:
                H = kr[nu,eta]*(Hext + Hu + Hdip + Hex + Hdmi)
            return H 
        
        # try_f("Heq", N_Heq, 0)
        
        D_Heq = self.matrix(N_Heq)
        return D_Heq



# with the matrix, calculates the dispersions and eigenmodes 
from numpy import zeros, linalg, real_if_close, where, real, vstack
from numpy import argsort, where
from joblib import Parallel, delayed

class Dispersion:  # class to calculate the magnetic interactions
    def __init__(self, film):
        self.film = film
        self.Nmag = MagneticTensor(film)
        self.mag_parmts = film.magnetic_parameters
        
    '''wv Calculates the eigenvalues and eigenvectors of the dynamic matrix D for a given wavevector k and external magnetic field Hext.'''    
    def wv(self, k, Hext, th_H0, phi_H0): 
        mag_parmts = self.mag_parmts
        Nmag = self.Nmag

        scale_factor = -1j*mag_parmts.gamma(0) * μ0
        D_tot = scale_factor*(Nmag.D_dip(k) + Nmag.D_ex(k) + Nmag.D_bdmi(k) + Nmag.D_bqex(k) + Nmag.D_Heq(Hext, th_H0, phi_H0))
            
        w, v = linalg.eig(D_tot)

        return w, v

    '''freq_ordering takes the eigenvalues and eigenvectors and orders them by abs value. this could be changed to order them by parity'''
    def freq_ordering(self, k, Hext, th_H0, phi_H0): 
        w, v = self.wv(k, Hext, th_H0, phi_H0)

        # Step 1: Filter positive eigenvalues and their indices
        positive_indices = where(real(w) >= 0)[0]

        # Step 2: Extract positive eigenvalues and corresponding eigenvectors
        w_pos = w[positive_indices]
        v_pos = v[:, positive_indices]

        # Step 3: Sort positive eigenvalues and reorder eigenvectors accordingly
        sorted_indices = argsort(w_pos)
        
        w_sorted = w_pos[sorted_indices]
        v_sorted = v_pos[:, sorted_indices]

        return w_sorted, v_sorted
    
    def w_vect(self, k_list, Hext, th_H0, phi_H0):
        def compute_for_k(k):
            w, v = self.freq_ordering(1e6 * k, Hext, th_H0, phi_H0)
            w_ghz = real_if_close(1e-9 * w / (2 * pi), tol=1e-6)
            return k, w_ghz, v

        # Run in parallel across all k values
        results = Parallel(n_jobs=-1)(delayed(compute_for_k)(k) for k in k_list)

        w_vect = []
        v_dic = {}

        for k, w, v in results:
            w_vect.append(w)
            v_dic[int(k)] = v

        # Filter to ensure consistent mode count
        mode_count = len(w_vect[0])
        w_vect = [w for w in w_vect if len(w) == mode_count]
        w_vect = vstack(w_vect)

        if len(w_vect) != len(k_list):
            print("Unstable equilibrium state.")

        return w_vect, v_dic
    # gives w_vect[k_index, mode_index] and v_dic[k][ :, mode_index]
    


# %matplotlib widget

from numpy import argmin
import matplotlib.pyplot as plt
import os



class Graphs:  # class to calculate the magnetic interactions
    def __init__(self):
        pass

    def dispersion(self, k_vals, w_list_vals, modes = [0], xlim = [-200, 200], ylim = [0,50], legends_list = None, v_k_mode= None):
        legend_counter = 0
        
        plt.figure(figsize=(8,5))
        
        for w_vals in w_list_vals: # plots frequencies form a list of w_vals arrays
            for mode in modes:
                plt.plot(k_vals, w_vals[:, mode], label = f"mode {mode} ({legends_list[legend_counter]})" if legends_list else f"mode {mode}")
            if v_k_mode:
                v_k, v_mode = v_k_mode
                k_index = argmin(abs(k_vals - v_k))
                plt.plot(*[v_k,w_vals[:,v_mode][k_index]], 'ro')  # 'r' for red, 'o' for circle marker
            legend_counter +=1
 

        plt.xlabel('k (1/μm)')
        plt.ylabel('f (GHz)')
        plt.title('f vs k')
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.legend()
        plt.grid(True)
        plt.show()

    def graph(self, x_vals, w_list_vals, modes = [0], xlim = [0, 1], ylim = [0,100], legends_list = None, v_k_mode= None, x_label = 'x variable'):
        legend_counter = 0
        
        plt.figure(figsize=(8,5))
        
        for w_vals in w_list_vals: # plots frequencies form a list of w_vals arrays
            for mode in modes:
                plt.plot(x_vals, w_vals[:, mode] , 'o', label = f"mode {mode} ({legends_list[legend_counter]})" if legends_list else f"mode {mode}")
            if v_k_mode:
                v_k, v_mode = v_k_mode
                k_index = argmin(abs(x_vals - v_k))
                plt.plot(*[v_k,w_vals[:,v_mode][k_index]], 'ro')  # 'r' for red, 'o' for circle marker
            legend_counter +=1
 

        plt.xlabel(x_label)
        plt.ylabel('f (GHz)')
        plt.title('f vs k')
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.legend()
        plt.grid(True)
        plt.show()

    def orbits(self, k, mode, v_vals, title = None, t0 = .1): # plots the orbits for a given k and mode
        v = v_vals[k][:,mode] 
        N = len(v) // 2
        t = np.linspace(0, 1, 200)
        arg = 1j*2*pi*t

        kap = 2/N
        desired_axes_size = (2*1.2, 2*N*kap)
        margin = 0.25  # inches

        fig_width = desired_axes_size[0] + 2 * margin
        fig_height = desired_axes_size[1] + 2 * margin

        fig = plt.figure(figsize=(fig_width, fig_height))
        ax = fig.add_axes([
            margin / fig_width,
            margin / fig_height,
            desired_axes_size[0] / fig_width,
            desired_axes_size[1] / fig_height
        ])

        for n in range(N):
            mx = v[2*n]
            my = v[2*n + 1]

            #print('layer', n, 'mx = ', round(real(((mx)*(np.conjugate(mx)))**(1/2)),4), ' my = ', round(real(((my)*(np.conjugate(my)))**(1/2)),4))
            eps = round(real(((mx/my)*(np.conjugate(mx/my)))**(1/2)),3)
            #print('ellipticity = ', eps)

            mx0 = real(mx*exp(1j*2*pi*t0))
            my0 = real(my*exp(1j*2*pi*t0)/kap)

            mx_t = real(exp(arg) * mx)
            my_t = real(exp(arg) * my/kap)

            y_offset = 1 * (n + 0.5)
            plt.plot(mx_t, my_t + y_offset, color='grey') #label=f"epsilon = {eps} to see how the orbits shape"
            plt.plot(mx0, my0 + y_offset, 'ro', ms=3, label = f"t0 = {t0}" if n==0 else None )  # 'r' for red, 'o' for circle marker

        plt.title(title)
        plt.xlabel(r"$m_x$")
        plt.ylabel(r"$m_y$")
        plt.legend()

        plt.xticks([])  # Remove x-axis numbers
        plt.yticks([])
        plt.grid(False)
        # plt.legend()
        plt.xlim(-.6,.6)
        plt.ylim(-.5,N+.5)
        plt.show()   

    def dmi(self, films_list): # plots the DMI profile for a list of films
        
        plt.figure(figsize=(6,4))
        max_value = 0
        min_value = 0
        for film in films_list:
            y_values = np.linspace(0, film.thickness*1e9, 50)
            mid_points = np.arange(film.layer_thickness*1e9/2, film.thickness*1e9, film.layer_thickness*1e9)
           
            D_vals_midpoint = 1e3*film.D_profile(mid_points*1e-9)
            D_vals = 1e3*film.D_profile(y_values*1e-9)
            max_value = max(max(D_vals),max_value)
            min_value = min(min(D_vals),min_value)
            
            plt.plot(D_vals_midpoint,mid_points, 'ro', ms=3)
            plt.plot(D_vals,y_values)
        
        film = films_list[0]
        N = film.number_of_layers
        for n in range(N):
            plt.hlines(film.layer_thickness*n*1e9, min_value, max_value, colors='gray', linestyles='dashed', linewidth=0.5)
        plt.ylabel('Position across thickness (nm)')
        plt.xlabel('D (mJ/m²)')
        plt.title('DMI Profile')
        plt.xlim(min_value,max_value)
        plt.ylim(0,film.thickness*1e9)
        plt.grid(False)
        # plt.legend()
        plt.show()