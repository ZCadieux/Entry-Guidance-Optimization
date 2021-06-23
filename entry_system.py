import numpy as np
import omtools.api as ot
#from entry_functions import EntryFunction as EF

class EntrySystem(ot.Group):
    def initialize(self):
        #declare non-state values, including constants and L, D, etc
        #L, D, sigma, Omega (at minimum)
        #what is number of nodes? are commented values necessary?
        self.options.declare('num_nodes', default=1, types=int)
        self.options.declare('D', default=1., types=(int, float))
        self.options.declare('L', default=1., types=(int, float))
        self.options.declare('Omega', default=1., types=(int, float))

    def setup(self):
        num = self.options['num_nodes']
        Omega = self.options['Omega']
        D = self.options['D']
        L = self.options['L']

        #declare state variables (inputs)
        #entry state variables would be r, V, gamma, psi, phi, theta
        r = self.declare_input('r', shape=(num,1)) #radial distance from center of Mars
        theta = self.declare_input('theta', shape=(num,1)) #longitude
        phi = self.declare_input('phi', shape=(num,1)) #latitude
        V = self.declare_input('V', shape=(num,1)) #Mars relative velocity
        gamma = self.declare_input('gamma', shape=(num,1)) #flight-path angle of Mars relative velocity
        psi = self.declare_input('psi', shape=(num,1)) #heading angle of Mars relative velocity

        #control variables, unit vector and thrust value
        #in entry case, control variable is sigma, bank angle
        sigma = self.declare_input('sigma', shape =(num, 1))

        #define values from initialize REDO THIS
        #vars =  EF.setup(r,V,self)
        #D = vars[3]
        #L = vars[2]

        #define state equations
        #based on dynamics
        dr_dt = V * ot.sin(gamma) #Equ 12
        dtheta_dt = V * ot.cos(gamma) * ot.sin(psi) / (r * ot.cos(phi)) #Equ 13
        dphi_dt = V * ot.cos(gamma) * ot.cos(psi) / r #Equ 14
        dV_dt = -D - (ot.sin(gamma) / r**2) + Omega**2 * r * ot.cos(phi) * \
            (ot.sin(gamma) * ot.cos(phi) - ot.cos(gamma) * ot.sin(phi) * ot.cos(psi)) #Equ 15
        dgamma_dt = (1./V)*(L * ot.cos(sigma) + (V**2 - 1./r)*(ot.cos(gamma)/r) + \
            2. * Omega * V*ot.cos(phi)*ot.sin(psi) + Omega**2*r*ot.cos(phi) * \
            (ot.cos(gamma)*ot.cos(phi) + ot.sin(gamma)*ot.cos(psi)*ot.sin(phi))) #Equ 16
        dpsi_dt = (1./V)*(L*ot.sin(sigma)/ot.cos(gamma) + (V**2/r)*ot.cos(gamma)* \
            ot.sin(psi)*ot.tan(phi) - 2.*Omega*V*(ot.tan(gamma)*ot.cos(psi)*ot.cos(phi) - \
            ot.sin(phi)) + (Omega**2*r/ot.cos(gamma))*ot.sin(psi)*ot.sin(phi)*ot.cos(phi)) #Equ 17
        

        #define outputs (derivatives)
        #time derivatives of state variables
        self.register_output('dr_dt', dr_dt)
        self.register_output('dtheta_dt', dtheta_dt)
        self.register_output('dphi_dt', dphi_dt)
        self.register_output('dV_dt', dV_dt)
        self.register_output('dgamma_dt', dgamma_dt)
        self.register_output('dpsi_dt', dpsi_dt)
