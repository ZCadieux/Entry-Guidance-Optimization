import numpy as np
import omtools.api as ot

class EntrySystem(ot.Group):
    def initialize(self):
        #declare non-state values, including constants and L, D, etc
        #L, D, sigma, Omega (at minimum)
        #what is number of nodes? are commented values necessary?
        self.options.declare('num_nodes', default=1, types=int)
        self.options.declare('R0', default=1., types=(int, float))
        self.options.declare('g0', default=1., types=(int, float))
        self.options.declare('vex', default=1., types=(int, float))
        #self.options.declare('w2', default=1., types=(int, float))
        self.options.declare('m_scale', default=1., types=(int, float))

        self.options.declare('D', default=1., types=(int, float))
        self.options.declare('L', default=1., types=(int, float))
        self.options.declare('sigma', default=1., types=(int, float))
        self.options.declare('flatten', default=1., types=(int, float))
        self.options.declare('density', default=1., types=(int, float))
        self.options.declare('Vsound', default=1., types=(int, float))

    def setup(self):
        #define values from initialize
        num = self.options['num_nodes']
        R0 = self.options['R0']
        g0 = self.options['g0']
        vex = self.options['vex']
        #w2 = self.options['w2']
        m_scale = self.options['m_scale']

        D = self.options['D'] #bring these in from function?
        L = self.options['L']
        Omega = self.options['Omega']
        #flatten = self.options['flatten']
        #density = self.options['density']
        #Vsound = self.options['Vsound']

        #declare state variables (inputs)
        #entry state variables would be r, V, gamma, psi, phi, theta
        r = self.declare_input('r', shape=(num,1)) #radial distance from center of Mars
        theta = self.declare_input('theta', shape=(num,1)) #longitude
        phi = self.declare_input('phi', shape=(num,1)) #latitude
        V = self.declare_input('V', shape=(num,1)) #Mars relative velocity
        gamma = self.declare_input('gamma', shape=(num,1)) #flight-path angle of Mars relative velocity
        psi = self.declare_input('psi', shape=(num,1)) #heading angle of Mars relative velocity

        #mass value (still necessary? carried over from ascent problem)
        m = self.declare_input('m', shape=(num,1))

        #control variables, unit vector and thrust value
        #in entry case, control variable is sigma, bank angle
        #ux = self.declare_input('ux', shape=(num,1))
        #uy = self.declare_input('uy', shape=(num,1))
        #uz = self.declare_input('uz', shape=(num,1))
        T = self.declare_input('T', shape=(num,1))
        sigma = self.declare_input('sigma', shape =(num, 1))


        #define state equations
        #based on dynamics
        dr_dt = V * ot.sin(gamma) #Equ 12
        dtheta_dt = V * ot.cos(gamma) * ot.sin(psi) / (r * ot.cos(phi)) #Equ 13
        dphi_dt = V * ot.cos(gamma) * ot.cos(psi) / r #Equ 14
        dV_dt = -D - (ot.sin(gamma) / r^2) + Omega^2 * r * ot.cos(phi) * \
            (ot.sin(gamma) * ot.cos(phi) - ot.cos(gamma) * ot.sin(phi) * ot.cos(psi)) #Equ 15
        dgamma_dt = (1./V)*(L * ot.cos(sigma) + (V^2 - 1./r)*(ot.cos(gamma)/r) + \
            2. * Omega * V*ot.cos(phi)*ot.sin(psi) + Omega^2*r*ot.cos(phi) * \
            (ot.cos(gamma)*ot.cos(phi) + ot.sin(gamma)*ot.cos(psi)*ot.sin(phi))) #Equ 16
        dpsi_dt = (1./V)*(L*ot.sin(sigma)/ot.cos(gamma) + (V^2/r)*ot.cos(gamma)* \
            ot.sin(psi)*ot.tan(phi) - 2.*Omega*V*(ot.tan(gamma)*ot.cos(psi)*ot.cos(phi) - \
            ot.sin(phi)) + (Omega^2*r/ot.cos(gamma))*ot.sin(psi)*ot.sin(phi)*ot.cos(phi)) #Equ 17
        
        #mass flow
        dm_dt = - (T / (m_scale * vex)) * np.sqrt(R0/g0)

        #define outputs (derivatives)
        #time derivatives of state variables
        self.register_output('dr_dt', dr_dt)
        self.register_output('dtheta_dt', dtheta_dt)
        self.register_output('dphi_dt', dphi_dt)
        self.register_output('dV_dt', dV_dt)
        self.register_output('dgamma_dt', dgamma_dt)
        self.register_output('dpsi_dt', dpsi_dt)

        self.register_output('dm_dt', dm_dt)

    def varEval(r, V, self):
        #correct setup and placement? Use of parameters?
        R0 = self.options['R0']
        flatten = self.options['flatten']
        area = 159.94

        #input r, V; output vector with density, Vsound, lift, drag
        #Step 1: Convert to geodetic, to get altitude (r and height in km)
        Req = R0/1000
        rinv = Req/r
        height = (r-Req)+Req*flatten*(0.5*(1.0-ot.cos(2*r)) + ((0.25*rinv-0.0625)*(1.0-ot.cos(4*r))*flatten))

        #Step 2, use fit polynomial for density
        A = -4.422345201370313; B = 0.1130629761317862
        C = 0.0003003890333711271; D = -2.148583579589203*10^-5
        E = 7.041319756467658*10^-8; F = -0.04399163896948651
        G = 0.0006555968430113561; H = -3.709114148189494*10^-6
        I = 8.137552734673534*10^-9
        vars[0] = ot.exp((A+r(B+r(C+r(D+r(E)))))/(1.0+r(F+r(G+r(H+r(I)))))) #density

        #Step 3, use fit polynomial for Vsound
        A = 3.90812820*10^2; B = -5.07169986*10^3
        C = 6.00818111*10^4; D = -3.71909017*10^5
        E = 1.27052434*10^6; F = -2.47999976*10^6
        G = 2.74880412*10^6; H = -1.60759656*10^6
        I = 3.84982694*10^5
        if height < 10000:
            x  = 10000
        elif height > 130000:
            x = 130000
        else:
            x = height
        x = x*1000/130081.410992
        vars[1] = A + x(B + x(C + x(D + x(E + x(F + x(G + x(H + I*x))))))) #speed of sound

        #Step 4, find CL, CD
        x = V/vars[1] #find mach number
        if x < 2:
            x = 2
        elif x > 12:
            x = 12
        A = -1.718628628451222; B = -0.1377158161003165
        C = 1.083905623513791; D = 7.069525902496836; E = -10.24433210397554
        C_L = A + B*x + C*ot.sqrt(x) + D*x^(-3/2) + E*ot.exp(-x)

        x = 1/x, A = 1.015355733524931; B = -0.04826676756876537
        C = 1.381399436898735; D = -0.9962918706041476
        C_D = A + x(B + x(C + x(D)))

        coeff = vars[1] * V^2 * area/2

        vars[2] = C_L * coeff #lift
        vars[3] = C_D * coeff #drag

        return vars