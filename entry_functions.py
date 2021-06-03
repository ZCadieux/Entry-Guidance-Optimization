import numpy as np
import omtools.api as ot

#change to same format, with initialize and setup functions
class EntryFunctions(ot.Group):
    def initialize(self):
        self.options.declare('num_nodes', default=1, types=int)
        self.options.declare('R0', default=1., types=(int, float))
        self.options.declare('flatten', default=1., types=(int, float))


    def setup(self):
        #only necessary if using geodetic
        #R0 = self.options['R0']
        #flatten = self.options['flatten']
        num = self.options['num_nodes']
        area = 159.94 #from vehicle data

        #instead of as parameters to setup function, should this be declared here?
        r = self.declare_input('r', shape=(num,1)) #radial distance from center of Mars
        V = self.declare_input('V', shape=(num,1)) #Mars relative velocity

        #input r, V; output vector with density, Vsound, lift, drag
        #Step 1: Convert to geodetic, to get altitude (r and height in km)
        #in beginning, using just geocentric
        #Req = R0/1000
        #rinv = Req/r
        #height = (r-Req)+Req*flatten*(0.5*(1.0-ot.cos(2*r)) + ((0.25*rinv-0.0625)*(1.0-ot.cos(4*r))*flatten))
        #update with geocentric assumption, assuming r = height

        #Step 2, use fit polynomial for density
        A = -4.422345201370313; B = 0.1130629761317862
        C = 0.0003003890333711271; D = -2.148583579589203e-5
        E = 7.041319756467658e-8; F = -0.04399163896948651
        G = 0.0006555968430113561; H = -3.709114148189494e-6
        I = 8.137552734673534e-9
        density = np.exp((A+r*(B+r*(C+r*(D+r*E))))/(1.0+r*(F+r*(G+r*(H+r*I))))) #density

        #Step 3, use fit polynomial for Vsound
        A = 3.90812820e2; B = -5.07169986e3
        C = 6.00818111e4; D = -3.71909017e5
        E = 1.27052434e6; F = -2.47999976e6
        G = 2.74880412e6; H = -1.60759656e6
        I = 3.84982694e5
        if r < 10000:
            x  = 10000.
        elif r > 130000:
            x = 130000.
        else:
            x = r
        x = x*1000/130081.410992
        Vsound = A + x*(B + x*(C + x*(D + x*(E + x*(F + x*(G + x*(H + I*x))))))) #speed of sound

        #Step 4, find CL, CD
        x = V/Vsound #find mach number
        if x < 2:
            x = 2
        elif x > 12:
            x = 12
        A = -1.718628628451222; B = -0.1377158161003165
        C = 1.083905623513791; D = 7.069525902496836; E = -10.24433210397554
        C_L = A + B*x + C*np.sqrt(x) + D*x**(-3/2) + E*np.exp(-x)

        xinv = 1./x; A = 1.015355733524931; B = -0.04826676756876537
        C = 1.381399436898735; D = -0.9962918706041476
        C_D = A + xinv*(B + xinv*(C + xinv*D))

        coeff = density * V**2 * area/2

        lift = C_L * coeff #lift
        drag = C_D * coeff #drag

        self.register_output('lift', lift)
        self.register_output('drag', drag)
        #instead of as a return, should this be outputted here?