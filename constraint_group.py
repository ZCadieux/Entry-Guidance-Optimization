import numpy as np
import omtools.api as ot
from numpy import linalg

class OrbitalConditionOne(ot.Group):
    def initialize(self):
        self.options.declare('a', default=1., types=float)
        self.options.declare('ecc', default=1., types=float)

    def setup(self):
        a = self.options['a']
        ecc = self.options['ecc']

        rx = self.declare_input('rx', shape=(1,))
        ry = self.declare_input('ry', shape=(1,))
        rz = self.declare_input('rz', shape=(1,))
        Vx = self.declare_input('Vx', shape=(1,))
        Vy = self.declare_input('Vy', shape=(1,))
        Vz = self.declare_input('Vz', shape=(1,))

        r = self.create_output('r', shape=(3,))
        V = self.create_output('V', shape=(3,))

        # print('REVIEW:')
        # print(rx)
        # print(r)

        r[0] = rx
        r[1] = ry
        r[2] = rz
        V[0] = Vx
        V[1] = Vy
        V[2] = Vz

        rxV = ot.cross(r,V,axis=0)

        C1 = ot.dot(rxV,1*rxV,axis=0) - a*(1 - ecc**2)
        # C1 = ot.pnorm(rxV,rxV) - a*(1 - ecc**2)

        self.register_output('C1', C1)


class OrbitalConditionTwo(ot.Group):

    def initialize(self):
        self.options.declare('a', default=1., types=float)

    def setup(self):
        a = self.options['a']

        rx = self.declare_input('rx')
        ry = self.declare_input('ry')
        rz = self.declare_input('rz')
        Vx = self.declare_input('Vx')
        Vy = self.declare_input('Vy')
        Vz = self.declare_input('Vz')

        r = self.create_output('r', shape=(3,))
        V = self.create_output('V', shape=(3,))

        r[0] = rx
        r[1] = ry
        r[2] = rz
        V[0] = Vx
        V[1] = Vy
        V[2] = Vz

        #r_norm = (rx**2 + ry**2 + rz**2)**0.5
        #V_norm = (Vx**2 + Vy**2 + Vz**2)**0.5

        # C2 = (((ot.sum(V**2))**0.5)**2)/2. - 1./((ot.sum(r**2))**0.5) + 1./(2*a)
        C2 = ((ot.pnorm(V))**2)/2. - 1./(ot.pnorm(r)) + 1./(2*a)
        #C2 = (V_norm**2)/2. - 1./(r_norm) + 1./(2*a)

        self.register_output('C2', C2)


class OrbitalConditionThree(ot.Group):

    def initialize(self):
        self.options.declare('inc', default=1., types=float)

    def setup(self):
        inc = self.options['inc']

        rx = self.declare_input('rx')
        ry = self.declare_input('ry')
        rz = self.declare_input('rz')
        Vx = self.declare_input('Vx')
        Vy = self.declare_input('Vy')
        Vz = self.declare_input('Vz')

        # z = self.create_indep_var('z', shape=(3,))
        # r = self.create_output('r', shape=(3,))
        # V = self.create_output('V', shape=(3,))

        # z = [0, 0, 1]
        #
        # r[0] = rx
        # r[1] = ry
        # r[2] = rz
        #
        # V[0] = Vx
        # V[1] = Vy
        # V[2] = Vz

        # rxV = ot.cross(r,V,axis=0)

        C3 = rx*Vy - ry*Vx - np.cos(inc)*((rx*Vy - ry*Vx)**2 \
            + (rx*Vz - rz*Vx)**2 + (ry*Vz - rz*Vy)**2)**(0.5)

        #C3 = ot.dot(z, rxV, axis=1) - ot.pnorm(rxV) * np.cos(inc)

        self.register_output('C3', C3)


class ThrustVectorConstraint(ot.Group):
    def initialize(self):
        self.options.declare('num_nodes', default=1, types=int)

    def setup(self):
        num = self.options['num_nodes']

        ux = self.declare_input('ux', shape=(num,1))
        uy = self.declare_input('uy', shape=(num,1))
        uz = self.declare_input('uz', shape=(num,1))

        CT = ux**2 + uy**2 + uz**2

        self.register_output('CT', CT)
