import numpy as np
import omtools.api as ot
#CURRENT CODE COPIED FROM ASCENT SYSTEM
class EntrySystem(ot.Group):
    def initialize(self):
        #declare non-state values, including constants and L, D, etc
        #L, D, sigma, Omega (at minimum)
        self.options.declare('num_nodes', default=1, types=int)
        self.options.declare('R0', default=1., types=(int, float))
        self.options.declare('g0', default=1., types=(int, float))
        self.options.declare('vex', default=1., types=(int, float))
        self.options.declare('w2', default=1., types=(int, float))
        self.options.declare('m_scale', default=1., types=(int, float))

    def setup(self):
        #define values from initialize
        num = self.options['num_nodes']
        g0 = self.options['g0']
        R0 = self.options['R0']
        vex = self.options['vex']
        w2 = self.options['w2']
        m_scale = self.options['m_scale']

        #declare state variables (inputs)
        #entry state variables would be r, V, gamma, psi, phi, theta
        rx = self.declare_input('rx', shape=(num,1))
        ry = self.declare_input('ry', shape=(num,1))
        rz = self.declare_input('rz', shape=(num,1))
        Vx = self.declare_input('Vx', shape=(num,1))
        Vy = self.declare_input('Vy', shape=(num,1))
        Vz = self.declare_input('Vz', shape=(num,1))
        m = self.declare_input('m', shape=(num,1))
        ux = self.declare_input('ux', shape=(num,1))
        uy = self.declare_input('uy', shape=(num,1))
        uz = self.declare_input('uz', shape=(num,1))
        T = self.declare_input('T', shape=(num,1))

        #define state equations
        #based on dynamics
        drx_dt = 1. * Vx
        dry_dt = 1. * Vy
        drz_dt = 1. * Vz
        dVx_dt = - w2 * rx + (T * ux) / (m * m_scale * g0)
        dVy_dt = - w2 * ry + (T * uy) / (m * m_scale * g0)
        dVz_dt = - w2 * rz + (T * uz) / (m * m_scale * g0)
        dm_dt = - (T / (m_scale * vex)) * np.sqrt(R0/g0)

        # Thrust acceleration used instead of thrust
        # dVx_dt = - w2 * rx + (T * ux) / (m * g0)
        # dVy_dt = - w2 * ry + (T * uy) / (m * g0)
        # dVz_dt = - w2 * rz + (T * uz) / (m * g0)
        # dm_dt = - (T / vex) * np.sqrt(R0/g0) #* (ux**2 + uy**2 + uz**2)

        #define outputs (derivatives)
        #time derivatives of state variables
        self.register_output('drx_dt', drx_dt)
        self.register_output('dry_dt', dry_dt)
        self.register_output('drz_dt', drz_dt)
        self.register_output('dVx_dt', dVx_dt)
        self.register_output('dVy_dt', dVy_dt)
        self.register_output('dVz_dt', dVz_dt)
        self.register_output('dm_dt', dm_dt)
        # self.create_indep_var('dm_dt',val=dm_dt, shape=(num,1))