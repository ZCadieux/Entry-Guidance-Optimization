#Example code, not used in final version

from openmdao.api import Problem
from omtools.api import Group
import omtools.api as ot
import numpy as np

class ProblemGroup(Group):
    def initialize(self):
        self.options.declare('r')

    def setup(self):
        r = self.options['r']

        x = self.declare_input('x', val = 1.)
        y = self.declare_input('y', val = 1.)

        z = r*(x+y)

        self.register_output('z', z)