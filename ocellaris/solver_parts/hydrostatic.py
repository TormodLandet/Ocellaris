import dolfin
from ocellaris.utils import ocellaris_error


class HydrostaticPressure(object):
    def __init__(self, rho, g, ph, zero_level, eps=1e-8):
        """
        Calculate the hydrostatic pressure

        The gravity vector g *must* be parallel to one of the axes
        """
        Vp = ph.function_space()
        p = dolfin.TrialFunction(Vp)
        q = dolfin.TestFunction(Vp)
        
        directions = set()
        for i, gi in enumerate(g.values()):
            if gi != 0:
                directions.add(i)
        
        if len(directions) == 0:
            self.active = False
            return
        elif len(directions) > 1:
            ocellaris_error('Error calculating hydrostatic pressure',
                            'Gravity vector %r is not parallel to an axis'
                            % g.values())
        
        self.active = True
        d = directions.pop()
        
        a = p.dx(d)*q.dx(d)*dolfin.dx
        L = g[d]*rho*q.dx(d)*dolfin.dx
        
        inside = lambda  x, on_boundary: zero_level - eps <= x[d] <= zero_level + eps
        self.zero_bc = dolfin.DirichletBC(Vp, 0.0, inside)
        self.func = ph
        self.tensor_lhs = dolfin.assemble(a)
        self.form_rhs = L
    
    def update(self):
        if not self.active:
            self.func.zero()
            return
        
        t = dolfin.Timer('Ocellaris update hydrostatic pressure')
        
        A = self.tensor_lhs
        b = dolfin.assemble(self.form_rhs)
        self.zero_bc.apply(A, b)
        dolfin.solve(A, self.func.vector(), b)
        
        t.stop()
