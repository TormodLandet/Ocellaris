# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from dolfin import Constant
from ocellaris.utils import timeit, ocellaris_interpolate
from . import Solver, register_solver


@register_solver('AnalyticalSolution')
class AnalyticalSolution(Solver):
    def __init__(self, simulation):
        """
        A Navier-Stokes "solver" using only analytical expressions to compute u and p
        """
        self.simulation = sim = simulation
        self.create_functions()

        # Solver control parameters
        sim.data['dt'] = Constant(simulation.dt)

        # Store number of iterations
        self.niters = None

    def create_functions(self):
        """
        Create functions to hold solutions
        """
        sim = self.simulation

        # Function spaces
        Vu = sim.data['Vu']
        Vp = sim.data['Vp']

        # Create segregated functions on component and vector form
        u_list, up_list, upp_list = [], [], []
        for d in range(sim.ndim):
            sim.data['u%d' % d] = u = dolfin.Function(Vu)
            sim.data['up%d' % d] = up = dolfin.Function(Vu)
            sim.data['upp%d' % d] = upp = dolfin.Function(Vu)
            u_list.append(u)
            up_list.append(up)
            upp_list.append(upp)
        sim.data['u'] = dolfin.as_vector(u_list)
        sim.data['up'] = dolfin.as_vector(up_list)
        sim.data['upp'] = dolfin.as_vector(upp_list)
        sim.data['p'] = dolfin.Function(Vp)
        sim.data['u_conv'] = sim.data['u']

    @timeit
    def apply_analytical_solution(self):
        """
        Read initial conditions for up0, up1 and p to update the fields
        """
        sim = self.simulation

        funcs = {
            'up0': sim.data['u0'] if 'u0' in sim.data else None,
            'up1': sim.data['u1'] if 'u1' in sim.data else None,
            'up2': sim.data['u2'] if 'u2' in sim.data else None,
            'p': sim.data['p'],
        }
        ic = sim.input.get_value('initial_conditions', {}, 'dict(string:dict)')
        for name, info in ic.items():
            name = str(name)

            if name not in funcs:
                continue
            func = funcs[name]

            if 'cpp_code' in info:
                # Initial condition given as a C++ code string
                cpp_code = str(info['cpp_code'])
                V = func.function_space()
                description = 'initial conditions for %r' % name

                # Update the function by running the C++ code
                ocellaris_interpolate(sim, cpp_code, description, V, func)

            else:
                # Initial condition given as a known field function
                field_name, func_name = info['function'].strip().split('/')
                field = self.simulation.fields[field_name]
                func.interpolate(field.get_variable(func_name))

    @timeit
    def run(self):
        """
        Run the simulation
        """
        sim = self.simulation
        sim.hooks.simulation_started()
        t = sim.time
        it = sim.timestep

        while True:
            # Get input values, these can possibly change over time
            dt = sim.input.get_value('time/dt', required_type='float')
            tmax = sim.input.get_value('time/tmax', required_type='float')

            # Check if the simulation is done
            if t + dt > tmax + 1e-6:
                break

            # Advance one time step
            it += 1
            t += dt
            self.simulation.data['dt'].assign(dt)
            self.simulation.hooks.new_timestep(it, t, dt)

            # Update the velocity and pressure fields
            self.apply_analytical_solution()

            # Move u -> up, up -> upp and prepare for the next time step
            for d in range(self.simulation.ndim):
                u_new = self.simulation.data['u%d' % d]
                up = self.simulation.data['up%d' % d]
                upp = self.simulation.data['upp%d' % d]
                upp.assign(up)
                up.assign(u_new)

            # Postprocess this time step
            sim.hooks.end_timestep()

        # We are done
        sim.hooks.simulation_ended(success=True)
