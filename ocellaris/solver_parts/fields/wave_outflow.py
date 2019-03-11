# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from ocellaris.utils import ocellaris_error, OcellarisCppExpression, verify_key
from . import register_known_field, KnownField


OPTIMISE_FLUXES = True


@register_known_field('WaveOutflow')
class WaveOutflowField(KnownField):
    description = 'Outflow from a wave tank'

    def __init__(self, simulation, field_inp):
        """
        This field is used to match the outflow in a wave tank to the
        inflow BCs in such a way that the mass flux of the two phases
        are perfectly balanced and no accumulation of mass will occur.

        The implementation assumes that inflow and outflow parts of
        the domain boundary are perpendicular to the wave traveling
        direction (x[0] direction) and that the outwards normal on the
        inlet and outlet are pointing in opposite directions
        """
        self.simulation = simulation
        self.read_input(field_inp)

        # Show the input data
        simulation.log.info('Creating a wave outflow field %r' % self.name)
        simulation.log.info('    Inflow BC field: %s' % self.inflow_field_name)

        # Get input from the incoming wave field
        self.V = self.inflow_field.V
        self.polydeg = self.inflow_field.polydeg
        self.stationary = self.inflow_field.stationary
        self.still_water_pos = self.inflow_field.still_water_pos

        # Our tunable parameters
        self.const_speed_above = dolfin.Constant(1.0)
        self.const_speed_below = dolfin.Constant(1.0)

        # Make the C++ code that computes the outflow
        self._cpp = {}
        self.construct_cpp_code()
        self._expressions = {}
        self._functions = {}

        # Tune factors used to make the outflow match the inflow perfectly
        # in terms of the volume fractions transported through the domain
        self.tune_factors(True)

        # Make the field call our update function when it itself has updated
        self.inflow_field.register_dependent_field(self)

    def read_input(self, field_inp):
        sim = self.simulation
        self.name = field_inp.get_value('name', required_type='string')
        self.inflow_region_name = field_inp.get_value(
            'inflow_region', required_type='string'
        )
        self.outflow_region_name = field_inp.get_value(
            'outflow_region', required_type='string'
        )

        def err(msg):
            return ocellaris_error(
                'Error in definition of WaveOutflow field %r' % self.name, msg
            )

        # Boundary conditions have not been constructed yet, so we need to read the input file
        # to get the correct BC input
        bcs = sim.input.get_value('boundary_conditions', required_type='list(dict)')

        # Get inflow boundary region index
        for i, bc_region in enumerate(bcs):
            if bc_region.get('name', None) == self.inflow_region_name:
                inflow_index = i
                self.inflow_ds_mark_id = inflow_index + 1
                break
        else:
            err('Did not find inflow boundary region %r' % self.inflow_bc_name)

        # Get outflow boundary region index
        for i, bc_region in enumerate(bcs):
            if bc_region.get('name', None) == self.outflow_region_name:
                outflow_index = i
                self.outflow_ds_mark_id = outflow_index + 1
                break
        else:
            err('Did not find outflow boundary region %r' % self.outflow_region_name)

        # Get the velocity BCs in the inflow region
        bc_inp = sim.input.get_value(
            'boundary_conditions/%d' % inflow_index, required_type='Input'
        )
        if not 'u' in bc_inp:
            err('Did not find velocity BC in region %r' % self.inflow_bc_name)

        # Get the field function used for the velocity BCs in the inflow region
        vel_inp = sim.input.get_value(
            'boundary_conditions/%d/u' % inflow_index, required_type='Input'
        )
        bc_type = vel_inp.get_value('type', '', 'string')
        if not bc_type == 'FieldFunction':
            err(
                'Velocity BC in region %r is not a FieldFunction!' % self.inflow_bc_name
            )
        func_vardef = vel_inp.get_value('function', required_type='string')

        # Get the name of the incomming wave field
        split_vardef = func_vardef.strip().split('/')
        if not len(split_vardef) == 2:
            err(
                'Velocity BC function in region %r is not a valid field function specifier!'
                % self.inflow_bc_name
            )
        self.inflow_field_name = split_vardef[0]

        # Get the incomming wave field object
        if self.inflow_field_name not in sim.fields:
            err(
                (
                    'Inflow field %s is not (yet) defined. It must be defined before this '
                    'WaveOutflow field'
                )
                % self.inflow_field_name
            )
        self.inflow_field = sim.fields[self.inflow_field_name]

    def construct_cpp_code(self):
        """
        This code runs once at init time and generates the C++ code used to
        define the outflow field
        """
        # Wrap the code in a lambda function that is evaluated immediately.
        # (to allow defining helper variables and write cleaner code)
        lamcode = '[&]() {\n  %s;\n}();'

        for name in 'uhoriz c'.split():
            lines = [
                'double val;',
                'double X = x[0];',
                'double Z0 = x[%d];' % (self.simulation.ndim - 1),
                'double swpos = %r;' % self.still_water_pos,
            ]

            if name == 'uhoriz':
                lines.append('if (Z0 <= swpos) {')
                lines.append('  val = u_below;')
                lines.append('} else {')
                lines.append('  val = u_above;')
                lines.append('}')
            elif name == 'c':
                lines.append('if (Z0 <= swpos) {')
                lines.append('  val = 1.0;')
                lines.append('} else {')
                lines.append('  val = 0.0;')
                lines.append('}')

            lines.append('return val;')
            self._cpp[name] = lamcode % ('\n  '.join(lines))

    def tune_factors(self, setup=False):
        if setup:
            # Incoming wave data
            u_inlet = self.inflow_field.get_variable('uhoriz')
            c_inlet = self.inflow_field.get_variable('c')
            ds_inlet = self.simulation.data['ds'](self.inflow_ds_mark_id)

            # Outlet data
            u_outlet = self._get_expression('uhoriz')
            c_outlet = self._get_expression('c')
            ds_outlet = self.simulation.data['ds'](self.outflow_ds_mark_id)

            # Compute unit fluxes
            self.const_speed_above.assign(dolfin.Constant(1.0))
            self.const_speed_below.assign(dolfin.Constant(1.0))
            self.uf_above = dolfin.assemble((1 - c_outlet) * u_outlet * ds_outlet)
            self.uf_below = dolfin.assemble(c_outlet * u_outlet * ds_outlet)

            self.form_inlet_above = dolfin.Form((1 - c_inlet) * u_inlet * ds_inlet)
            self.form_inlet_below = dolfin.Form(c_inlet * u_inlet * ds_inlet)
            self.form_inlet_total = dolfin.Form(u_inlet * ds_inlet)

            self.form_outlet_above = dolfin.Form((1 - c_outlet) * u_outlet * ds_outlet)
            self.form_outlet_below = dolfin.Form(c_outlet * u_outlet * ds_outlet)
            self.form_outlet_total = dolfin.Form(u_outlet * ds_outlet)

        import time

        t1 = time.time()

        # Compute flux above and below at the inlet
        inlet_flux_above = dolfin.assemble(self.form_inlet_above)
        inlet_flux_below = dolfin.assemble(self.form_inlet_below)
        inlet_flux_total = dolfin.assemble(self.form_inlet_total)
        is_motion = abs(inlet_flux_above) > 0 or abs(inlet_flux_below) > 0

        # Estimate the fluxes at the outlet (this gets within 99.99% of total zero flux)
        self.const_speed_above.assign(dolfin.Constant(inlet_flux_above / self.uf_above))
        self.const_speed_below.assign(dolfin.Constant(inlet_flux_below / self.uf_below))

        if OPTIMISE_FLUXES and is_motion:
            # Use root finding to adjust the speed above the free surface
            # in order to get zero total volume flux in the domain

            def func_to_minimise(vel_above):
                "Compute the difference in total flux"
                # Adjust outflow speed above the FS
                self.const_speed_above.assign(dolfin.Constant(vel_above))
                outlet_flux_total = dolfin.assemble(self.form_outlet_total)
                return (outlet_flux_total - inlet_flux_total) ** 2

            # Starting values - use the basic unit fluxes
            x1 = inlet_flux_above / self.uf_above
            x0 = x1 * 0.99

            xN, niter = find_root_secant(x0, x1, func_to_minimise)

        if False and abs(inlet_flux_below) > 0:
            outlet_flux_above = dolfin.assemble(self.form_outlet_above)
            outlet_flux_below = dolfin.assemble(self.form_outlet_below)
            outlet_flux_total = dolfin.assemble(self.form_outlet_total)
            print('uf ab & be', self.uf_above, self.uf_below)
            print(
                'flux_above',
                inlet_flux_above,
                outlet_flux_above,
                outlet_flux_above / inlet_flux_above,
            )
            print(
                'flux_below',
                inlet_flux_below,
                outlet_flux_below,
                outlet_flux_below / inlet_flux_below,
            )
            print(
                'flux_total',
                inlet_flux_total,
                outlet_flux_total,
                outlet_flux_total / inlet_flux_total,
            )
            print(xN - x0)
            print(niter)
            print('Took %g seconds' % (time.time() - t1))

    def update(self, timestep_number, t, dt):
        """
        Called by the incomming wave field on the start of each time step
        """
        if self.stationary:
            return

        # Compute new outlet fluxes
        with dolfin.Timer('Ocellaris tune wave outflow factors'):
            self.tune_factors()

        # Update C++ expressions
        for name, func in self._functions.items():
            if name == 'c':
                # The c field is stationary
                continue
            expr, updater = self._expressions[name]
            updater(timestep_number, t, dt)
            func.interpolate(expr)

    def _get_expression(self, name):
        keys = list(self._cpp.keys()) + ['u', 'uvert']
        verify_key('variable', name, keys, 'Wave outflow field %r' % self.name)

        params = dict(u_above=self.const_speed_above, u_below=self.const_speed_below)

        if name not in self._expressions:
            expr, updater = OcellarisCppExpression(
                self.simulation,
                self._cpp[name],
                'Wave outflow field %r' % name,
                self.polydeg,
                update=False,
                return_updater=True,
                params=params,
            )
            self._expressions[name] = expr, updater
        return self._expressions[name][0]

    def get_variable(self, name):
        zero = dolfin.Constant(0.0)
        if name == 'u':
            # Assume that the waves are traveling in x-direction
            if self.simulation.ndim == 2:
                return dolfin.as_vector([self.get_variable('uhoriz'), zero])
            else:
                return dolfin.as_vector([self.get_variable('uhoriz'), zero, zero])
        elif name == 'uvert':
            return zero

        if name not in self._functions:
            expr = self._get_expression(name)
            self._functions[name] = dolfin.interpolate(expr, self.V)
        return self._functions[name]


def find_root_secant(x0, x1, func, niter=100, tolerance=1e-14):
    """
    Root finding using the secant method
    https://en.wikipedia.org/wiki/Secant_method

    Iterates niter times and returns the last value of x or
    the value of x that makes func less than the tolerance
    if this is found first. Also returns the number of
    iterations performed.

    @return: (x_root, num_iters_performed)
    """
    f0 = func(x0)
    f1 = func(x1)

    for i in range(niter):
        x_new = (x0 * f1 - x1 * f0) / (f1 - f0)
        f_new = func(x_new)
        if abs(f_new) < tolerance:
            break
        x0, f0 = x1, f1
        x1, f1 = x_new, f_new
    return x_new, i + 1
