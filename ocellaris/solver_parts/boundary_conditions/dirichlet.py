# Copyright (C) 2015-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import dolfin
from . import register_boundary_condition, BoundaryConditionCreator
from ocellaris.utils import (
    CodedExpression,
    OcellarisCppExpression,
    OcellarisError,
    verify_field_variable_definition,
)


class OcellarisDirichletBC(dolfin.DirichletBC):
    def __init__(
        self, simulation, V, value, subdomain_marker, subdomain_id, updater=None
    ):
        """
        A simple storage class for Dirichlet boundary conditions
        """
        super().__init__(
            V, value, subdomain_marker, subdomain_id, method='geometric'
        )
        self.simulation = simulation
        self._value = value
        self.subdomain_marker = subdomain_marker
        self.subdomain_id = subdomain_id
        self._updater = updater

    def func(self):
        """
        The boundary value derivative function
        """
        return self._value

    def ds(self):
        """
        Returns the ds measure of the subdomain
        """
        return self.simulation.data['ds'](self.subdomain_id)

    def copy_and_change_function_space(self, V):
        """
        Return a copy with a new function space. Used when converting from
        BCs for a segregated solver (default) to BCs for a coupled solver
        """
        return OcellarisDirichletBC(
            self.simulation, V, self._value, self.subdomain_marker, self.subdomain_id
        )

    def update(self):
        """
        Update the time and other parameters used in the BC.
        This is used every timestep and for all RK substeps
        """
        if self._updater:
            self._updater(
                self.simulation.timestep, self.simulation.time, self.simulation.dt
            )

    def __repr__(self):
        return '<OcellarisDirichletBC on subdomain %d>' % self.subdomain_id


@register_boundary_condition('ConstantValue')
class ConstantDirichletBoundary(BoundaryConditionCreator):
    description = 'A prescribed constant value Dirichlet condition'

    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Dirichlet condition with constant value
        """
        self.simulation = simulation
        if var_name[-1].isdigit():
            # A var_name like "u0" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name[:-1]]
        else:
            # A var_name like "u" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name]

        value = inp_dict.get_value('value', required_type='any')
        if isinstance(value, list):
            assert len(value) == simulation.ndim
            for d in range(simulation.ndim):
                name = '%s%d' % (var_name, d)
                self.register_dirichlet_condition(
                    name, value[d], subdomains, subdomain_id
                )
        else:
            self.register_dirichlet_condition(var_name, value, subdomains, subdomain_id)

    def register_dirichlet_condition(self, var_name, value, subdomains, subdomain_id):
        """
        Add a Dirichlet condition to this variable
        """
        if not isinstance(value, (float, int)):
            raise OcellarisError(
                'Error in ConstantValue BC for %s' % var_name,
                'The value %r is not a number' % value,
            )
        df_value = dolfin.Constant(value)

        # Store the boundary condition for use in the solver
        bc = OcellarisDirichletBC(
            self.simulation, self.func_space, df_value, subdomains, subdomain_id
        )
        bcs = self.simulation.data['dirichlet_bcs']
        bcs.setdefault(var_name, []).append(bc)

        self.simulation.log.info('    Constant value %r for %s' % (value, var_name))


@register_boundary_condition('CodedValue')
class CodedDirichletBoundary(BoundaryConditionCreator):
    description = 'A coded Dirichlet condition'

    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Dirichlet condition with coded value
        """
        self.simulation = simulation
        if var_name[-1].isdigit():
            # A var_name like "u0" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name[:-1]]
        else:
            # A var_name like "u" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name]

        # Make a dolfin Expression object that runs the code string
        code = inp_dict.get_value('code', required_type='any')

        if isinstance(code, list):
            assert len(code) == simulation.ndim
            for d in range(simulation.ndim):
                name = '%s%d' % (var_name, d)
                description = 'coded value boundary condition for %s' % name
                sub_code = inp_dict.get_value('code/%d' % d, required_type='string')
                expr = CodedExpression(simulation, sub_code, description)
                self.register_dirichlet_condition(name, expr, subdomains, subdomain_id)
        else:
            description = 'coded value boundary condition for %s' % var_name
            expr = CodedExpression(simulation, code, description)
            self.register_dirichlet_condition(var_name, expr, subdomains, subdomain_id)

    def register_dirichlet_condition(self, var_name, expr, subdomains, subdomain_id):
        """
        Store the boundary condition for use in the solver
        """
        bc = OcellarisDirichletBC(
            self.simulation, self.func_space, expr, subdomains, subdomain_id
        )
        bcs = self.simulation.data['dirichlet_bcs']
        bcs.setdefault(var_name, []).append(bc)
        self.simulation.log.info('    Coded value for %s' % var_name)


@register_boundary_condition('CppCodedValue')
class CppCodedDirichletBoundary(BoundaryConditionCreator):
    description = 'A C++ coded Dirichlet condition'

    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Dirichlet condition with C++ coded value
        """
        self.simulation = simulation
        if var_name[-1].isdigit():
            # A var_name like "u0" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name[:-1]]
        else:
            # A var_name like "u" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name]

        # Make a dolfin Expression object that runs the code string
        code = inp_dict.get_value('cpp_code', required_type='any')

        if isinstance(code, list):
            assert len(code) == simulation.ndim
            for d in range(simulation.ndim):
                name = '%s%d' % (var_name, d)
                sub_code = inp_dict.get_value('cpp_code/%d' % d, required_type='string')
                self.register_dirichlet_condition(
                    name, sub_code, subdomains, subdomain_id
                )
        else:
            self.register_dirichlet_condition(var_name, code, subdomains, subdomain_id)

    def register_dirichlet_condition(
        self, var_name, cpp_code, subdomains, subdomain_id
    ):
        """
        Store the boundary condition for use in the solver
        """
        description = 'boundary condititon for %s' % var_name
        P = self.func_space.ufl_element().degree()
        expr, updater = OcellarisCppExpression(
            self.simulation, cpp_code, description, P, return_updater=True
        )

        bc = OcellarisDirichletBC(
            self.simulation,
            self.func_space,
            expr,
            subdomains,
            subdomain_id,
            updater=updater,
        )
        bcs = self.simulation.data['dirichlet_bcs']
        bcs.setdefault(var_name, []).append(bc)
        self.simulation.log.info('    C++ coded value for %s' % var_name)


@register_boundary_condition('FieldFunction')
class FieldFunctionDirichletBoundary(BoundaryConditionCreator):
    description = 'A Dirichlet condition with values from a field function'

    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Dirichlet boundary condition with value from a field function
        """
        self.simulation = simulation
        if var_name[-1].isdigit():
            # A var_name like "u0" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name[:-1]]
        else:
            # A var_name like "u" was given. Look up "Vu"
            self.func_space = simulation.data['V%s' % var_name]

        # Get the field function expression object
        vardef = inp_dict.get_value('function', required_type='any')
        description = 'boundary condititon for %s' % var_name
        if isinstance(vardef, list):
            assert len(vardef) == simulation.ndim
            exprs = [
                verify_field_variable_definition(simulation, vd, description)
                for vd in vardef
            ]
        else:
            expr = verify_field_variable_definition(simulation, vardef, description)
            if expr.ufl_shape != ():
                assert expr.ufl_shape == (
                    simulation.ndim,
                ), 'Expected shape %r got %r' % ((simulation.ndim,), expr.ufl_shape)
                exprs = [expr[d] for d in range(simulation.ndim)]
            else:
                exprs = [expr]

        # Register BCs
        if len(exprs) > 1:
            for d in range(simulation.ndim):
                name = '%s%d' % (var_name, d)
                self.register_dirichlet_condition(
                    name, exprs[d], subdomains, subdomain_id
                )
        else:
            self.register_dirichlet_condition(
                var_name, exprs[0], subdomains, subdomain_id
            )

    def register_dirichlet_condition(self, var_name, expr, subdomains, subdomain_id):
        """
        Store the boundary condition for use in the solver
        """
        assert expr.ufl_shape == ()
        bc = OcellarisDirichletBC(
            self.simulation, self.func_space, expr, subdomains, subdomain_id
        )
        bcs = self.simulation.data['dirichlet_bcs']
        bcs.setdefault(var_name, []).append(bc)
        self.simulation.log.info('    Field function value for %s' % var_name)


@register_boundary_condition('FieldVelocityValve')
class FieldVelocityValveDirichletBoundary(BoundaryConditionCreator):
    description = 'A Dirichlet condition that compensates for non-zero total flux of a known velocity field'

    def __init__(self, simulation, var_name, inp_dict, subdomains, subdomain_id):
        """
        Dirichlet boundary condition with value from a field function
        """
        self.simulation = simulation

        # A var_name like "u0" should be given. Look up "Vu"
        self.func_space = simulation.data['V%s' % var_name[:-1]]

        # Get the field function expression object
        vardef = inp_dict.get_value('function', required_type='any')
        description = 'boundary condititon for %s' % var_name
        self.velocity = verify_field_variable_definition(
            simulation, vardef, description
        )
        field = simulation.fields[vardef.split('/')[0]]

        # The expression value is updated as the field is changed
        inp_dict.get_value('function', required_type='any')
        field.register_dependent_field(self)
        self.flux = dolfin.Constant(1.0)

        # Create the
        bc = OcellarisDirichletBC(
            self.simulation, self.func_space, self.flux, subdomains, subdomain_id
        )
        bcs = self.simulation.data['dirichlet_bcs']
        bcs.setdefault(var_name, []).append(bc)
        self.simulation.log.info('    Field velocity valve for %s' % var_name)

        # Compute the region area, then update the flux
        mesh = simulation.data['mesh']
        self.area = dolfin.assemble(self.flux * bc.ds()(domain=mesh))
        self.region_names = inp_dict.get_value('regions', required_type='list(string)')
        self.update()

    def update(self, timestep_number=None, t=None, dt=None):
        """
        The main field has changed, update our flux to make the total sum to zero
        """
        regions = self.simulation.data['boundary']
        mesh = self.simulation.data['mesh']
        n = dolfin.FacetNormal(mesh)
        flux = 0
        count = 0
        for region in regions:
            if region.name in self.region_names:
                f = dolfin.dot(self.velocity, n) * region.ds()
                flux += dolfin.assemble(f)
                count += 1
        assert count == len(self.region_names)
        # FIXME: assumes n is pointing outwards along the axis in the positive
        #        direction in this boundary region
        self.flux.assign(dolfin.Constant(-flux / self.area))
