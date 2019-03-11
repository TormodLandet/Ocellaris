# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

from .vof_dg0_surface_locator import FreeSurfaceLocatorImplDG0


def get_free_surface_locator(simulation, name, function, value):
    """
    Get a free surface locator. Typical input:

        name = 'c'
        function = <dolfin.Function DG0 colour function>
        value = 0.5

    The locator can in principle be used for non-density fields to locate
    other iso-surfaces than the free surface. The locator is cached in
    simulation.iso_surface_locators to be able to reuse the locator in
    multiple parts of Ocellaris without these parts being aware of each
    other.
    """
    key = (name, value)
    cache = simulation.iso_surface_locators
    if key not in cache:
        cache[key] = FreeSurfaceLocator(simulation, function, name, value)
    return cache[key]


class FreeSurfaceLocator:
    def __init__(self, simulation, density, density_name, free_surface_value):
        """
        Function to find the location of the free surface from a density
        (or VOF) function that is probably DG, and a free surface value
        (which is 0.5 for VOF)

        Defines an attribute ``crossing_points`` that is a dictionary
        mapping cell index to a list of crossing points in that cell.
        Cells without any crossing points will not be contained in the
        dictionary.
        """
        self.simulation = simulation
        self.density = density
        self.density_name = density_name
        self.free_surface_value = free_surface_value
        self.installed_hooks = set()
        self.callbacks = []

        V = density.function_space()
        deg = V.ufl_element().degree()
        fam = V.ufl_element().family()
        if fam == 'Discontinuous Lagrange' and deg == 0:
            self.impl = FreeSurfaceLocatorImplDG0(simulation, density, free_surface_value)
        else:
            raise NotImplementedError(
                'FreeSurfaceLocator not implemented for %r degree %r' % (fam, deg)
            )

        # Set need_update to True to prepare for queries
        self.update()

    @property
    def crossing_points(self):
        """
        A dictionary mapping cell index to list of crossing point 
        coordinates (tuples)
        """
        if self._needs_update:
            # Update the position of the free surface based on an updated
            # density field. This can be an expensive operation
            self._crossing_points = self.impl.compute_crossing_points()
        return self._crossing_points

    def add_update_hook(self, hook_name, callback=None, description=None):
        """
        Add a hook point to connect to in order to be notified whenever
        the underlying field changes. In this way the user of the free
        surface locator does not have to wory about calling .update() and
        it hence avoids calling update multiple times when multiple parts
        of Ocellaris use the same free surface locator
        """
        # Avoid installing multiple hooks, wo only need one notification
        if hook_name in self.installed_hooks:
            return

        def hook(*argv, **kwargs):
            return self.update()

        self.simulation.hooks.add_custom_hook(
            hook_name,
            hook,
            'Update FreeSurfaceLocator for %s at %r' % (self.density_name, self.free_surface_value),
        )
        self.installed_hooks.add(hook_name)

        if callback is not None:
            self.callbacks.append((callback, description))

    def update(self):
        self._needs_update = True
        self._crossing_points = None
        for cb, descr in self.callbacks:
            if descr is None:
                descr = 'FreeSurfaceLocator callback'
            self.simulation.hooks.call_after(cb, descr)
