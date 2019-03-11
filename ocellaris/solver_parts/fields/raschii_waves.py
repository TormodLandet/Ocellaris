# Copyright (C) 2018-2019 Tormod Landet
# SPDX-License-Identifier: Apache-2.0

import time
from raschii import get_wave_model
from ocellaris.utils import ocellaris_error
from . import register_known_field, DEFAULT_POLYDEG
from .base_wave_field import BaseWaveField, COLOUR_PROJECTION_DEGREE


@register_known_field('RaschiiWaves')
class RaschiiWaveField(BaseWaveField):
    description = 'Nonlinear regular waves'

    def __init__(self, simulation, field_inp):
        """
        Nonlinear regular waves using the Raschii Python library
        """
        super().__init__(simulation, field_inp)
        self.read_input(field_inp)
        simulation.log.info('Creating nonlinear regular wave field %r' % self.name)
        simulation.log.info('    Wave depth below: %r' % self.h)
        simulation.log.info('    Wave depth above: %r' % self.h_above)
        simulation.log.info('    Blending height: %r' % self.blending_height)
        simulation.log.info('    Pos. free surface: %r' % self.still_water_pos)
        simulation.log.info('    Vertical comp. gravity: %r' % self.g)
        simulation.log.info('    Wave length: %r' % self.wave_length)
        simulation.log.info('    Wave height: %r' % self.wave_height)
        simulation.log.info('    Wave approximation order: %r' % self.order)
        simulation.log.info('    Polynomial degree: %r' % self.polydeg)
        simulation.log.info(
            '    Colour proj. degree: %r' % self.colour_projection_degree
        )
        simulation.log.info('    Raschii wave model: %s' % self.wave_model)
        simulation.log.info('    Raschii air model: %s' % self.air_model)
        self.construct_cpp_code()

    def read_input(self, field_inp):
        sim = self.simulation
        self.wave_model = field_inp.get_value('wave_model', 'Fenton', 'string')
        self.air_model = field_inp.get_value('air_model', 'FentonAir', 'string')
        self.order = field_inp.get_value('model_order', 5, 'int')

        # Get global physical constants
        g = abs(sim.data['g'].values()[-1])
        if g == 0:
            ocellaris_error(
                'Waves require gravity', 'The vertical component of gravity is 0'
            )
        h = field_inp.get_value('depth', required_type='float')
        if h <= 0:
            ocellaris_error(
                'Waves require a still water depth', 'The still water depth is %r' % h
            )
        self.wave_length = field_inp.get_value('wave_length', required_type='float')
        self.wave_height = field_inp.get_value('wave_height', required_type='float')
        self.stationary = self.wave_height == 0
        self.ramp_time = field_inp.get_value('ramp_time', 0, required_type='float')

        self.still_water_pos = field_inp.get_value(
            'still_water_position', required_type='float'
        )
        self.current_speed = field_inp.get_value(
            'current_speed', 0, required_type='float'
        )
        self.wind_speed = field_inp.get_value('wind_speed', 0, required_type='float')
        self.polydeg = field_inp.get_value(
            'polynomial_degree', DEFAULT_POLYDEG, required_type='int'
        )
        self.g = g
        self.h = h
        h_above = 3 * self.wave_height
        self.h_above = field_inp.get_value(
            'depth_above', h_above, required_type='float'
        )
        self.blending_height = field_inp.get_value('blending_height', None, 'float')

        # Project the colour function to DG0 (set degree to -1 to prevent this)
        self.colour_projection_degree = field_inp.get_value(
            'colour_projection_degree', COLOUR_PROJECTION_DEGREE, 'int'
        )
        self.colour_projection_form = None

    def construct_cpp_code(self):
        """
        This code runs once at setup time and constructs the C++ code that
        defines the analytical Airy wave field below and above the free surface

        Above the free surface the velocities are made continuous by reversing
        the Wheeler streatching a distance h_above up (depth_above input param)
        """
        # Get the requested classes and set up the required input
        WaveClass, AirClass = get_wave_model(self.wave_model, self.air_model)
        wave_args = {
            'height': self.wave_height,
            'depth': self.h,
            'length': self.wave_length,
        }
        if AirClass is not None:
            wave_args['air'] = AirClass(self.h_above, self.blending_height)
        if 'N' in WaveClass.required_input:
            wave_args['N'] = self.order

        # Initialize the Raschii wave model. This may take som time for stream
        # function waves where an optimalization loop needs to run as a part
        # of the wave initialization process
        sim = self.simulation
        sim.log.info('    Creating Raschii wave ...')
        t0 = time.time()
        self.raschii_wave = WaveClass(**wave_args)
        sim.log.info('    Wave created in %.2f seconds' % (time.time() - t0))
        k = self.raschii_wave.k
        c = self.raschii_wave.c
        sim.log.info('    Wave number: %r' % k)
        sim.log.info('    Phase speed: %r' % c)
        sim.log.info('    Wave period: %r' % (self.wave_length / c))

        # Construct the C++ code
        cpp_e = self.raschii_wave.elevation_cpp()
        cpp_u, cpp_w = self.raschii_wave.velocity_cpp(all_points_wet=False)
        self._cpp['elevation'] = cpp_e
        self._cpp['uhoriz'] = cpp_u
        self._cpp['uvert'] = cpp_w
        self._cpp['c'] = 'x[2] <= (%s) ? 1.0 : 0.0' % cpp_e

        # Adjust the z-coordinate such that the bottom is at z=0
        zdiff = self.still_water_pos - self.h
        for k, v in list(self._cpp.items()):
            self._cpp[k] = v.replace('x[2]', '(x[2] - %r)' % zdiff)

        # Adjust the C++ code z-coordinate if the simulation is 2D (then z -> y)
        if self.simulation.ndim == 2:
            for k, v in list(self._cpp.items()):
                self._cpp[k] = v.replace('x[2]', 'x[1]')
