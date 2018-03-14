from math import pi, tanh, sqrt
import dolfin
from ocellaris.utils import ocellaris_error, OcellarisCppExpression, verify_key
from . import register_known_field, KnownField, DEFAULT_POLYDEG


@register_known_field('AiryWaves')
class AiryWaveField(KnownField):
    description = 'Linear airy waves'
    
    def __init__(self, simulation, field_inp):
        """
        A linear Airy wave field (sum of linear sine wave components)
        """
        self.simulation = simulation
        self.read_input(field_inp)
        
        # Show the input data
        simulation.log.info('Creating a linear Airy wave field %r' % self.name)
        simulation.log.info('    Wave depth below: %r' % self.h)
        simulation.log.info('    Wave depth above: %r' % self.h_above)
        simulation.log.info('    Pos. free surface: %r' % self.still_water_pos)
        simulation.log.info('    Vertical comp. gravity: %r' % self.g)
        simulation.log.info('    Wave frequencies: %r' % self.omegas)
        simulation.log.info('    Wave periods: %r' % self.periods)
        simulation.log.info('    Wave lengths: %r' % self.wave_lengths)
        simulation.log.info('    Wave numbers: %r' % self.wave_numbers)
        simulation.log.info('    Wave phases: %r' % self.thetas)
        simulation.log.info('    Wave amplitudes: %r' % self.amplitudes)
        simulation.log.info('    Current speed: %r' % self.current_speed)
        simulation.log.info('    Wind speed: %r' % self.wind_speed)
        simulation.log.info('    Polynomial degree: %r' % self.polydeg)
        
        self._cpp = {}
        self.construct_cpp_code()
        self._expressions = {}
        self._functions = {}
        self.V = dolfin.FunctionSpace(simulation.data['mesh'], 'CG', self.polydeg)
        simulation.hooks.add_pre_timestep_hook(self.update, 'Update Airy wave field %r' % self.name)
        self._dependent_fields = []
    
    def read_input(self, field_inp):
        sim = self.simulation
        self.name = field_inp.get_value('name', required_type='string')
        
        # Get global physical constants
        g = abs(sim.data['g'].values()[-1])
        if g == 0:
            ocellaris_error('Airy waves require gravity',
                            'Cannot compute Airy waves when the vertical component of gravity is 0')
        h = field_inp.get_value('depth', required_type='float')
        if h <= 0:
            ocellaris_error('Airy waves require a still water depth',
                            'Cannot compute Airy waves when the still water depth is %r' % h)
        
        # Get user specified wave data (user must specify one and only one of these)
        omegas = field_inp.get_value('omegas', None, required_type='list(float)')
        periods = field_inp.get_value('periods', None, required_type='list(float)')
        wave_lengths = field_inp.get_value('wave_lengths', None, required_type='list(float)')
        wave_numbers = field_inp.get_value('wave_numbers', None, required_type='list(float)')
        
        # Compute the missing data
        self.omegas, self.periods, self.wave_lengths, self.wave_numbers = \
            get_airy_wave_specs(g, h, omegas, periods, wave_lengths, wave_numbers)
        Nwave = len(self.omegas)
        self.stationary = Nwave == 0
        
        self.still_water_pos = field_inp.get_value('still_water_position', required_type='float')
        self.current_speed = field_inp.get_value('current_speed', 0, required_type='float')
        self.wind_speed = field_inp.get_value('wind_speed', 0, required_type='float')
        self.polydeg = field_inp.get_value('polynomial_degree', DEFAULT_POLYDEG, required_type='int')
        self.g = g
        self.h = h
        self.thetas = field_inp.get_value('wave_phases', [0] * Nwave, required_type='list(float)')
        if not len(self.thetas) == Nwave:
            ocellaris_error('Error with wave phase in Airy wave field input',
                            'The length of the wave phase list does not match the number '
                            'of waves specified, %d != %d' % (len(self.thetas), Nwave))
        self.amplitudes = field_inp.get_value('amplitudes', required_type='list(float)')
        if not len(self.amplitudes) == Nwave:
            ocellaris_error('Error with wave amplitudes in Airy wave field input',
                            'The length of the wave amplitude list does not match the number '
                            'of waves specified, %d != %d' % (len(self.amplitudes), Nwave))
        max_ampl = sum(abs(a) for a in self.amplitudes)
        self.h_above = field_inp.get_value('depth_above', 2*max_ampl, required_type='float')
    
    def update(self, timestep_number, t, dt):
        """
        Called by simulation.hooks on the start of each time step
        """
        if self.stationary:
            return
        
        # Update C++ expressions
        for name, func in self._functions.items():
            expr, updater = self._expressions[name]
            updater(timestep_number, t, dt)
            func.interpolate(expr)
        
        # Update dependet
        for f in self._dependent_fields:
            f.update(timestep_number, t, dt)
    
    def register_dependent_field(self, field):
        """
        We may have dependent fields, typically still water outflow
        fields that use our functions in order to compute their own
        functions 
        """
        self._dependent_fields.append(field)
    
    def construct_cpp_code(self):
        """
        This code runs once at setup time and constructs the C++ code that
        defines the analytical Airy wave field below and above the free surface
        
        Above the free surface the velocities are made continuous by reversing
        the Wheeler streatching a distance h_above up (depth_above input param) 
        """
        rho_min, rho_max = self.simulation.multi_phase_model.get_density_range()
        
        for name in 'elevation c rho uhoriz uvert pdyn pstat ptot'.split():
            # C++ code for still water
            # The fact that the first element of below_cpp is the static
            # value is used in the Wheeler stretching of velocities above
            # the free surface. Beware when refactoring this code!
            above_cpp = None
            blend_sign = None
            if name == 'elevation':
                below_cpp = ['swpos']
            elif name == 'c':
                below_cpp = ['1']
                above_cpp = ['0']
            elif name == 'rho':
                below_cpp = ['rho_max']
                above_cpp = ['rho_min']
            elif name == 'uhoriz':
                blend_sign = '-'
                below_cpp = [repr(self.current_speed)]
                above_cpp = [repr(self.wind_speed)]
            elif name == 'uvert':
                blend_sign = '+'
                below_cpp = ['0']
                above_cpp = ['0']
            elif name == 'pdyn':
                below_cpp = ['0']
                above_cpp = ['0']
            elif name in ('pstat', 'ptot'):
                below_cpp = ['-1 * rho_max * g * D']
                above_cpp = ['0']
            
            # Construct C++ code to compute the named variable
            Nwave = len(self.omegas)
            for i in range(Nwave):
                params = dict(a=self.amplitudes[i],
                              w=self.omegas[i],
                              k=self.wave_numbers[i],
                              theta=self.thetas[i])
                cppb = None
                cppa = None
                if name == 'elevation':
                    cppb = '{a} * sin({w} * t - {k} * X + {theta})'.format(**params)
                elif name == 'uhoriz':
                    cppb = '{w} * {a} * cosh({k} * (Zp + h)) / sinh({k} * h) * sin({w} * t - {k} * X + {theta})'.format(**params)
                elif name == 'uvert':
                    cppb = '{w} * {a} * sinh({k} * (Zp + h)) / sinh({k} * h) * cos({w} * t - {k} * X + {theta})'.format(**params)
                elif name in ('pdyn', 'ptot'):
                    cppb = 'rho_max * g * {a} * cosh({k} * (Zp + h)) / cosh({k} * h) * sin({w} * t - {k} * X + {theta})'.format(**params)
                
                if cppb is not None:
                    below_cpp.append(cppb)
                if cppa is not None:
                    above_cpp.append(cppa)
            
            # Wrap the code in a lambda function that is evaluated immediately.
            # (to allow defining helper variables and write cleaner code)
            lamcode = '[&]() {\n  %s;\n}();'
            
            # Lines inside the lambda function. Some will be unused, but
            # the C++ compiler should be able to remove these easily
            lines = ['double val;',
                     'double X = x[0];',
                     'double Z0 = x[%d];' % (self.simulation.ndim - 1),
                     'double swpos = %r;' % self.still_water_pos,
                     'double Z = Z0 - swpos;',
                     'double h = %r;' % self.h,
                     'double h_above = %r;' % self.h_above,
                     'double g = %r;' % self.g,
                     'double rho_min = %r;' % rho_min,
                     'double rho_max = %r;' % rho_max]
            
            if name != 'elevation':
                # Compute the vertical distance D to the free surface
                # (negative below the free surface, positive above) 
                elev_cpp = self._cpp['elevation'].replace('  ', '    ')
                lines.append('double elev = %s;' % elev_cpp)
                lines.append('double D = Z0 - elev;')
                # Wave elevation relative to still water
                lines.append('double eta = elev - swpos;')
                # Wheeler stretching below the free surface
                # Zp is between -h and 0 when Z is between -d and eta
                lines.append('double Zp = h * (Z + h) / (eta + h) - h;')
            
            full_code_below = ' + '.join(below_cpp)
            if above_cpp is None:
                # No special treatment of values above the free surface
                lines.append('val = %s;' % full_code_below)
            elif blend_sign is not None:
                # Separate between values above and below the free surface
                # The potential flow solution is that the velocities have
                # opposite sign above the free surface
                full_code_above = ' + '.join(above_cpp)
                if Nwave == 0:
                    full_code_below_above = '0'
                else:
                    # Wave values without the static component
                    full_code_below_above = ' + '.join(below_cpp[1:])
                lines.append('if (D <= 0) {')
                lines.append('  val = %s;' % full_code_below)
                lines.append('} else if (D <= h_above and h_above > 0) {')
                lines.append('  val = %s;' % full_code_above)
                # Wheeler stretching above the free surface
                # Zp is between -h and 0 when Z is between eta and h_above
                lines.append('  Zp = (Z * h - eta * h) / (eta - h_above);')
                lines.append('  val %s= %s;' % (blend_sign, full_code_below_above))
                lines.append('} else {')
                lines.append('  val = %s;' % full_code_above)
                lines.append('}')
            else:
                # Separate between values above and below the free surface
                full_code_above = ' + '.join(above_cpp)
                lines.append('if (D <= 0) {')
                lines.append('  val = %s;' % full_code_below)
                lines.append('} else {')
                lines.append('  val = %s;' % full_code_above)
                lines.append('}')
            lines.append('return val;')
            self._cpp[name] = lamcode % ('\n  '.join(lines))
    
    def _get_expression(self, name):
        keys = list(self._cpp.keys()) + ['u']
        verify_key('variable', name, keys, 'Airy wave field %r' % self.name)
        if name not in self._expressions:
            expr, updater = OcellarisCppExpression(self.simulation, self._cpp[name],
                                                   'Airy wave field %r' % name,
                                                   self.polydeg, update=False,
                                                   return_updater=True)
            self._expressions[name] = expr, updater
        return self._expressions[name][0]
    
    def get_variable(self, name):
        if name == 'u':
            # Assume uhoriz == u0 and uvert = uD where D is 2 or 3
            if self.simulation.ndim == 2:
                return dolfin.as_vector([self.get_variable('uhoriz'),
                                         self.get_variable('uvert')])
            else:
                return dolfin.as_vector([self.get_variable('uhoriz'),
                                         dolfin.Constant(0.0),
                                         self.get_variable('uvert')])
        
        if name not in self._functions:
            expr = self._get_expression(name)
            self._functions[name] = dolfin.interpolate(expr, self.V)
        return self._functions[name]


def get_airy_wave_specs(g, h, omegas=None, periods=None, wave_lengths=None, wave_numbers=None):
    """
    Give one of omegas, periods, wave_lengths or wave_numbers. Leave
    the others as None. The input must be a list. This function will
    return the lists omegas, periods, wave_lengths or wave_numbers
    which are equal length and correspond to the same waves according
    to linear Airy wave theory
    
    Parameters:
    
    - g is the magnitude of the acceleration of gravity
    - h is the depth, 
    
    Linear wave relations:
    
        period = 2 * pi / omega     (omega is the wave frequency in rad/s)
        length = 2 * pi / k         (k is the wave number in m^-1)
        omega²/g = k * tanh(k * h)  (the linear dispersion relation)
    
    """
    def err_inp(name1, name2):
        ocellaris_error('Airy wave input error',
                        'You have given both %s and %s, please specify only one!'
                        % (name1, name2))
    
    # Check input and make sure either omegas or wave_numers is defined
    inp = dict(omegas=omegas, periods=periods, wave_lengths=wave_lengths, wave_numbers=wave_numbers)
    if omegas is not None:
        for name, val in inp.items():
            if name is not 'omegas' and val is not None:
                err_inp('omegas', name)
    elif periods is not None:
        for name, val in inp.items():
            if name is not 'periods' and val is not None:
                err_inp('periods', name)
        omegas = [2*pi/t for t in periods]
    elif wave_numbers is not None:
        for name, val in inp.items():
            if name is not 'wave_numbers' and val is not None:
                err_inp('wave_numbers', name)
    elif wave_lengths is not None:
        for name, val in inp.items():
            if name is not 'wave_lengths' and val is not None:
                err_inp('wave_lengths', name)
        wave_numbers = [2*pi/lam for lam in wave_lengths]
        
    # Define remaining variables
    if omegas is None:
        omegas = [sqrt(g * k * tanh(k * h)) for k in wave_numbers]
    if periods is None:
        periods = [2*pi/w for w in omegas]
    if wave_numbers is None:
        wave_numbers = [calc_wave_number(g, h, w) for w in omegas]
    if wave_lengths is None:
        wave_lengths = [2*pi/k for k in wave_numbers]
    
    return omegas, periods, wave_lengths, wave_numbers


def calc_wave_number(g, h, omega, relax=0.5, eps=1e-15):
    """
    Relaxed Picard iterations to find k when omega is known
    """
    k0 = omega**2 / g
    for _ in range(100):
        k1 = omega**2 / g / tanh(k0 * h)
        if abs(k1 - k0) < eps:
            break
        k0 = k1 * relax + k0 * (1 - relax)
    else:
        ocellaris_error('calc_wave_number did not converge',
                        'Input g=%r h=%r omega=%r, tolerance=%e'
                        % (g, h, omega, eps))
    return k1
