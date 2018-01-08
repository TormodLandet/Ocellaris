import sys
import os
import numpy
import yaml
from io import StringIO
from .files import get_result_file_type

if sys.version[0] != '2':
    unicode = lambda s: s


class Results(object):
    def __init__(self, file_name, derived=True, inner_iterations=True):
        """
        Represents the results from an Ocellaris simulation
        
        This is supposed to work as a back end for post processing
        use cases that are not covered by Paraview and similar
        programs and is hence intended for lighter data such as
        time step reports etc
        
        The file name given can be either a simulation log file
        (for ongoing simulations) or an Ocellaris restart file 
        """
        # Main info
        self.input = None
        self.log = None
        self.reports = None
        self.reports_x = None
        self.surfaces = None
        self.input = None
        
        # Auxillary info
        self.ndofs = None
        self.ncpus = None
        
        self.reload(file_name, derived, inner_iterations)
    
    def reload(self, file_name=None, derived=True, inner_iterations=True):
        if file_name is None:
            file_name = self.file_name
        
        self.file_name = os.path.abspath(file_name)
        self.file_type = get_result_file_type(file_name)
        if self.file_type == 'h5':
            read_h5_data(self)
        elif self.file_type == 'log':
            read_log_data(self)
        else:
            raise IOError('Unknown result file type of file %r' % file_name)
        
        # Add derived reports
        reps = self.reports
        if derived:
            if 'Ep' in reps and 'Ek' in reps and 'Et' not in reps:
                reps['Et'] = reps['Ek'] + reps['Ep']
            if 'mass' in reps:
                m, t = reps['mass'], reps['timesteps']
                dm = numpy.zeros_like(m)
                dm[1:] = (m[1:] - m[:-1])/(t[1:] - t[:-1])
                reps['mass change'] = dm
        
        # Set the time to be on the x axis for the report plots
        self.reports_x = {}
        if self.reports:
            for report_name in reps:
                self.reports_x[report_name] = reps['timesteps']
            del self.reports['timesteps']
        
        # Add inner iteration reports from the log
        if inner_iterations:
            read_iteration_reports(self)
        
        if self.surfaces is None:
            read_surfaces(self)
        else:
            for surf in self.surfaces.values():
                surf.reload()
    
    def get_file_path(self, name, check=True):
        """
        Try to get the path of an output file based on
        the information in the input file and the location
        of the restart/log file
        """
        prefix = self.input.get('output', {}).get('prefix', '')
        prefix = self._process_inp(prefix)
        fn = prefix + name
        
        loc = self.file_name.rfind(prefix)
        pth = './' + fn
        if prefix and loc != -1:
            pth = self.file_name[:loc] + fn
        
        if check:
            if not os.path.exists(pth):
                bd = os.path.split(self.file_name)[0]
                pth = os.path.join(bd, fn)
            
            if not os.path.exists(pth):
                raise IOError('Could not find file %r when prefix is %r and result file is %r'
                              % (name, prefix, self.file_name))
        
        return pth
    
    def _process_inp(self, value):
        """
        Simplified handling of Python coded input fields
        """
        if not isinstance(value, str) or not 'py$' in value:
            return value
        value = value.strip()
        assert value.startswith('py$')
        return self._eval(value[3:])
    
    def _eval(self, code):
        consts = self.input.get('user_code', {}).get('constants')
        return eval(code, globals(), consts)


class IsoSurfaces(object):
    def __init__(self, name, field_name, value, file_name):
        self.name = name
        self.field_name = field_name
        self.value = value
        self.file_name = file_name
        self._cache = None
        
    def reload(self):
        self._cache = None
    
    def get_surfaces(self, cache=True):
        if cache and self._cache is not None:
            return self._cache
        
        timesteps = []
        data = []
        
        with open(self.file_name, 'rt') as f:
            description = f.readline()[1:].strip()
            value = float(f.readline().split()[-1])
            dim = int(f.readline().split()[-1])
            
            line = f.readline()
            while line:
                wds = line.split()
                try:
                    time = float(wds[1])
                    nsurf = int(wds[3])
                except Exception:
                    break
                
                if nsurf == 0:
                    timesteps.append(time)
                    data.append([])
                    line = f.readline()
                    continue
                
                datalines = [f.readline() for _ in range(nsurf*3)]
                if not datalines[-1]:
                    break
                timesteps.append(time)
                data.append([])
                for i in range(nsurf):
                    xvals = [float(v) for v in datalines[i*3+0].split()]
                    yvals = [float(v) for v in datalines[i*3+1].split()]
                    zvals = [float(v) for v in datalines[i*3+2].split()]
                    data[-1].append((xvals, yvals, zvals))
                    
                line = f.readline()
        
        res = (description, value, dim, numpy.array(timesteps), data)
        if cache:
            self._cache = res
        return res


def read_h5_data(results):
    """
    Read metadata and reports from a restart file (HDF5 format)
    """
    import h5py
    hdf = h5py.File(results.file_name, 'r')
    
    string_datasets = 'input_file' in hdf['/ocellaris']
    
    # Read the input file
    if string_datasets:
        input_str = hdf['/ocellaris/input_file'].value
    else:
        input_str = hdf['/ocellaris'].attrs['input_file']
    results.input = yaml.load(input_str)
    
    # Read reports
    reps = {}
    N = 1e100
    for rep_name in hdf['/reports']:
        arr = numpy.array(hdf['/reports'][rep_name])
        reps[rep_name] = arr
        N = min(N, len(arr))
    
    # Ensure equal length arrays
    for key in list(reps.keys()):
        reps[key] = reps[key][:N]
    
    # Read log
    if string_datasets:
        log = hdf['/ocellaris/full_log'].value
    else:
        log = []
        i = 0
        while True:
            logname = 'full_log_%d' % i
            i += 1
            if not logname in hdf['/ocellaris'].attrs:
                break
            log.append(hdf['/ocellaris'].attrs[logname])
        log = ''.join(log)
    
    results.reports = reps
    results.log = log


def read_log_data(results):
    """
    Read metadata and reports from a log file (ASCII format)
    """
    INP_START = '----------------------------- configuration begin -'
    INP_END = '------------------------------ configuration end -'
    in_input_section = False
    input_strs = []
    data = {}
    
    # Read input and timestep reports from log file
    with open(results.file_name, 'rt') as f:
        for line in f:
            if line.startswith(INP_START):
                in_input_section = True
            elif line.startswith(INP_END):
                in_input_section = False
            elif in_input_section:
                input_strs.append(line)
            elif line.startswith('Reports for timestep'):
                parts = line[12:].split(',')
                for pair in parts:
                    try:
                        key, value = pair.split('=')
                        key = key.strip()
                        value = float(value)
                        data.setdefault(key, []).append(value)
                    except:
                        break
        f.seek(0)
        log = f.read()
    if data:
        del data['timestep']
    
    # Read the input section
    if input_strs:
        input_str = ''.join(input_strs)
        results.input = yaml.load(input_str)
    else:
        results.input = {}
    
    reps = {}
    N = 1e100
    for key, values in data.items():
        arr = numpy.array(values)
        if key == 'time':
            key = 'timesteps'
        reps[key] = arr
        N = min(N, len(arr))
    
    # Ensure equal length arrays in case of partially written 
    # time steps on the log file
    for key in list(reps.keys()):
        reps[key] = reps[key][:N]
    
    results.reports = reps
    results.log = log


def read_iteration_reports(results):
    """
    Read less inportant reports that are on the log file, but not
    saved specifically to the h5 file format. This is mainly inner
    iteration data, but also some other tibits like number of
    degrees of freedom, number of CPUs, ...
    """
    iter_reps = {}
    for line in StringIO(unicode(results.log)):
        if line.startswith('Running simulation on'):
            results.ncpus = int(line.split()[3])
        elif line.startswith('Degrees of freedom'):
            results.ndofs = int(line.split()[3])
        elif not ('iteration' in line and 'Krylov' in line):
            continue
        
        # Parse the inner iteration line
        try:
            parts = line.split(' - ')
            iter_num = int(parts[0].split()[-1])
            iter_reps.setdefault('iteration', []).append(iter_num)
            for part in parts[1:]:
                if 'Krylov' in part:
                    continue
                wds = part.split()
                if wds[0] in ('u', 'p') and len(wds) == 2:
                    value = int(wds[1])
                    name = 'solver iterations %s' % wds[0]
                    iter_reps.setdefault(name, []).append(value)
                else:
                    name = ' '.join(wds[:-1])
                    value = float(wds[-1])
                    iter_reps.setdefault(name, []).append(value)
        
        except:
            pass
    
    if not iter_reps:
        return
    
    # Get minimum length
    Nmin = 1e100
    for name, value in iter_reps.items():
        Nmin = min(Nmin, len(value))
    
    # Store reports and crop reports to same length
    xaxis = numpy.arange(Nmin)
    for name, value in iter_reps.items():
        name2 = 'Inner iteration: %s' % name
        results.reports[name2] = numpy.array(value[:Nmin])
        results.reports_x[name2] = xaxis


def read_surfaces(res):
    inp = res.input
    res.surfaces = {}
    if not 'probes' in inp:
        return
    
    for probe in inp['probes']:
        if not (probe.get('enabled', True) and
                probe.get('type', '') == 'IsoSurface'):
            continue
        name = probe['name']
        field_name = probe['field']
        value = probe['value']
        file_name_postfix = probe['file_name']
        file_name = res.get_file_path(file_name_postfix)
        isosurf = IsoSurfaces(name, field_name, value, file_name) 
        res.surfaces[name] = isosurf
