import os
import numpy
import yaml


class Results(object):
    def __init__(self, file_name, derived=True):
        """
        Represents the results from an Ocellaris simulation
        
        This is supposed to work as a back end for post processing
        use cases that are not covered by Paraview and similar
        programs and is hence intended for lighter data such as
        time step reports etc
        
        The file name given can be either a simulation log file
        (for ongoing simulations) or an Ocellaris restart file 
        """
        self.input = None
        self.log = None
        self.reports = None
        self.surfaces = None
        self.input = None
        self.reload(file_name, derived)
    
    def reload(self, file_name=None, derived=True):
        if file_name is None:
            file_name = self.file_name
        
        self.file_name = os.path.abspath(file_name)
        if file_name.endswith('.h5'):
            read_h5_data(self, derived)
        elif file_name.endswith('.log'):
            read_log_data(self, derived)
        
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


def read_h5_data(results, derived):
    import h5py
    hdf = h5py.File(results.file_name, 'r')
    
    # Read the input file
    input_str = hdf['/ocellaris'].attrs['input_file']
    results.input = yaml.load(input_str)
    
    # Read reports
    reps = {}
    for rep_name in hdf['/reports']:
        reps[rep_name] = numpy.array(hdf['/reports'][rep_name])
    
    if derived:
        if 'Ep' in reps and 'Ek' in reps and 'Et' not in reps:
            reps['Et'] = reps['Ek'] + reps['Ep']
            
    # Read log
    log = []
    i = 0
    while True:
        logname = 'full_log_%d' % i
        i += 1
        if not logname in hdf['/ocellaris'].attrs:
            break
        log.append(hdf['/ocellaris'].attrs[logname])
    log = ''.join(log).split('\n')  
    
    results.reports = reps
    results.log = log


def read_log_data(results, derived):
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
        
    # Read the input section
    input_str = ''.join(input_strs)
    results.input = yaml.load(input_str)
    
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
    for key in reps.keys():
        reps[key] = reps[key][:N]
    
    if derived:
        if 'Ep' in reps and 'Ek' in reps and 'Et' not in reps:
            N = min(reps['Ek'].size, reps['Ep'].size)
            reps['Et'] = reps['Ek'][:N] + reps['Ep'][:N]
    
    results.reports = reps
    results.log = log


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
