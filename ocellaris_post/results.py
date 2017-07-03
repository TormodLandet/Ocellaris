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
        self.file_name = os.path.abspath(file_name)
        self.input = None
        self.reports = None
        
        if file_name.endswith('.h5'):
            read_h5_data(self, derived)
        elif file_name.endswith('.log'):
            read_log_data(self, derived)


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
    
    results.reports = reps


def read_log_data(results, derived):
    INP_START = '----------------------------- configuration begin -'
    INP_END = '------------------------------ configuration end -'
    in_input_section = False
    input_strs = []
    data = {}
    
    # Read input and timestep reports from log file
    for line in open(results.file_name, 'rt'):
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
