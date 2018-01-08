import os, sys
import dolfin


ALWAYS_WRITE = 1e10
NO_COLOR = '%s'
RED = '\033[91m%s\033[0m'    # ANSI escape code Bright Red
YELLOW = '\033[93m%s\033[0m' # ANSI escape code Bright Yellow


class Log(object):
    # Names for the available log levels in Dolfin and Ocellaris
    AVAILABLE_LOG_LEVELS = {'all': ALWAYS_WRITE,
                            'critical': dolfin.LogLevel.CRITICAL,
                            'error': dolfin.LogLevel.ERROR,
                            'warning': dolfin.LogLevel.WARNING,
                            'info': dolfin.LogLevel.INFO,
                            'progress': dolfin.LogLevel.PROGRESS,
                            'debug': dolfin.LogLevel.DEBUG}
    
    def __init__(self, simulation):
        self.simulation = simulation
        self.log_level = dolfin.LogLevel.INFO
        self.simulation.hooks.register_custom_hook_point('flush')
        self.simulation.hooks.add_custom_hook('flush', self.flush, 'Flush log file')
        self.write_log = False
        self.write_stdout = False
        self._the_log = []
    
    def write(self, message, msg_log_level=ALWAYS_WRITE, color=NO_COLOR):
        """
        Write a message to the log without checking the log level
        """
        message = str(message) # hypotetically a collective operation ...
        
        if self.log_level <= msg_log_level:
            if self.write_log:
                self.log_file.write(message + '\n')
            if self.write_stdout:
                print(color % message)
        
        # Store all messages irrespective of the log level
        self._the_log.append(message)
    
    def set_log_level(self, log_level):
        """
        Set the Ocellaris log level
        (not the dolfin log level!)
        """
        self.log_level = log_level
    
    def error(self, message):
        "Log an error message"
        self.write(message, dolfin.LogLevel.ERROR, RED)
    
    def warning(self, message=''):
        "Log a warning message"
        self.write(message, dolfin.LogLevel.WARNING, YELLOW)
    
    def info(self, message=''):
        "Log an info message"
        self.write(message, dolfin.LogLevel.INFO)
    
    def debug(self, message=''):
        "Log a debug message"
        self.write(message, dolfin.LogLevel.DEBUG)
    
    def setup(self):
        """
        Setup logging to file if requested in the simulation input
        """
        # Ensure that the output directory exist
        tmp = self.simulation.input.get_output_file_path('BOGUS', 'xxxxx')
        output_dir = os.path.split(tmp)[0]
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        log_name = self.simulation.input.get_output_file_path('output/log_name', '.log')
        log_on_all_ranks = self.simulation.input.get_value('output/log_on_all_ranks', False, 'bool')
        log_enabled = self.simulation.input.get_value('output/log_enabled', True, 'bool')
        log_append_existing = self.simulation.input.get_value('output/log_append_to_existing_file', True, 'bool')
        stdout_on_all_ranks = self.simulation.input.get_value('output/stdout_on_all_ranks', False, 'bool')
        stdout_enabled = self.simulation.input.get_value('output/stdout_enabled', True, 'bool')
        rank = self.simulation.rank
        
        self.write_stdout = (rank == 0 or stdout_on_all_ranks) and stdout_enabled
        self.write_log = False
        if log_enabled:
            if log_on_all_ranks and rank > 0:
                log_name = '%s.%d' % (log_name, self.simulation.rank)
            
            if rank == 0 or log_on_all_ranks:
                self.write_log = True
                self.log_file_name = log_name
                if log_append_existing:
                    self.log_file = open(self.log_file_name, 'at')
                    self.log_file.write('\n\n')
                else:
                    self.log_file = open(self.log_file_name, 'wt')
        
        # Set the Ocellaris log level
        log_level = self.simulation.input.get_value('output/ocellaris_log_level', 'info')
        self.simulation.log.set_log_level(self.AVAILABLE_LOG_LEVELS[log_level])
        
        # Set the Dolfin log level
        df_log_level = self.simulation.input.get_value('output/dolfin_log_level', 'warning')
        dolfin.set_log_level(self.AVAILABLE_LOG_LEVELS[df_log_level])
    
    def flush(self):
        """
        The simulation has started, flush to make sure
        input values are shown
        
        or
        
        The simulation is done. Make sure the output
        file is flushed, but keep it open in case 
        some more output is coming
        """
        if self.write_log:
            self.log_file.flush()
        if self.write_stdout:
            sys.stdout.flush()
    
    def get_full_log(self):
        """
        Get the contents of all logged messages as a string
        """
        return '\n'.join(self._the_log)
