import dolfin


class Log(object):
    # Names for the available log levels in Dolfin and Ocellaris
    AVAILABLE_LOG_LEVELS = {'critical': dolfin.CRITICAL,
                            'error': dolfin.ERROR,
                            'warning': dolfin.WARNING,
                            'info': dolfin.INFO,
                            'progress': dolfin.PROGRESS,
                            'debug': dolfin.DEBUG}
    
    def __init__(self, simulation):
        self.simulation = simulation
        self.log_level = dolfin.INFO
        self.simulation.hooks.add_post_simulation_hook(lambda success: self.end_of_simulation(), 'Flush log file')
        self.write_log = False
        self.write_stdout = False
    
    def write(self, message):
        """
        Write a message to the log without checking the log level
        """
        message = str(message) # hypotetically a collective operation ...
        if self.write_log:
            self.log_file.write(message + '\n')
        if self.write_stdout:
            print message
    
    def set_log_level(self, log_level):
        """
        Set the Ocellaris log level
        (not the dolfin log level!)
        """
        self.log_level = log_level
    
    def error(self, message):
        "Log an error message"
        if self.log_level <= dolfin.ERROR:
            self.write(message)
    
    def warning(self, message=''):
        "Log a warning message"
        if self.log_level <= dolfin.WARNING:
            self.write(message)
    
    def info(self, message=''):
        "Log an info message"
        if self.log_level <= dolfin.INFO:
            self.write(message)
    
    def debug(self, message=''):
        "Log a debug message"
        if self.log_level <= dolfin.DEBUG:
            self.write(message)
    
    def setup(self):
        """
        Setup logging to file if requested in the simulation input
        """
        log_name = self.simulation.input.get_output_file_path('output/log_name', None)
        log_on_all_ranks = self.simulation.input.get_value('output/log_on_all_ranks', False, 'bool')
        stdout_on_all_ranks = self.simulation.input.get_value('output/stdout_on_all_ranks', False, 'bool')
        rank = self.simulation.rank
        
        self.write_stdout = (rank == 0 or stdout_on_all_ranks)
        self.write_log = False
        if log_name is not None:    
            if log_on_all_ranks and rank > 0:
                log_name = '%s.%d' % (log_name, self.simulation.rank)
            
            if rank == 0 or log_on_all_ranks:
                self.write_log = True
                self.log_file_name = log_name
                self.log_file = open(self.log_file_name, 'wt')
        
        # Set the Ocellaris log level
        log_level = self.simulation.input.get_value('output/ocellaris_log_level', 'info')
        self.simulation.log.set_log_level(self.AVAILABLE_LOG_LEVELS[log_level])
        
        # Set the Dolfin log level
        df_log_level = self.simulation.input.get_value('output/dolfin_log_level', 'warning')
        dolfin.set_log_level(self.AVAILABLE_LOG_LEVELS[df_log_level])
    
    def end_of_simulation(self):
        """
        The simulation is done. Make sure the output
        file is flushed, but keep it open in case 
        some more output is coming
        """
        if self.write_log:
            self.log_file.flush()
