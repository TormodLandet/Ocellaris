import traceback
import dolfin
from ocellaris.utils import timeit, verify_key


class Hooks(object):
    def __init__(self, simulation):
        """
        This class allows registering functions to run at
        given times during the simulation, e.g. to update
        some values for the next time step, report something
        after each time step or clean up after the simulation
        """
        self.simulation = simulation
        
        # Main hooks
        self._pre_simulation_hooks = []
        self._pre_timestep_hooks = []
        self._post_timestep_hooks = []
        self._post_simulation_hooks = []
        
        # More specialized hooks
        self._matrix_ready_hooks = []
        self._custom_hooks = {}
        
        # TODO: consider making all other hooks wrappers over custom hooks
        
    def register_custom_hook_point(self, hook_point):
        """
        Add a new custom hook points to which hooks can be added
        """
        assert hook_point not in self._custom_hooks
        self._custom_hooks[hook_point] = []
    
    # ------------------------------------------
    # Hook adders:
    
    def add_pre_simulation_hook(self, hook, description):
        """
        Add a function that will run before the simulation starts
        """
        self._pre_simulation_hooks.append((hook, description))
    
    def add_pre_timestep_hook(self, hook, description, timer_name=None):
        """
        Add a function that will run before the solver in each time step
        """
        if timer_name is None:
            timer_name = description
        self._pre_timestep_hooks.append((hook, description, timer_name))
    
    def add_post_timestep_hook(self, hook, description, timer_name=None):
        """
        Add a function that will run after the solver in each time step
        """
        if timer_name is None:
            timer_name = description
        self._post_timestep_hooks.append((hook, description, timer_name))
    
    def add_post_simulation_hook(self, hook, description):
        """
        Add a function that will run after the simulation is done
        """
        self._post_simulation_hooks.append((hook, description))
    
    def add_matrix_ready_hook(self, hook, description):
        """
        Add a function that will run after matrix assembly
        """
        self._matrix_ready_hooks.append((hook, description))
    
    def add_custom_hook(self, hook_point, hook, description):
        """
        Add other type of hook, must give a string name for
        the hook name
        """
        verify_key('custom hook point', hook_point, self._custom_hooks)
        self._custom_hooks[hook_point].append((hook, description))
    
    # ------------------------------------------
    # Hook runners:
    
    def simulation_started(self):
        """
        Called by the solver when the simulation starts
        
        Will run all pre simulation hooks in the reverse
        order they have been added
        """
        for hook, description in self._pre_simulation_hooks[::-1]:
            try:
                hook()
            except:
                self.simulation.log.error('Got exception in hook: %s' % description)
                self.simulation.log.error(traceback.format_exc())
                raise
        # Flush open files
        self.simulation.hooks.run_custom_hook('flush')
    
    @timeit.named('all hooks: new_timestep')
    def new_timestep(self, timestep_number, t, dt):
        """
        Called by the solver at the beginning of a new time step
        
        Will run all pre timestep hooks in the reverse
        order they have been added 
        """
        self.simulation._at_start_of_timestep(timestep_number, t, dt)
        for hook, description, timer_name in self._pre_timestep_hooks[::-1]:
            with dolfin.Timer('Ocellaris hook %s' % timer_name):
                try:
                    hook(timestep_number=timestep_number, t=t, dt=dt)
                except:
                    self.simulation.log.error('Got exception in hook: %s' % description)
                    self.simulation.log.error(traceback.format_exc())
                    raise
    
    @timeit.named('all hooks: end_timestep')
    def end_timestep(self):
        """
        Called by the solver at the end of a time step
        
        Will run all post timestep hooks in the reverse
        order they have been added
        """
        for hook, description, timer_name in self._post_timestep_hooks[::-1]:
            with dolfin.Timer('Ocellaris hook %s' % timer_name):
                try:
                    hook()
                except:
                    self.simulation.log.error('Got exception in hook: %s' % description)
                    self.simulation.log.error(traceback.format_exc())
                    raise
        self.simulation._at_end_of_timestep()
    
    def simulation_ended(self, success):
        """
        Called by the solver when the simulation is done
        
        Will run all post simulation hooks in the reverse
        order they have been added
        
        Arguments:
            success: True if nothing went wrong, False for
            diverging solution and other problems
        """
        self.simulation.success = success
        for hook, description in self._post_simulation_hooks[::-1]:
            try:
                hook(success=success)
            except:
                self.simulation.log.error('Got exception in hook: %s' % description)
                self.simulation.log.error(traceback.format_exc())
                raise
        # Flush open files
        self.simulation.hooks.run_custom_hook('flush')
    
    @timeit.named('all hooks: matrix_ready')
    def matrix_ready(self, Aname, A, b=None):
        """
        Called by the solver after assembly and before a linear
        solve. Can be used i.e for studies of condition numbers
        and reporting matrix sizes etc
        """
        for hook, description in self._matrix_ready_hooks[::-1]:
            t = dolfin.Timer('Ocellaris hook %s' % description)
            try:
                hook(Aname=Aname, A=A, b=b)
            except:
                self.simulation.log.error('Got exception in hook: %s' % description)
                self.simulation.log.error(traceback.format_exc())
                raise
            finally:
                t.stop()
                
    @timeit.named('all hooks: custom_hook')
    def run_custom_hook(self, hook_point, *args, **kwargs):
        """
        Called by the solver at a custom point
        
        Will run all custom hooks in the reverse order they have
        been added
        """
        verify_key('custom hook point', hook_point, self._custom_hooks)
        for hook, description in self._custom_hooks[hook_point][::-1]:
            try:
                hook(*args, **kwargs)
            except:
                self.simulation.log.error('Got exception in hook: %s' % description)
                self.simulation.log.error(traceback.format_exc())
                raise
    
    def show_hook_info(self):
        """
        Show all registered hooks
        """
        show = self.simulation.log.info
        all_hooks = [('Pre-simulation', self._pre_simulation_hooks),
                     ('Pre-timestep', self._pre_timestep_hooks),
                     ('Post-timestep:', self._post_timestep_hooks),
                     ('Post-simulation', self._post_simulation_hooks),
                     ('Matrix ready', self._matrix_ready_hooks)]
        all_hooks.extend(('Custom hook "%s"' % name, hooks) for name, hooks in self._custom_hooks.items())
        
        show('\nRegistered hooks:')
        for hook_type, hooks in all_hooks:
            show('    %s:' % hook_type)
            for info in hooks[::-1]:
                show('        - %s' % info[1])
