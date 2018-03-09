.. _sec-simulation:

Simulation classes
==================

The simulation classes holds the simulation data (velocity, pressure, function 
spaces, mesh etc) and is also responsible for most utility functionality such
as plugins (hooks), logging, reporting, plotting and input file handling. 

.. autoclass:: ocellaris.Simulation

    .. attribute:: input
        :annotation:
        
        An :class:`ocellaris.simulation.input.Input` object that holds the input from
        the input file provided by the user
    
    .. attribute:: io
        :annotation:
        
        An :class:`ocellaris.simulation.io.InputOutputHandling` object that holds 
        supports writing VTK, XDMF, restart HDF5 files etc
    
    .. attribute:: hooks
        :annotation:
        
        An :class:`ocellaris.simulation.hooks.Hooks` object that keeps track of
        functions that should run at certain times during the simulation

    .. attribute:: reporting
        :annotation:
        
        An :class:`ocellaris.simulation.reporting.Reporting` object that helps report
        summary values each time step

    .. attribute:: log
        :annotation:
        
        An :class:`ocellaris.simulation.log.Log` object that helps with logging
        messages to screen while the simulation is running
    
    .. attribute:: data
        :annotation: a dictionary
        
        In a typical simulation the :attr:`data` attribute will contain contents
        such as:
        
        .. list-table:: 
            :header-rows: 1
            :widths: 30, 70
            
            * - Name
              - Value
            * - boundary
              - 
            * - boundary_marker
              -
            * - cell_info
              -
            * - connectivity_CC
              - 
            * - connectivity_CF
              -
            * - connectivity_FC
              -
            * - constrained_domain
              - A SubDomain object that reporesents the periodic boundary conditions, or None if no periodic boundary conditions are applied
            * - dirichlet_bcs
              - A dictionary of lists of boundary conditions {'u0': [...], 'p': [...], ...}
            * - ds
              - A measure of the boundaries that has been marked by the boundary conditions 
            * - facet_info
              - 
            * - g
              - The acceleration of gravity (vector constant)
            * - mesh
              - The FEniCS mesh object
            * - neumann_bcs
              - A dictionary of lists of boundary conditions {'u0': [...], 'p': [...], ...}
            * - nu
              - The kinematic viscosity (constant or scalar valued function)
            * - p
              - The pressure at the current time step 
            * - p_hat
              - The pressure correction used to get to the current time step
            * - rho
              - The density (constant or scalar valued function)
            * - u, u0, u1, u2
              - The velocity at the current time stepas a vector and as components 
            * - u_conv, u_conv0, u_conv1, u_conv2
              - The extrapolated convecting velocity at the current time step as a vector and as components
            * - u_star, u_star0, u_star1, u_star2
              - The predicted velocity at the current time step as a vector and as components
            * - up, up0, up1, up2
              - The velocity at the previous time step as a vector and as components
            * - upp, upp0, upp1, upp2
              - The velocity at the previous-previous time step as a vector and as components
            * - Vp
              - The function space of velocities
            * - Vp
              - The function space of pressures
            

.. autoclass:: ocellaris.simulation.input.Input

.. autoclass:: ocellaris.simulation.io.InputOutputHandling

.. autoclass:: ocellaris.simulation.hooks.Hooks

.. autoclass:: ocellaris.simulation.reporting.Reporting

.. autoclass:: ocellaris.simulation.log.Log
