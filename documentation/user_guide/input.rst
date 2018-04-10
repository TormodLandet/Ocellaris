Ocellaris input file description
--------------------------------

.. contents:: Contents
    :local:

File format
...........

Ocellaris uses the YAML format for input files. The input file is divided
into separate sections dealing with geometry, boundary conditions, solver
parameters etc. The different sections are described below. Multiple demos
are provided along with Ocellaris and it is recommended to start with one
of the demo input files and use the below documentation as an aid to change
the demo input file into one that is describing your specific simulation.

Note that since JSON is a valid subset of YAML, you can also write the input
file in JSON format. JSON has no simple support for multi-line strings and
comments, so YAML is the format used by the Ocellaris demos and also in the
descriptions below.

Common errors
~~~~~~~~~~~~~

Some errors that are easy to make when writing a YAML input file:

- Boleans in YAML are written all lower case  (:code:`true, false`) unlike
  in Python where the first letter is upper case (:code:`True, False`). It
  can be easier to use the alternatives :code:`on` or :code:`off` so this
  confusion is avoided.
- The value ``5e-3`` is a string in YAML while ``5.0e-3`` is a float.
- Indentation is significant, just like in Python
- Misspellings are not checked!


Header
......

The input file **must** start with the following header:

.. code-block:: yaml

    ocellaris:
        type: input
        version: 1.0

You can *optionally* specify some metadata if you feel like it. This is not
required, but can be useful for explainations and later references.

.. code-block:: yaml

    metadata:
        author: Tormod Landet
        date: 2015-02-13
        description: |
            Free form text description of the input
            It can be quite usefull to have some text to 
            describe the purpose of the simulation etc for
            future reference 

Here you also see the syntax for multi-line strings in YAML.


Templates
~~~~~~~~~

You can specify a list of base input files that will be read first and used
as a basis for the input. Any values given in an input file will then extend
the template basis. This support is limited to key-value mappings. It is not
possible to replace parts of a list. Changing a list must be done by changing
the whole list in the derived input file.

Example base input file, ``base.inp``:

.. code-block:: yaml

    ocellaris:
        type: input
        version: 1.0
    
    user_code:
        constants:
            A: 2
    
    some_section:
        D: py$ A/B

The derived input file can use values defined in the base and extend it with
further data—you may need to read the `User constants and code`_ section to
fully  understand this example:

.. code-block:: yaml

    ocellaris:
        type: input
        version: 1.0
        bases:
        -   base.inp
    
    user_code:
        constants:
            B: 4
    
    some_section:
        C: py$ A*B

Ocellaris will interpret the input as:

.. code-block:: yaml

    ocellaris:
        type: input
        version: 1.0
    
    some_section:
        D: 0.5
        C: 8  


User constants and code
.......................

You can specify constants that can be used in subsequent sections to make
the input file easily configurable. You can also specify some code that
will run right after the input file has been read, before any of the 
simulation setup such as loading the mesh has been done. You can even
change the input by accessing the ``simulation.input`` object since no
parts of Ocellaris has accessed the input yet.

.. code-block:: yaml

    user_code:
        constants:
            L: 200       # channel length
            theta: 30    # angle
        code: |
            import subprocess
            subprocess.call(['command', 'to', 'generate', 'mesh'])
    
Example of using the constants in later sections of the input file:

.. code-block:: yaml

    some:
        section:
            param1: 4.3
            param2: py$ 2.3 * L * sin(theta)
            cpp_code: 'x[0] + L * sin(theta)' 

Any value (except inside the ``user_code/constants`` block) can be given as
a string starting with ``py$``. Ocellaris will then execute the given Python
code to produce the value to be used in Ocellaris just as if you had written
the value directly into the input file. The Python code you give can evaluate
to a list, string, number...

Code given as strings in the input file, either Python or C++ can also use
the constants as can be seen in the example. These are typically expressions
defining initial or boundary values. You can even combine these functions:

.. code-block:: yaml

    some-section:
        cpp_code: py$ 'x[0] + L * sin(theta)'.replace('theta', 'theta + L') 

This can be handy if you give the C++ code to compute the value of a field
as a user constant string, and then you can use python code to replace the
variable  ``t`` in the string with ``(t - dt)`` in order to specify the two
initial conditions, both at ``t=0`` and ``t=0-dt`` without having to repeat
the C++ code. This can, e.g., be used to describe a Taylor-Green vortex in
such a way that the time stepping can be second order from the first time
step (normally the first time setp is first order accurate since only one
initial condition is specified:


.. code-block:: yaml

    user_code:
        constants:
            u0a: '-sin(pi*x[1])*cos(pi*x[0])*exp(-2*pi*pi*nu*t)'
            u1a: ' sin(pi*x[0])*cos(pi*x[1])*exp(-2*pi*pi*nu*t)'
    
    initial_conditions:
        up0:
            cpp_code: py$ u0a
        up1:
            cpp_code: py$ u1a
        upp0:
            cpp_code: py$ u0a.replace('*t)', '*(t - dt))')
        upp1:
            cpp_code: py$ u1a.replace('*t)', '*(t - dt))')


Physical properties
...................

You will need to specify some physical constants. A simple example: 

.. code-block:: yaml

    physical_properties:
        g: [0, 0, 0]
        nu: 0.001
        rho: 1.0

.. describe:: g

    The acceleration of gravity given as a list of numbers. The length of the
    list must match the number of spatial directions, e.g. 2 or 3.
    Use ``[0, -9.81]`` in 2D and ``[0, 0, -9.81]`` in 3D for "standard" gravity.


Single phase properties
~~~~~~~~~~~~~~~~~~~~~~~

.. describe:: nu

    The kinematic viscosity

.. describe:: rho

    The density of the fluid, defaults to ``1.0``.


VOF two phase properties
~~~~~~~~~~~~~~~~~~~~~~~~

.. describe:: nu0, rho0

    The kinematic viscosity and density of fluid 0 

.. describe:: nu1, rho1

    The kinematic viscosity and density of fluid 1

For a water/air simulation fluid 0 is typically water and corresponds to
VOF colour function value ``1.0`` while fluid 1 is typically air and
corresponds to VOF colour function value ``0.0``. 


Variable density properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. describe:: nu

    The kinematic viscosity of both fluids (single value) 

.. describe::  rho_min, rho_max

    The range of allowable densities. Give one number for each of these settings.


Mesh
....

You can specify simple geometries using FEniCS DOLFIN built in mesh generators,
and also load a mesh from file. For realistic cases using something like gmsh
to generate meshes is recommended. The meshio_ program can be used to convert 
between different mesh file formats and also loading these formats directly,
see below.

.. _meshio: https://github.com/nschloe/meshio


Simple geometries
~~~~~~~~~~~~~~~~~

Example: 2D rectangle

.. code-block:: yaml
        
    mesh:
        type: Rectangle
        Nx: 64
        Ny: 64
        diagonal: left/right  # defaults to 'right'
        startx: 0             # defaults to 0
        endx:   2             # defaults to 1
        # you can also give starty and endy

Example: 3D box

.. code-block:: yaml
        
    mesh:
        type: Box
        Nx: 64
        Ny: 64
        Nz: 15
        startx: 0  # defaults to 0
        endx:   2  # defaults to 1
        # you can also give starty and endy, startz and endz

Example: 2D disc

.. code-block:: yaml
        
    mesh:
        type: UnitDisc
        N: 20
        degree: 1  # defaults to 1 (degree of mesh elements)


Mesh file formats
~~~~~~~~~~~~~~~~~

Example: legacy DOLFIN XML format

.. code-block:: yaml
        
    mesh:
        type: XML
        mesh_file: mesh.xml
        facet_region_file: regions.xml  # not required

Ocellaris will look for the xml files first as absolute paths, then as paths
relative to the current working directory and last as paths relative to the
directory of the input file. If it cannot find the file in any of these
places you will get an error message and Ocellaris will quit.

A sample mesh xml file and facet marker file is included in the ``demo/files``
directory. The mesh ``ocellaris_mesh.xml.gz`` and the facet regions
``ocellaris_facet_regions.xml.gz``. You can load these files without unzipping
them. The *flow around Ocellaris* demo shows how it is done.

Example: XDMF format

.. code-block:: yaml
        
    mesh:
        type: XDMF
        mesh_file: mesh.xdmf

Example: Ocellaris HDF5 restart file format

.. code-block:: yaml
        
    mesh:
        type: HDF5
        mesh_file: ocellaris_savepoint000010.h5

This will only load the mesh and (possibly) facet regions. You can also start
the simulation from a restart file instead of an input file. Then the mesh *and*
the function values from that save point are used, allowing you to restart the
simulation more or less like it was never stopped.

Example: using meshio_ to load all its supported formats

.. code-block:: yaml
        
    mesh:
        type: meshio
        mesh_file: mesh.msh
        meshio_type: gmsh

The supported formats (as of April 2018) can be found `in this list 
<https://github.com/nschloe/meshio/blob/f37dd29d67a35ec49a02bd131d1bf11ffde9ff85/meshio/helpers.py#L104>`_
in the meshio source on github.


Moving the mesh
~~~~~~~~~~~~~~~

Ocellaris can move the mesh right after it has been created or read from file.
To move the mesh in order to refine, skew, scale, rotate or translate it you
must specify a C++ description of the mesh *displacement* from the initial
position (which was specified in the input file or in the loaded mesh file).

An example is the following 140 meter long 2D wave tank which is 10 m high. To
refine the mesh in the y-direction such that it is finest around ``x[1] = 7``
meters—where the free surface is to be located—a function is specified which
is zero on the boundaries (to avoid changing the domain size) and non-zero in
the interior in order to move the nodes closer to the free surface. No refinement
is performed in the x-direction (``x[0]``).

.. code-block:: yaml
        
    mesh:
        type: Rectangle
        Nx: 140
        Ny: 20
        endx: 140
        endy: 20
        move: ['0', '0.0297619048*pow(x[1], 3) - 0.520833333*pow(x[1], 2) + 2.23214286*x[1] + 3.55271368e-15']

In order to develop and check the mesh refinement function it can be beneficial
to generate and plot it, e.g., using matplotlib in jupyter or using similar
interactive tools. The above refinement was developed using polynomial fitting
in numpy::

    from matplotlib import pyplot
    import numpy
    
    # Find a polynomial that refines the mesh
    y_target = [0, 4, 7.5, 10]
    dy_target = [0, 2.5, 0, 0]  # zero at the boundary
    P = numpy.polyfit(y_target, dy_target, 3)
    
    # Realise the polynomial
    y = numpy.linspace(0, 10, 20)
    dy = numpy.polyval(P, y)
    
    # Plot the results
    for ypos in (y + dy):
        pyplot.plot([0, 1], [ypos, ypos], '-k', lw=1)'
    pyplot.axhline(7, c='b', ls=':')
    pyplot.axhline(6, c='b', ls=':', lw=1)
    pyplot.axhline(8, c='b', ls=':', lw=1)

For more complicated meshes it is recommended to perform mesh grading and other
mesh operation in an external mesh generator such as gmsh. 
There is also some (not much used, hence possibly buggy) support for ALE where
the mesh moves every timestep, but that is not covered by the ``mesh`` section
of the input file.


Boundary conditions
...................

You need a list of boundary conditions for your problem. For each region of the
boundary you first need to tell Ocellaris how to find this region and then the
boundary conditions to apply to each of the variables (velocity and pressure for
a single phase simulation).

You can select constant Dirichlet boundary conditions (``ConstantValue``) or
constant Neumann conditions (``ConstantGradient``). You can also have coded
boundary conditions where you give a source code snippet that is executed to
calculate the boundary condition value, either in Python (type ``CodedValue``)
or in C++ (type ``CppCodedValue``). 

How to mark different areas of the boundary is explained below. For the lid
driven cavity the boundary conditions are as follows:

.. code-block:: yaml
                
    boundary_conditions:
    -   name: walls    
        selector: code
        inside_code: on_boundary
        u:
            type: ConstantValue
            value: [0, 0]
        p:
            type: ConstantGradient
            value: 0
    -   name: lid
        selector: code
        inside_code: on_boundary and x[1] >= 1.0 - 1e-8
        u:
            type: ConstantValue
            value: [1, 0]
        p:
            type: ConstantGradient
            value: 0

Note that the ``-`` in front of the ``name: ...`` lines marks the start of a
list item. The boundary conditions should be given as a list of boundary
regions. Each region specifies boundary conditions for all variables on the
selected boundary. 

The boundary conditions for the velocity components can also be broken up and
written per component. This allows you to apply different boundary conditions
types for each component. In this case it can be written (for the lid):
 
.. code-block:: yaml
    
    u0:
        type: ConstantValue
        value: 1
    u1:
        type: ConstantValue
        value: 0

Available options 
~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "key", "Default value", "Description"

    "boundary_conditions/[i]/name", "**required input**", "The name of the region. For more helpful error messages etc."
    "boundary_conditions/[i]/selector", "**required input**", "How the region is selected. Supported methods are ``code`` and ``mesh_facet_region``."
    "boundary_conditions/[i]/inside_code", "**required** when the selector is ``code``", "Python code to mark facets as inside the region or not"
    "boundary_conditions/[i]/mesh_facet_regions", "**required** when the selector is ``mesh_facet_region``", "List of identificator numbers of the facet regions from the mesh. See below."
    "boundary_conditions/[i]/map_code", "**required** when using periodic boundary conditions", "Code for mappinc coordinates when using periodic boundary conditions. See below."
    "boundary_conditions/[i]/var_name", "", "Boundary conditions for var_name. See below."

The boundary condition for each variable is given in a sub-dictionary that has
the following options:

.. csv-table::
   :header: "key", "Default value", "Description"

    "../var_name/type", "**required input**", "What type of BC to apply. Currently the following are available: ``ConstantValue``, ``ConstantGradient``, ``CodedValue`` and ``CppCodedValue``"
    "../var_name/value", "**required** when using ConstantXxxxx", "The value to apply. Either a scalar or a list of scalars."
    "../var_name/code", "**required** when using CodedXxxx", "Python code to calculate the value. Must be a multiline string that assigns to the value[i] variable (see below)"
    "../var_name/cpp_code", "**required** when using CppCodedXxxx", "C++ expression to calculate the value. Must evaluate to the requested value."

Selecting regions by code
~~~~~~~~~~~~~~~~~~~~~~~~~

You can select regions of the boundary by code in the same format as in FEniCS.
Ocellaris will run the Python code provided in the ``inside_code`` input key in
a statement equivalent to:

.. code-block:: python

    def boundary(x, on_boundary):
        return YOUR_REGION_CODE
        
if you give a single line expression, or

.. code-block:: python

    def boundary(x, on_boundary):
        YOUR_REGION_CODE
        return inside

if you give a multi line expression. In this case you need to assign a boolean
value to the name :code:`inside`.

How the inside_code works is that any facet where your code evaluates to
``True`` will be marked. As you can se above it is possible to mark everything
as is done for the walls and then overwrite this mark for parts of the boundary
as is done for the lid. The above will have walls everywhere below y=1 and lid
on y≥1. The FEniCS / dolfin syntax is used so ``x[0]`` is the x-component and 
``x[1]`` is the y-component.

Selecting regions from XML input 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you load the mesh along with a facet region file you can select boundary
regions by referencing their number given in the facet region file. You can
select one or more mesh facet region per Ocellaris boundary region. In the
demo calculating flow around the 2D outline of an Ocellaris clownfish the
selection of the top and bottom wall is done as follows. Here 2 and 4 are the
numbers given to the top and bottom wall respectively in the Gmsh preprocessor
using :code:`Physical Line(2) =  {...}; Physical Line(4) =  {...};`:

.. code-block:: yaml

    boundary_conditions:
    -   name: Top and bottom
        selector: mesh_facet_region
        mesh_facet_regions: [2, 4]
        u1:
            type: ConstantValue
            value: 0
        p:
            type: ConstantGradient
            value: 0

The above code applies a free-slip boundary condition on these two horisontal
walls. No boundary condition is applied in the tangential, ``u0``, direction.
Here it was necessary to split the velocity boundary condition into per
component boundary conditions.

Coded boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~

An example of coded boundary conditions can be seen in the the following which
applies the analytical Taylor-Green vortex solution as Dirichlet conditions:

.. code-block:: yaml

    boundary_conditions:
    -   name: walls
        selector: code
        inside_code: on_boundary
        u:
            type: CodedValue
            code:
            -   value[0] = -sin(pi*x[1]) * cos(pi*x[0]) * exp(-2*pi*pi*nu*t)
            -   value[0] =  sin(pi*x[0]) * cos(pi*x[1]) * exp(-2*pi*pi*nu*t)
        p:
            type: CodedValue
            code: value[0] = -(cos(2*pi*x[0]) + cos(2*pi*x[1])) * exp(-4.*pi*pi*nu*t)/4

Notice that there is a list of two code blocks for the velocity. Both are
evaluated as scalar fields and must assign to the zeroth component of the
:code:`value[]` array that is provided by FEniCS in order to set the Dirichlet
value at the boundary.

Boundary conditions can also be written in C++. If you write the boundary
conditions in C++ instead of Python it will normally be *significantly faster*.

The same example as above would be:

.. code-block:: yaml

    boundary_conditions:
    -   name: walls
        selector: code
        inside_code: on_boundary
        u:
            type: CppCodedValue
            cpp_code:
            -   -sin(pi*x[1]) * cos(pi*x[0]) * exp(-2*pi*pi*nu*t)
            -    sin(pi*x[0]) * cos(pi*x[1]) * exp(-2*pi*pi*nu*t)
        p:
            type: CppCodedValue
            cpp_code: -(cos(2*pi*x[0]) + cos(2*pi*x[1])) * exp(-4.*pi*pi*nu*t)/4

Note that there is no assignment to the :code:`value[]` array. All math
functions from ``<cmath>`` are available as well as scalars like the time "t",
the timestep "dt", time index "it" and number of geometrical dimensions "ndim".
For single phase simulations "nu" and "rho" are also available.


Initial conditions
..................

In the lid driven cavity test case both the velocity and the pressure fields
start from zero, so no initial values need to be given. The following is an
example of how to specify initial values for the Taylor-Green vortex on a 2D
square with side lengths equal to 2.0:

.. code-block:: yaml

    initial_conditions:
        up0:
            cpp_code: -sin(pi*x[1])*cos(pi*x[0])*exp(-2*pi*pi*nu*t)
        up1:
            cpp_code:  sin(pi*x[0])*cos(pi*x[1])*exp(-2*pi*pi*nu*t)
        p:
            cpp_code: -(cos(2*pi*x[0]) + cos(2*pi*x[1])) * exp(-4*pi*pi*nu*t)/4

.. csv-table::
   :header: "key", "Default value", "Description"

    "initial_conditions/var_name/cpp_code", "**required input**", "C++ code that gives the value of the field at each point. Variables ``rho``, ``nu`` and ``t`` are available"


Timestepping
............

This section sets the end time and time step. Currently only fixed time step is
available, though the time step can be altered in user coding at the expense of
slight errors in the treatment of the convecting velocity at the two time steps
following the change in time step:

.. code-block:: yaml
                     
    time:
        dt: 0.01
        tmax: 60.0

Example user code that changes the time step. See details under hooks below:

.. code-block:: yaml

    hooks:
        pre_timestep:
        -   name: decrease time step
            code: |
                if t > 10:
                    simulation.input['time']['dt'] = 0.005

Output control
..............

All the following parameters have sensible defaults and can be left out. The
output prefix can be useful to control in which directory the output files end
up. The final file name of all output files will be 
``output_prefix + file name``.

.. code-block:: yaml
        
    output:
        prefix: lid_driven_cavity_flow
        log_name: .log
        dolfin_log_level: warning
        ocellaris_log_level: info


.. csv-table::
   :header: "key", "Default value", "Description"

    "...", "**required input**", "FIXME: finish this table"


Linear solvers
..............

All equation systems that require global solves, like the velocity, pressure
and potentially multi phase models, will have their own optional definition of
the linear solver. These can be described in two ways, the simple FEniCS DOLFIN
based setup where some limited configuration is possible, or the full PETSc KSP
setup where all of the PETSc options are configurable plus a few options added
by Ocellaris.

It is recommended to use the KSP setup. It is the default, it is more powerfull
and it can do everything supported by the FEniCS DOLFIN setup. The DOLFIN setup
is kept for comparison and to be able to test the exact same setup used by
"normal" FEniCS codes.


PETSC KSP solver setup (use_ksp = yes)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This linear solver setup is used by most linear solvers inside Ocellaris. Most
solvers set reasonable defaults. Use these as starting points for your own
experimentations. The Ocellaris log file shows the setup which is used for the
different linear solvers in your simulation.

.. code-block:: yaml
        
    solver:
        u:
            use_ksp: yes
            petsc_ksp_type: gmres
            petsc_pc_type: asm
            petsc_ksp_initial_guess_nonzero: yes
            inner_iter_rtol: [1.0e-15, 1.0e-15, 1.0e-15]
            inner_iter_atol: [1.0e-15, 1.0e-15, 1.0e-15]
            inner_iter_max_it: [100, 100, 100]

.. describe:: use_ksp: yes

    Signal that we want to use the KSP solver setup (this is default in most
    situations).

.. describe:: petsc_XXXX

    Any PETSc parameter. Examples: ``ksp_type`` sets the solver name and
    ``pc_type`` sets the preconditioner name. Look at the PETSc documentation
    for the full list of tunable parameters, or give ``petsc_help: 'ENABLED'``
    to get a dump of possible parameters (the program will exit after giving
    the parameter listing). 

.. describe:: inner_iter_control

    The number of iterations and tolerances in the Krylov solver can be set for
    three categories of solves. The first X inner iterations (pressure
    correction iterations in the Navier-Stokes solver), the last Y inner
    iterations and the rest of the iterations (the middle number). The numbers
    X and Y are set by ``inner_iter_control: [X, Y]``. The default values are
    ``X=Y=3``.

.. describe:: inner_iter_rtol, inner_iter_atol, inner_iter_max_it

    The relative and absolute tolerances in the Krylov solver (default values
    are typically ``rtol = 1.0e-10`` and ``atol = 1.0e-15``). The maximum
    number of Krylov iterations is by default ``100`` for most solvers. If the
    solution is not converged the procedure will just continue, it is not
    always necessary to fully converge when applying an iterative solver, at
    least not in the inner first iterations (see below note on iterations).

.. note::

    Inner iterations refer to the main iterations inside each time step,
    typically pressure correction iterations (implemented in code inside
    Ocellaris). Krylov iterations refer to iterations inside the linear
    equation solver (provided by PETSc). The Krylov iterations are nested
    inside the inner iterations which are nested inside the time loop.


FEniCS DOLFIN solver setup (use_ksp = no)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml
        
    solver:
        u:
            use_ksp: no
            solver: gmres
            preconditioner: additive_schwarz
            parameters:
                any_parameter_supported_by_dolfin: valid_value

.. describe:: use_ksp: no

    Signal that we want to use the simplified setup

.. describe::  solver, preconditioner

    The names of the preconditioner and linear solver. Any values (string)
    supported by FEniCS DOLFIN are supported. The default values in FEniCS
    are used if none are specified (bad idea for large systems)

.. describe::  parameters

    Any parameter keys and values supported by FEniCS DOLFIN. See the DOLFIN
    documentation for these.


Navier-Stokes solvers
.....................

Some parameters are shared between all the available velocity-pressure solvers.
All the parameters in the following example have sensible defaults except for
the solver ``type`` which you must set. For the other parameters in the example
the values shown are the default values. The possible values for solver type
(IPCS-A, SIMPLE, PISO etc) are listed in the sections below with a brief
description. 

.. code-block:: yaml
    
    solver:
        type: IPCS-A
        num_inner_iter: 10
        allowable_error_inner: 1.0e-10
        polynomial_degree_pressure: 1
        polynomial_degree_velocity: 2
        function_space_pressure: DG
        function_space_velocity: DG
        u:
            # see linear solver documentation above
        p:
            # see linear solver documentation above

The inner iterations (pressure correction iterations) will run a maximum of
``num_inner_iter`` times for each time step, but the iterations will exit early
if the :math:`l^2` error of the difference between the predicted and corrected
velocity field is less than the given value ``allowable_error_inner``.

Some control parameters exist outside the common ones shown above, but none of
these are of the type that a normal user would probably need to change, so they
are only documented in the source code of the individual solvers.

The following parameters are relevant for under-relaxed solver implementations
(SIMPLE, PISO, PIMPLE):

.. describe::  relaxation_u, relaxation_p

    Relaxation factors. A value of 1.0 means no relaxation, 0.0 means no update
    at all (pointless). A value of 0.5 means that the result is an even blend
    of the computed value and the previous iteration value
    
.. describe::  relaxation_u_last_iter, relaxation_p_last_iter

    Some solvers will differentiate the last inner iteration from all other
    iterations. These parameters default to 1.0 in order to perform a "propper"
    update at the end of a time step with no relaxation applied.

IPCS-A
~~~~~~

Incremental Pressure Correction Scheme on Algebraic form. This is an iterative
Chorin/Temam type pressure correction solver.

IPCS-D
~~~~~~

Incremental Pressure Correction Scheme on Differential form. This is an
iterative Chorin/Temam type pressure correction solver where the pressure
correction Poisson equation is assembled from an elliptic operator and not
algebraicly from matrices. The divergence of the velocity field is hence not
very low and the method is not so strongly recommended for DG FEM, but it is
one of the most common solvers for the Navier-Stokes equations outside of DG
FEM and it has a smaller numerical stencil and may be faster than the IPCS-A
method.

SIMPLE
~~~~~~

Semi-Implicit Method for Pressure-Linked Equations. The implementation of the
algorithm is based on Klein, Kummer, Keil & Oberlack (2015).

PISO
~~~~

The pressure correction method by Issa (1986), Pressure-Implicit with Splitting
of Operators. PISO adds an additional correction step to the SIMPLE algorithm.

PIMPLE
~~~~~~

A Navier-Stokes solver based on the PIMPLE algorithm as implemented in OpenFOAM
and partially described in the PhD thesis of Jasak (1996; the PISO loop only).

.. describe::  num_pressure_corr

    The number of PISO iterations for each PIMPLE loop (the number of PIMPLE
    loops is controlled by the standard ``num_inner_iter`` parameter).

Coupled
~~~~~~~

Solves the velocity-pressure saddle point block-matrix equation system coupled.
Do not use this solver for large meshes. Even when using the multi-cpu
distributed multi frontal MUMPS or SuperLU_dist direct solvers there is a quite
small (perhaps around 1 million on a recent workstation?) limit to how many
degrees of freedom can be computed. For very small examples it may be faster
than using pressure-correction iterations and there is no resulting splitting
error which makes it great for testing and benchmarking the split solvers.

No block-system preconditioners are available in Ocellaris for the coupled
Navier-Stokes solver, so iterative linear solvers will either not converge or
perhaps "converge" to nonsensical solutions. Only use with direct solvers!

Analytical
~~~~~~~~~~

Use the initial condition C++ code (possibly containing the time variable ``t``
which will be updated for each time step) to define the velocity and pressure
for all time steps. This can be usefull for testing other parts of the 
Ocellaris solution framework with a known Navier-Stokes solution.


Multi phase solver
..................

If you are creating a two fluid simulation you will have to specify some
parameters of the multi-phase solver. For the lid driven cavity we can leave
the multi phase solver specification out of the input file. The default value 
of this section is:

.. code-block:: yaml

    multiphase_solver:
        type: SinglePhase

When using the multi phase VOF solver by specifying :code:`type: BlendedAlgebraicVOF`
the following parameters can be specified:

.. csv-table::
   :header: "key", "Default value", "Description"

    "multiphase_solver/function_space_colour", "DG", "CG for continuous Galerkin, DG for discontinuous Galerkin"
    "multiphase_solver/polynomial_degree_colour", "0", "The degree of the approximating polynomials"

In addition you will have to specify a convection scheme for the VOF colour
function in order to keep the free surface sharp. For specifying the convection
scheme, see below.


Convection
..........

Convecting fluxes have to be specified for all DG fields that are operated on
by a convection operator.

.. code-block:: yaml
                
    convection:
        u:
            convection_scheme: Upwind

.. csv-table::
   :header: "key", "Default value", "Description"

    "...", "**required input**", "FIXME: finish this table"

FIXME: describe HRIC/ CICSAM etc

Probes
......

Line probes can be added to sample the solution at each time step or at regular
intervals. Ocellaris can also show a plot of the sampled probe values that it
will update while it is running so that you can visually inspect the solution.

.. code-block:: yaml
        
    probes:
    -   name: u-vel center
        type: LineProbe
        field: u0
        startpos: [0.5, 0]
        endpos: [0.5, 1]
        Npoints: 100
        file_name: _uprobe.txt
        show_interval: 1
        write_interval: 10
        target_name: Ghia et al
        target_abcissa: [1.0, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5,
                         0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0]
        target_ordinate: [1, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719, 0.05702, -0.0608,
                          -0.10648, -0.27805, -0.38289, -0.2973, -0.2222, -0.20196, -0.18109, 0]
        
      
    -   name: v-vel center
        type: LineProbe
        field: u1
        startpos: [0, 0.5]
        endpos: [1, 0.5]
        Npoints: 100
        
        file_name: _vprobe.txt
        write_interval: 10
        
        target_abcissa: [1.0, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5,
                         0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0]
        target_name: Ghia et al
        target_ordinate: [0, -0.21388, -0.27669, -0.33714, -0.39188, -0.5155, -0.42665, -0.31966,
                          0.02526, 0.32235, 0.33075, 0.37095, 0.32627, 0.30353, 0.29012, 0.27485, 0.0]


.. csv-table::
   :header: "key", "Default value", "Description"

    "...", "**required input**", "FIXME: finish this table"


User code / hooks
.................

TODO: describe this. See example under timestepping above for now.

.. csv-table::
   :header: "key", "Default value", "Description"

    "...", "**required input**", "FIXME: finish this table"

The example below shows that each hook gets it's own dictionary ``hook_data``
to store whatever it wants between calls. The example also shows how to read
the input file parameters in a hook that is defined in the same input file, and
how to perform output to file in a configurable manner:

.. code-block:: yaml

    -   name: save colour function field
        enabled: yes
        code: |
            if not 'cf' in hook_data:
                prefix = simulation.input.get_value('output/prefix')
                hook_data['cf'] = File(prefix + '_c.pvd')
            if t > 1:
                hook_data['cf'] << (c, t)
