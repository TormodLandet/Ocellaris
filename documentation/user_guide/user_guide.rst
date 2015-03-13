User guide
==========

The following is a quick introduction to setting up an Ocellaris simulation.
Note that many (MANY) possible simulations can be set up that Ocellaris will
fail to solve properly. You must validate the code for your own purposes before
trusting the results.

.. contents:: :local:


The input file
---------------

To run Ocellaris you must create an input file. The input file is on YAML
format and allows you to control most of the solution process. The different
sections of the input file are given below. 

Note that since JSON is a valid subset of YAML you can also specify the input
file in JSON format. JSON has no easy support for multi-line strings and
comments, so YAML is the format used by the Ocellaris demos and also in the
descriptions below.

Gotchas
.......

Some errors that are easy to make when writing a YAML input file:

- Boleans in YAML are written all lower case  (:code:`true, false`) unlike
  in Python where the first letter is upper case (:code:`True, False`).
- The value ``5e-3`` is a string in YAML while ``5.0e-3`` is a float.
- Indentation is significant, just like in Python

Header
......

The input file must start with the following header:

.. code-block:: yaml

    ocellaris:
        type: input
        version: 1.0

You can optionally specify some metadata if you feel like it. This is not
required, but can be useful for explainations and later references.

.. code-block:: yaml

    metadata:
        author: Tormod Landet
        date: 2015-02-13
        description: |
            Implements the lid driven cavity flow test case. This benchmark case
            has no analytical solution, but is a nice way to test the solver in
            presence of viscous shear and discontinuities.
            
            Comparison data is included. This is from the paper by V. Ghia, 
            K.N. Ghia and C.T. Shin: "High-Re solutions for incompressible flow
            using the Navier-Stokes equations and a multi-grid method" in 
            J.Comp.Phys., v.48, 1982, pp.387-411. 

Here you also see the syntax for multi-line strings in YAML.

Physical properties
...................

You will need to specify some physical constants:

.. code-block:: yaml

    physical_properties:
        g: [0, 0]
        nu0: 0.001
        rho0: 1.0

The postfix ``0`` is there to  allow for more than one fluid in one simulation.
The next fluid would use postfix ``1`` and so on.

The computational mesh
......................

You will want to load or create a mesh. Currently only 2D rectangle meshes are
supported, but it is planned that ``type`` can be set to some other value to for
instance load a mesh from a file.

.. code-block:: yaml
        
    mesh:
        type: Rectangle
        Nx: 64
        Ny: 64

You can also specify ``startx``, ``starty``, ``endx`` and ``endy``. 

Boundary conditions
...................

You need boundary conditions on your mesh. Only constant value boundary
conditions are implemented currently, but this is very easy to extend. You can
select constant Dirichlet boundary conditions (``ConstantValue``) or constant
Neumann conditions (``ConstantGradient``). You must provide a vector for the
velocity and a scalar for the pressure.

How to mark different areas of the boundary is explained below.
 
.. code-block:: yaml
                
    boundary_conditions:
    -   name: walls    
        selector: region
        inside_code: on_boundary and True
        u:
            type: ConstantValue
            value: [0, 0]
        p:
            type: ConstantGradient
            value: 0
    -   name: lid
        selector: region
        inside_code: on_boundary and x[1] >= 1.0 - 1e-8
        u:
            type: ConstantValue
            value: [1, 0]
        p:
            type: ConstantGradient
            value: 0

.. warning:: 

    The following description is outdated and need an update

Currently you can only select regions of the boundary by code on the same
format as in FEniCS. Ocellaris will run the Python code provided in the
``region_code`` input key in a statement equivalent to:

.. code-block:: python

    def boundary(x, on_boundary):
        return on_boundary and YOUR_REGION_CODE

This means that any part of the boundary where your code evaluates to ``True``
will be marked. As you can se above it is possible to mark everything as is
done for the walls and then overwrite this mark for parts of the boundary as is
done for the lid. The above will have walls everywhere below y=1 and lid on 
y≥1. The FEniCS / dolfin syntax is used so ``x[0]`` is the x-component and 
``x[1]`` is the y-component.


Timestepping
............

The following should be quite self-explanatory:

.. code-block:: javascript
                     
    time:
        dt: 0.01
        tmax: 60.0

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

The solver
..........

All the following parameters have sensible defaults. They all control the 
solution process in one way or the other. See the FEniCS documentation for the
available selection of solvers and preconditioners.

The inner iterations will run maximum ``num_inner_iter`` times, but will exit
early if the :math:`L^\infty` error of the difference between the predicted and
corrected velocity field is less than a given value ``allowable_error_inner``.

.. code-block:: yaml
    
    solver:
        type: IPCS
        num_inner_iter: 20
        allowable_error_inner: 5.0e-3
        polynomial_degree_pressure: 1
        polynomial_degree_velocity: 2
        function_space_pressure: DG
        function_space_velocity: DG
        timestepping_method: BDF
    
Convection
..........

Convecting fluxes have to be specified for all DG fields that are operated on
by a convection operator.

.. code-block:: yaml
                
    convection:
        u:
            convection_scheme: Upwind


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


Complete input file example
---------------------------

The following is the input file of the lid driven cavity demo at the time this
documentation was generated:

.. literalinclude:: ../../demos/lid_driven_cavity_flow.inp
   :language: yaml
