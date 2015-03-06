User guide
==========

The following is a quick introduction to setting up an Ocellaris simulation.
Note that many (MANY) possible simulations can be set up that Ocellaris will
fail to solve properly. You must validate the code for your own purposes before
trusting the results.

.. contents:: :local:


The input file
---------------

To run Ocellaris you must create an input file. The input file is on JSON
format and allows you to control most of the solution process.

All strings must use double quotes ``"`` not single quotes ``'`` in JSON.
The JSON parser in Python is quite strict, so it is a high likelihood that it 
will complain about missing or extraneous commas or parenthesis the first time 
you try. Another gotcha is that booleans in JSON are written all lower case 
(:code:`true, false`) unlike in Python where the first letter is upper case 
(:code:`True, False`).

Don't forget to close the any brackets and remember that there must be no 
comma after the last item in a lists or dictionary.


Header
......

The input file must start with the following header:

.. code-block:: javascript

    {
        "program": "ocellaris",
        "type": "input",
        "version": 1.0,
        
Physical properties
...................


You may then specify some physical constants. The postfix ``0`` is there to 
allow for more than one fluid in one simulation. The next fluid would use
postfix ``1`` and so on.

.. code-block:: javascript

        "physical_properties": {
            "rho0": 1.0,
            "nu0": 0.001,
            "g": [0, 0]
        },

You will want to load a mesh. Currently only 2D rectangle meshes are supported,
but it is planned that ``type`` can be set to some other value to for instance
load a mesh from a file.

.. code-block:: javascript
        
        "mesh": {
            "type": "Rectangle",
            "Nx": 64,
            "Ny": 64
        },

You can also specify ``startx``, ``endx`` and ``starty``, ``endy``. 

Boundary conditions
...................

You need boundary conditions on your mesh. Only constant value boundary
conditions are implemented currently, but this is very easy to extend. You can
select constant Dirichlet boundary conditions (``ConstantValue``) or constant
Neumann conditions (``ConstantGradient``). You must provide a vector for the
velocity and a scalar for the pressure.

How to mark different areas of the boundary is explained below.
 
.. code-block:: javascript
                
        "boundary_conditions": [
            {
                "name": "walls",
                "selector": "region",
                "region_code": "True",
                
                "p": {
                    "type": "ConstantGradient",
                    "value": 0
                },
                "u": {
                    "type": "ConstantValue",
                    "value": [0, 0]
                }
            },
            {
                "name": "lid",
                "selector": "region",
                "region_code": "x[1] >= 1.0 - 1e-8",
                
                "p": {
                    "type": "ConstantGradient",
                    "value": 0
                },
                "u": {
                    "type": "ConstantValue",
                    "value": [1, 0]
                }
            }
        ],
 
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
                     
        "time": {
            "dt": 0.01,
            "tmax": 60.0
        },


Output control
..............

All the following parameters have sensible defaults and can be left out. The
output prefix can be useful to control in which directory the output files end
up. The final file name of all output files will be 
``output_prefix + file name``.

.. code-block:: javascript

        "output": {
            "prefix": "precond_DG_",
            "log_name": "log.txt",
            "plot_at_end": false,
            "dolfin_log_level": "warning",
            "ocellaris_log_level": "info"
        },


The solver
..........

All the following parameters have sensible defaults. They all control the 
solution process in one way or the other. See the FEniCS documentation for the
available selection of solvers and preconditioners.

The inner iterations will run maximum ``num_inner_iter`` times, but will exit
eraly if the :math:`L^\inf` error of the difference between the predicted and
correcte dvelocity field is less than a given value ``allowable_error_inner``.

.. code-block:: javascript
        
        "solver": {
            "type": "IPCS",
            "num_inner_iter": 20,
            "allowable_error_inner": 5e-3,
            "function_space_velocity": "DG",
            "function_space_pressure": "DG",
            "polynomial_degree_velocity": 2,
            "polynomial_degree_pressure": 1,
            "calculate_L2_norms": false,
            "p": {
                "solver": "gmres",
                "preconditioner": "hypre_amg"
            },
            "u": {
                "solver": "bicgstab",
                "preconditioner": "bjacobi"
            }
        },


Convection
..........

Convecting fluxes have to be specified for all DG fields that are operated on
by a convection operator.

.. code-block:: javascript
                
        "convection": {
            "u": {
                "convection_scheme": "Upwind"
            }
        },


Probes
......

Line probes can be added to sample the solution at each time step or at regular
intervals. Ocellaris can also show a plot of the sampled probe values that it
will update while it is running so that you can visually inspect the solution.

.. code-block:: javascript
    
        "probes": [
            {
                "name": "u-vel center",
                "type": "LineProbe",
                "startpos": [0.5, 0],
                "endpos": [0.5, 1],
                "Npoints": 100,
                "field": "u0",
                "show_interval": 10,
                "write_interval": 10,
                "file_name": "uprobe.txt",
                "target_name": "Ghia et al",
                "target_abcissa": [1.0, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5,
                                   0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0],
                "target_ordinate": [1, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719, 0.05702, -0.06080, -0.10648,
                                    -0.27805, -0.38289, -0.29730, -0.22220, -0.20196, -0.18109, 0]
            },
            {
                "name": "v-vel center",
                "type": "LineProbe",
                "startpos": [0, 0.5],
                "endpos": [1, 0.5],
                "Npoints": 100,
                "field": "u1",
                "file_name": "vprobe.txt",
                "write_interval": 10,
                "target_name": "Ghia et al",
                "target_abcissa": [1.0, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5,
                                   0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0],
                "target_ordinate": [0, -0.21388, -0.27669, -0.33714, -0.39188, -0.51550, -0.42665, -0.31966, 
                                    0.02526, 0.32235, 0.33075, 0.37095, 0.32627, 0.30353, 0.29012, 0.27485, 0.0]
            }
        ],


Complete input file example
---------------------------

The following is the input file of the lid driven cavity demo at the time this
documentation was generated:

.. literalinclude:: ../../demos/lid_driven_cavity_flow.inp
   :language: python
