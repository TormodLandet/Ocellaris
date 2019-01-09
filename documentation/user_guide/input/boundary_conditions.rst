.. _inp_boundary_conditions:

Boundary conditions
===================

You need a list of boundary conditions for your problem. For each region of the
boundary you first need to tell Ocellaris how to find this region and then the
boundary conditions to apply to each of the variables. You must specify BCs for
the velocity for a single phase simulation and also for the pressure if you are
using a non-algebraic pressure correction method such as IPCS. Coupled methods
or algebraic pressure correction methods such as IPCS-A does not use boundary
conditions for the pressure at the same location as you give boundary
conditions for the velocity. You can of course use a outlet type boundary where
the pressure is given but not the velocity. Refer to a textbook about solving
the Navier-Stokes equations for details about what types of boundary conditions
are reasonable to expect to work.

You can select constant Dirichlet boundary conditions (``ConstantValue``) or
constant Neumann conditions (``ConstantGradient``). You can also have coded
boundary conditions where you give a source code snippet that is executed to
calculate the boundary condition value, either in Python (type ``CodedValue``)
or in C++ (type ``CppCodedValue``).

.. contents:: Contents
    :local:

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
        p:         # <-- only used in non-algebraic pressure correction methods
            type: ConstantGradient
            value: 0
    -   name: lid
        selector: code
        inside_code: on_boundary and x[1] >= 1.0 - 1e-8
        u:
            type: ConstantValue
            value: [1, 0]
        p:         # <-- only used in non-algebraic pressure correction methods
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

.. describe:: name

    The name of the region. For more helpful error messages etc.

.. describe:: selector

    How the region is selected. Supported methods are ``code`` and
    ``mesh_facet_region``.

.. describe:: inside_code

    Required when the selector is ``code``: Python code to mark facets as
    inside the region or not

.. describe:: mesh_facet_regions

    Required when the selector is ``mesh_facet_region``: List of identificator
    numbers of the facet regions from the mesh. See below.

.. describe:: map_code

    Required when using periodic boundary conditions: Code for mapping
    coordinates when using periodic boundary conditions. See below.
    Since periodic DG function spaces are currently not supported by FEniCS,
    the support for periodic boundary conditions may have bitrotted. It used to
    work for CG function spaces, but may not any more.


Selecting regions by code
-------------------------

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


Selecting regions from mesh facet regions
-----------------------------------------

If you load the mesh along with a facet region file you can select boundary
regions by referencing their number given in the facet region file. You can
select one or more mesh facet region per Ocellaris boundary region. In the
demo calculating flow around the 2D outline of an Ocellaris clownfish the
selection of the top and bottom wall is done as follows. Here 2 and 4 are the
numbers given to the top and bottom wall respectively in the Gmsh preprocessor
using :code:`Physical Line(3) =  {...}; Physical Line(5) =  {...};`:

.. code-block:: yaml

    boundary_conditions:
    -   name: Top and bottom
        selector: mesh_facet_region
        mesh_facet_regions: [3, 5]
        u:
            type: FreeSlip


Common to all boundary conditions
----------------------------------

The boundary condition for each variable is given in a sub-dictionary that has
the following options:

.. describe:: type

    What type of BC to apply. Currently the following are available:

    * ``ConstantValue``
    * ``ConstantGradient``
    * ``CodedValue``
    * ``CppCodedValue``
    * ``FieldFunction``
    * ``FreeSlip``
    * ``OpenOutletBoundary``, An open outlet zero stress boundary condition
    * ``ConstantRobin``
    * ``SlipLength``
    * ``InterfaceSlipLength``

.. describe:: enforce_zero_flux

    For Neumann and Robin boundaries where the value is not prescribed, but you
    want to ensure that nothing of the variable you are describing flows
    through the wall. This can be useful when the mesh should be a plane
    normal to the direction you are describing, but the mesh is not a perfect
    plane, but has some innacurracies causing the normal to be slightly off


BC type Constant
----------------

ConstantValue or ConstantGradient boundary conditions

.. describe:: value

    The value to apply when using ConstantXxxx. Either a scalar or a list of
    scalars.


BC type Python coded
--------------------

CodedValue or CodedGradient boundary conditions

.. describe:: code

    Python code to calculate the value

An example of coded boundary conditions can be seen in the the following which
applies the Taylor-Green vortex solution as Dirichlet boundary conditions:

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

Notice that there is a list of two code blocks for the velocity. Both are
evaluated as scalar fields and must assign to the zeroth component of the
:code:`value[]` array that is provided by FEniCS in order to set the Dirichlet
value at the boundary.

Note: If you write the boundary conditions in C++ instead of Python it will
normally be *significantly faster*.


BC type C++ coded
-----------------

CppCodedValue or CppCodedGradient boundary conditions

.. describe:: cpp_code

    C++ expression to calculate the value

An example of C++ boundary conditions can be seen in the the following which
applies the Taylor-Green vortex solution as Dirichlet boundary conditions:

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

Note that there is no assignment to the :code:`value[]` array. All math
functions from ``<cmath>`` are available as well as scalars like the time "t",
the timestep "dt", time index "it" and number of geometrical dimensions "ndim".
For single phase simulations "nu" and "rho" are also available.

You can use C++ lambda expressions to write multi-line expressions:

.. code-block:: yaml

    # ...
    cpp_code: |
        // This is the same as 'x[2] <= H ? 1.0 : 0.0'
        [&]() {
            bool is_above = x[2] > H;
            if (is_above) {
                return 0.0;
            } else {
                return 1.0;
            }
        }()


BC type known field function
----------------------------

Dirichlet BCs given by a known function

.. describe:: function

    The name of a known field function

Example from a wave simulation


.. code-block:: yaml

    boundary_conditions:
    -   name: Inlet
        selector: code
        inside_code: 'on_boundary and x[0] < 1e-5'
        u0:
            type: FieldFunction
            function: waves/uhoriz
        u1:
            type: ConstantValue
            value: 0
        u2:
            type: FieldFunction
            function: waves/uvert
        c:
            type: FieldFunction
            function: waves/c


BC type Constant Robin
-----------------------

Robin condition with constant values

    d(VAR)/dn = 1/blend (dval - VAR) + nval

.. describe:: blend

    A constant blending factor

.. describe:: dval

    The Dirichlet value

.. describe:: nval

    The Neumann value


BC type Slip length
--------------------

Wall slip length (Navier) boundary condition with constant value

.. describe:: value

    The value to prescribe, default 0

.. describe:: slip_length

    The slip length (the distance into the wall where the value is prescribed).


BC type Interface slip length
-----------------------------

Wall slip length (Navier) boundary condition where the slip length is
multiplied by a slip factor ∈ [0, 1] that varies along the domain boundary
depending on the distance to an interface (typically a free surface between two
fluids).

.. describe:: value, slip_length

    Same as above

.. describe:: slip_factor_function

    The function the blends from 0 slip length to full slip length. Typically
    the name of a :ref:`inp_fields_freesurfacezone` is given.
