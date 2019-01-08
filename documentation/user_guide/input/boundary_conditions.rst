.. _inp_boundary_conditions:

Boundary conditions
===================

You need a list of boundary conditions for your problem. For each region of the
boundary you first need to tell Ocellaris how to find this region and then the
boundary conditions to apply to each of the variables (velocity and pressure
for a single phase simulation).

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
-----------------

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

Selecting regions from XML input
--------------------------------

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
-------------------------

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
