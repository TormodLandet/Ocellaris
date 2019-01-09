Ocellaris input file description
================================


.. contents:: Contents
   :local:


Getting started
---------------

Starting from one of the demos is recommended. Most of the input options used
in the demo input files (and most other input options) are explained in
:ref:`inp_file_sections` below.

As a general rule, based on what parts of Ocellaris are most used and best
tested, the following options are good starting points for creating your own
simulation input files:

* Use a simple mesh (unit square or cube etc) or load your mesh with the meshio
  reader, see :ref:`inp_mesh`.
* Use the IPCS-A solver. If you are running on only one CPU then consider the
  Coupled solver (but it lacks some functionality comapred to IPCS-A), see
  :ref:`inp_solver` for the full list.
* Use single phase or the standard algebraic VOF multi phase solver, see
  ``BlendedAlgebraicVOF`` under :ref:`inp_multiphase_solver` for details.
* Write restart files once in a while to be able to restart your long running
  simulations should something hang, see :ref:`inp_output` for details on this.
  You should also consider using the ``orun`` script to "babysit" your
  simulations, see :ref:`script_orun`.


File format
-----------

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
.............

Some errors that are easy to make when writing a YAML input file:

- Boleans in YAML are written all lower case  (:code:`true, false`) unlike
  in Python where the first letter is upper case (:code:`True, False`). It
  can be easier to use the alternatives :code:`on` or :code:`off` so this
  confusion is avoided.
- The value ``5e-3`` is a string in YAML while ``5.0e-3`` is a float.
- Indentation is significant, just like in Python
- Misspellings are not errors!

The input you specify is validated against a schema, but you only get a warning
in the top of the Ocellaris simulation log file and no errors since the schema
is not guaranteed to be perfect. The input file YSchema_ file
(``ocellaris/input_file_schema.yml`` in the source code) is an approximate
representation of all valid and invalid Ocellaris input files and it is not
trusted to reject files that it believes are not valid. The warning can still
help catching some misspelling errors, but you must read the log file to see
them. When performing a new type of simulation it is always **strongly**
recommended to thoroughly read through the start of the log file for warnings
from the schema validator and from the many other modules of Ocellaris that
will warn you about fishy input.

.. _Yschema: https://bitbucket.org/trlandet/yschema


Header
------

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
.........

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
further data—you may need to read the :ref:`inp_user_code` section to fully
understand this example:

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


.. _inp_file_sections:

Input file sections
-------------------

.. toctree::

    user_code
    physical_properties
    mesh
    boundary_conditions
    initial_conditions
    time
    solver
    multiphase_solver
    convection
    slope_limiter
    fields
    forcing_zones
    output
    probes
    hooks
