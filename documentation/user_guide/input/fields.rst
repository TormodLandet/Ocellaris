.. _inp_fields:

Known fields
============

Known fields can be used to specify boundary and initial conditions. If they
are the same then one field can be used for both avoiding repeated definitions.
Known fields can also be used to specify custom fields that cannot easily be
described in C++, like the regular wave fields that are implemented. Known
fields can also be used to specify damping zones in a domain.

.. describe:: name

    Name of the field. Used to refer to the field from other places

.. describe:: type

    The field type, one of

    * AiryWaves
    * BlendedField
    * FreeSurfaceZone
    * MaxField
    * MinField
    * RaschiiWaves
    * ScalarField
    * VectorField
    * WaterBlock
    * WaveOutflow

.. describe:: stationary

    Set this to true/on (default is off) if you want to stop the field
    recomputing its value every time

.. describe:: plot

    Most fields support plotting their value to a file which can be enabled
    by this setting


AiryWaves
---------

Standard linear waves. This field exposes the following functions:

* field_name/``elevation``
* field_name/``c``
* field_name/``rho``
* field_name/``uhoriz``
* field_name/``uvert``
* field_name/``u`` - vector field, either [uhoriz, uvert] or [uhoriz, 0, uvert]
* field_name/``pdyn``
* field_name/``pstat``
* field_name/``ptot``

.. describe:: depth

    The water depth

.. describe:: depth_above

    Height of the air phase above the still water

.. describe:: still_water_position

    The location of the the still water plane. It will be equal to the depth if
    the origin is at the bottom.

.. describe:: omegas, periods, wave_lengths, wave_numbers

    List of numbers. You must specify one and only one of these

.. describe:: wave_phases

    The phase of the waves, list of numbers, default is 0 for all waves

.. describe:: amplitudes

    The wave amplitudes. List of numbers, no default values

.. describe:: current_speed, wind_speed

    A superinposed current in the water phase / air phase respectively

.. describe:: ramp_time

    Ramp up the amplitudes over a given time interval, default 0.

.. describe:: ramp_time

    Ramp up the amplitudes over a given time interval, default 0.

.. describe:: colour_projection_degree

    Project the colour function to DG0 using quadrature of this degree,
    default 6 (set degree to -1 to prevent this projection and just use
    interpolation instead)


RaschiiWaves
------------

Construct regular waves by use of the Rascii_ Python library to construct C++
code that describes the wave field

.. _Rascii: https://bitbucket.org/trlandet/raschii

.. describe:: wave_model

    The Raschii wave model to use, default is ``Fenton``. Available models as
    of January 2019 are:

    * Airy
    * Fenton
    * Stokes

.. describe:: air_model

    The Raschii air model to use, default is ``FentonAir``. Available models as
    of January 2019 are:

    * ConstantAir
    * FentonAir

.. describe:: wave_length, wave_height

    Wave parameters

.. describe:: model_order

    The order of the wave model, default is 5. Using 1 will always give Airy
    waves. Stokes waves are implemented up to order 5, Fenton waves can be
    calculated for arbitrary order, but order 10 is normally a good compromise
    between accuracy and speed.

.. describe:: depth, depth_above, ramp_time, still_water_position, current_speed

    Same as for AiryWaves

.. describe:: blending_height

    Distance used for blending water and air stream functions above the free
    surface in Raschii, see the `Rascii documentation
    <https://raschii.readthedocs.io/en/latest/index.html#documentation>`_ for
    more details.


ScalarField
-----------

.. describe:: variable_name

    The name of the scalar function that will be exposed, default ``phi``.

.. describe:: polynomial_degree

    The polynomial degree of the Continuous Lagrange function used to interpolate
    the field. Default 2

.. describe:: cpp_code

    The field description. You can use the coordinate vector, ``x[i]``, the
    time, ``t``, and any constants given in ``user_constants``, see
    :ref:`inp_user_code`.


VectorField
-----------

Same as scalar field, but ``cpp_code`` must be a list of C++ expressions and
the length of the list must be 2 in 2D and 3 in 3D.


BlendedField
------------

All variables, such as phi, are expressed as

.. math::

    \phi = (1 - b) \phi_0 + b \phi_1

where b is the scalar blending function with values [0, 1]

.. describe:: field0, field1

    Names of the two known fields to be blended (not field functions, just
    the name of the field). All field functions of the respective fields are
    blended.

.. describe:: blending_function

    The name of the scalar function used for blending

.. code-block:: yaml

    fields:
    -   name: x squared
        type: ScalarField
        cpp_code: pow(x[0], 2)
        stationary: true
    -   name: t squared
        type: ScalarField
        cpp_code: pow(t, 2)
    -   name: xt-blend
        type: BlendedField
        field0: x squared
        field0: t squared
        blending_function: x squared/phi


MaxField and MinField
---------------------

Return the max or min of two fields

.. describe:: field0, field1

    Names of the two known field functions

.. _inp_fields_freesurfacezone:

FreeSurfaceZone
---------------

A field that is 1.0 near the interface and 0.0 away from the interface

.. describe:: variable_name

    The name of the scalar function that will be exposed, default ``phi``.

.. describe:: radius

    The resulting field will be 1.0 inside a distance ``radius`` from the free
    surface and 0.0 outside two times the radius. Between these points there
    will be a smooth transition

        transition = 2 * r ** 3 - 9 * r ** 2 + 12 * r - 4


WaterBlock
----------

A block of water, projected to obtain the best representation of the colour
field c given a mesh that does not conform to the block. Inside the block the
field is 1.0, outside the block it is 0.0. The free surface is projected to be
as sharp as possible

.. describe:: variable_name

    The name of the scalar function that will be exposed, default ``c``.

.. describe:: xmin, xmax, ymin, ymax, zmin, zmax

    The extents of the block, all default to 0

.. describe:: polynomial_degree

    The polynomial degree of the resulting field, default 0, which gives a
    piecewise constant field

.. describe:: colour_projection_degree

    Project the colour function to DG0 using quadrature of this degree,
    default 6 (set degree to -1 to prevent this projection and just use
    interpolation instead)
