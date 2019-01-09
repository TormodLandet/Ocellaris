.. _inp_forcing_zones:

Forcing zones
=============

Forcing zones can be used to damp waves near the inlet and outlet to avoid
reflections. They work by specifying a zone (a known field), a target result
(a known field) and a penalty to be added to the equation system in order to
pull the unknown solutiuon towards the known field solution in the forcing
zone. See, e.g., :cite:`peric_amm_2016` and :cite:`peric_amm_2018` for more
details on forcing zones and the selection of penalty parameters.

An example from a wave tank:

.. code-block:: yaml

    forcing_zones:
    -   name: outlet velocity damping
        type: MomentumForcing
        zone: outlet zone/beta
        target: waves/u
        penalty: 10
        plot: no

.. describe:: name

    Used for log output and better error messages only

.. describe:: type

    One of ``MomentumForcing`` or ``ScalarForcing``

.. describe:: variable

    If you use ScalarForcing then you must give the name of the variable, e.g.
    ``c`` to force the colour field. For MomentumForcing this defaults to the
    velocity field ``u``.

.. describe:: penalty

    The penalty used to "nunge" the solution towards the target

.. describe:: zone

    The name of the known field function that is 1.0 inside the zone and 0.0
    outside. The transition is typically smooth. Often a
    :ref:`inp_fields_scalar` is used for this purpose.

.. describe:: target

    The name of a known field function that we want to "nunge" our solution
    towards using the forcing zone.

.. describe:: plot

    Show the forcing zone (plots to file). Default off.
