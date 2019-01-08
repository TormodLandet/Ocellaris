.. _inp_convection:

Convection
==========

Some numerical methods, most notaby the algebraic Volume of Fluid method for
interface capturing in free surface flows, require the selection of a
convective flux.


Upwind
------

The simplest flux for the colour function is pure upwind. This flux leads to
a very diffusive free surface when used with a piecewise constant VOF colour
function.

.. code-block:: yaml

    convection:
        c:
            convection_scheme: Upwind


HRIC
----

High Resolution Interface Capturing. A compressive flux designed to be used in
combination with a piecewise constant VOF colour function.

Three versions are implemented

* Standard HRIC (**HRIC**, default, the only one that is well tested)

  "A Two-Fluid Navie-Stokes Solver to Simulate Water Entry",
  Twenty-Second Symposium on Naval Hydrodynamics, 1999
  S. Muzaferija, M. Peric, P. Sames, t. Schellin

* Modified HRIC (**MHRIC**, practically untested code, use with caution)

  "Fluent Theory Guide v15.0"
  2013
  ANSYS

* Refined HRIC (**RHRIC**, practically untested code, use with caution)

  "Free surface flow analysis based on improved HRIC VOF method",
  Journal of the Society of Naval Architects of Korea,
  Vol. 47, No 3, pp. 279-290, June 2010
  Il-Ryong Park, Kwang-Soo Kim, Jin Kim and Suak-Ho Van

.. code-block:: yaml

    convection:
        c:
            convection_scheme: HRIC
            HRIC_version: HRIC          # <-- not needed, HRIC is default

CICSAM
------

Compressive Interface Capturing Scheme for Arbitrary Meshes. A compressive flux
designed to be used in combination with a piecewise constant VOF colour
function.

  "Numerical prediction of two fluid systems with sharp interfaces",
  Imperial College, London, 1997,
  Onno Ubbink

.. code-block:: yaml

    convection:
        c:
            convection_scheme: CICSAM
