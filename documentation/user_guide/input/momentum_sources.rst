.. inp_momentum_sources:

Momentum sources
================

You can add additional momentum sources to the simulations which can be useful
for MMS testing (Method of Manufactured Solutions). An example from the
variable density disk test case included in the Ocellaris source code
repository can be seen below (this MMS test was taken from
:cite:`guermond_salgado_2011`).

.. code-block:: yaml

    momentum_sources:
    -   name: MMS
        degree: 2
        cpp_code:
        -   Q*nu*sin(sin(t))*cos(t) - (x[0]*cos(t)*cos(t) - x[1]*sin(t))*(Q*sqrt(x[0]*x[0] + x[1]*x[1])*cos(sin(t) - atan2(x[1], x[0])) + 1 + Q) + sin(t)*sin(x[1])*cos(x[0])
        -   -Q*nu*cos(t)*cos(sin(t)) - (x[0]*sin(t) + x[1]*cos(t)*cos(t))*(Q*sqrt(x[0]*x[0] + x[1]*x[1])*cos(sin(t) - atan2(x[1], x[0])) + 1 + Q) + sin(t)*sin(x[0])*cos(x[1])

.. describe:: degree

    The polynomial degree of the function used to represent the momentum source
    C++ expression
