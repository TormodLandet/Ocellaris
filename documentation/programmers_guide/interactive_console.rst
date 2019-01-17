
.. _sec-interactive-console:

Interactive console
===================

At the end of each time step Ocellaris will optionally open an interactive
console so that you can inspect the internal state of the simulation. To
access this pres :kbd:`d` then :kbd:`Enter` ("d" for debug). At the end of the
next time step the console should open and you will have full access to the
internal variables. The variables are listed so that you can get a head start.

.. note:: This only works when running interactively on 1 CPU

Most of the variables are described in the :ref:`simulation-api` documentation
under the :attr:`ocellaris.Simulation.data` attribute.

If you press :kbd:`Ctrl+d` inside the interactive console Ocellaris will
continue running the time loop. If you type ``exit()`` or  ``quit()`` you will
stop Ocellaris and return to the command line immediately.

It is also possible to specify that the console should open at the end of the
simulation. If you want this put the following on the input file:

.. code-block:: yaml

    console_at_end: true

This can be very useful for ad-hoc postprocessing of the simulation results.
