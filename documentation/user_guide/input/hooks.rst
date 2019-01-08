.. _inp_hooks:

Hooks - run custom code at given times during a simulation
==========================================================

Hooks allows you to run code at given times during a simulation

* ``pre_simulation``
* ``post_simulation``
* ``pre_timestep``
* ``post_timestep``
* Custom hooks, most notably ``MultiPhaseModelUpdated``.

.. csv-table::
   :header: "key", "Default value", "Description"

    "...", "**required input**", "FIXME: finish this table"


Maintaining state with hook_data
--------------------------------

The example below shows that each hook gets it's own dictionary ``hook_data``
to store whatever it wants between calls. The example also shows how to read
the input file parameters in a hook that is defined in the same input file, and
how to perform output to file in a configurable manner:

.. code-block:: yaml

    -   name: save colour function field
        enabled: yes
        code: |
            if not 'cf' in hook_data:
                prefix = simulation.input.get_value('output/prefix')
                hook_data['cf'] = File(prefix + '_c.pvd')
            if t > 1:
                hook_data['cf'] << (c, t)
