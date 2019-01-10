.. _inp_hooks:

Hooks - run custom code at given times during a simulation
==========================================================

Hooks allows you to run code at given times during a simulation

* ``pre_simulation``
* ``post_simulation``
* ``pre_timestep``
* ``post_timestep``

.. * Custom hooks, most notably ``MultiPhaseModelUpdated``.

Each hook is described be the following

.. describe:: name

    The name of the hook, used for better log and error messages

.. describe:: enabled

    Flag to turn the hook on or off. Default on

.. describe:: code

    The Python code to run. You can access the ``simulation`` object which
    gives access to all fields and every part of the simulation that the normal
    Ocallaris code can access. You also have all user code constants instantly
    available just like all other code in the input file.

The example below shows that each hook gets it's own dictionary ``hook_data``
to store whatever it wants between calls. The example also shows how to read
the input file parameters in a hook that is defined in the same input file, and
how to perform output to file in a configurable manner:

.. code-block:: yaml

    hooks:
        post_timestep:
        -   name: save colour function field
            enabled: yes
            code: |
                if not 'cf' in hook_data:
                    prefix = simulation.input.get_value('output/prefix')
                    hook_data['cf'] = File(prefix + '_c.pvd')
                if t > 1:
                    hook_data['cf'] << (c, t)

The ``hook_data`` dictionary is saved to restart files and the contents are
brought back as long as it consists of basic data types (lists, dicts, strings,
numbers) since the data is internally serialized to ASCII YAML format before
saving to HDF5 (which is the binary format of the restart file).
