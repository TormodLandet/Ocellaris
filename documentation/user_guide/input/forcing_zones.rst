.. _inp_forcing_zones:

Forcing zones
=============

Forcing zones can be used to damp waves near the inlet and outlet to avoid
reflections. They work by specifying a zone (a known field), a target result
(a known field) and a penalty to be added to the equation system in order to
pull the unknown solutiuon towards the known field solution in the forcing
zone.
