from ocellaris.solver_parts.fields.airy_waves import get_airy_wave_specs


def test_get_airy_wave_specs():
    # Round all numbers in an iterable and return a tuple
    rt = lambda q: tuple(round(v, 7) for v in q)
    # Run rt() on all items in an iterable and return a tuple
    tt = lambda q: tuple(rt(v) for v in q)

    omegas = [0.2, 0.0209439510239, 1.7340744614782326, 17.1551741466]
    periods = [31.4159265359, 300, 3.623365343736975, 0.366255990962]
    wave_lengths = [259.096496708, 2485.89122837, 20, 0.209439510239]
    wave_numbers = [0.02425036767, 0.00252753830717, 0.3141592653589793, 30]
    d1 = tt(get_airy_wave_specs(9.81, 7.0, omegas=omegas))
    d2 = tt(get_airy_wave_specs(9.81, 7.0, periods=periods))
    d3 = tt(get_airy_wave_specs(9.81, 7.0, wave_lengths=wave_lengths))
    d4 = tt(get_airy_wave_specs(9.81, 7.0, wave_numbers=wave_numbers))
    assert d1 == d2 == d3 == d4
