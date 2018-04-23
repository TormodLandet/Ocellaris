import dolfin


def check_vector_value_histogram(vec, expected, comm=None, round_digits=0):
    """
    Given a dict of expected values (keys) and the number of this value expected
    to be in the vector (value), check if the histogram (value count) of the
    vector vec corresponds to the expected or not
    """
    if comm is None:
        comm = dolfin.MPI.comm_world
    
    hist = get_vector_value_histogram(vec, comm)
    print(hist, expected)
    assert len(hist) == len(expected)
    for k, v in hist.items():
        if round_digits:
            k = round(k, round_digits)
        assert k in expected
        assert v == expected[k]


def get_vector_value_histogram(vec, comm=None):
    """
    Return a dict with keys corresponding to the values of the vector vec and
    where the dictionary values is the corresponding number of occurances
    """
    if comm is None:
        comm = dolfin.MPI.comm_world
    
    hist = {}
    for v in vec.get_local():
        if v in hist:
            hist[v] += 1
        else:
            hist[v] = 1
    
    all_hist = comm.gather(hist)
    hist2 = None
    if comm.rank == 0:
        hist2 = {}
        for hist in all_hist:
            for k, v in hist.items():
                if k in hist2:
                    hist2[k] += v
                else:
                    hist2[k] = v
    return comm.bcast(hist2)
