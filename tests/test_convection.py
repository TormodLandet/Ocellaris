from dolfin import UnitSquareMesh, FunctionSpace, Function, Expression, as_vector
from ocellaris import Simulation
from ocellaris.solver_parts.convection import get_convection_scheme, VelocityDGT0Projector


def mk_scheme(N, Vname, Vorder, cpp_expr, expr_args, convection_inp):
    mesh = UnitSquareMesh(N, N)
    V = FunctionSpace(mesh, Vname, Vorder)
    C = Function(V)
    e = Expression(cpp_expr, element=V.ufl_element(), **expr_args)
    C.interpolate(e)
    
    D = Function(V)
    D.assign(C)
    
    sim = Simulation()
    sim.set_mesh(mesh)
    sim.data['constrained_domain'] = None
    sim.data['C'] = C
    for key, value in convection_inp.items():
        sim.input.set_value('convection/C/%s' % key, value)
    
    scheme_name = convection_inp['convection_scheme']
    return get_convection_scheme(scheme_name)(sim, 'C')


def mk_vel(sim, Vname, Vorder, cpp_exprs):
    mesh = sim.data['mesh']
    V = FunctionSpace(mesh, Vname, Vorder)
    
    vel = []
    for cpp in cpp_exprs:
        u = Function(V)
        u.interpolate(Expression(cpp, element=V.ufl_element()))
        vel.append(u)
    
    return as_vector(vel)


def mk_blending_factor(convection_inp):
    # Create convected function and convection scheme
    N = 10
    scheme = mk_scheme(N, Vname='DG', Vorder=0, 
                       cpp_expr='A + A*sin(B*pi*x[0])*sin(B*pi*x[1])',
                       expr_args={'A': 0.5, 'B': 2},
                       convection_inp=convection_inp)
    sim = scheme.simulation
    C = sim.data['C']
    beta = scheme.blending_function
    
    # Filter C to make it more like a free surface simulation
    Carr = C.vector().get_local()
    Carr[Carr > 0.8] = 1
    Carr[Carr < 0.2] = 0
    C.vector().set_local(Carr)
    
    # Make convecting velocity
    vel = mk_vel(sim, 'DG', 2, ['1', '1'])
    proj = VelocityDGT0Projector(sim, vel)
    proj.update()
    vel_dgt0 = proj.velocity
    
    # Run the convection scheme
    h = 1/N
    dt = h * 0.3 / 2**0.5
    if scheme.need_alpha_gradient:
        scheme.initialize_gradient()
        scheme.gradient_reconstructor.reconstruct()
    scheme.update(dt, vel_dgt0)
    
    if False:
        import dolfin
        from matplotlib import pyplot
        
        pyplot.figure()
        patches = dolfin.plot(C)
        pyplot.colorbar(patches)
        pyplot.savefig('test_conv_C.png')
        
        from ocellaris.utils.plotting_trace import plot_matplotlib_dgt
        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        patches = plot_matplotlib_dgt(beta)
        fig.colorbar(patches, ax=ax)
        pyplot.savefig('test_conv_beta.png')
    
    return beta, sim


def test_hric_cpp():
    """
    Run HRIC with and without C++ implementation to  verify the C++ 
    implementation correctness
    """
    # regression test value
    ANS1 = 7.85326907878949600
    ANS2 = 0.42426406871193789
    
    cinp = {'convection_scheme': 'HRIC', 
            'HRIC_version': 'HRIC'}
    
    betas = []
    Cofs = []
    for use_cpp in (False, True):
        cinp['use_cpp'] = use_cpp
        beta, sim = mk_blending_factor(cinp)
        betas.append(beta)
        bn = beta.vector().norm('l2')
        
        Cof_max = sim.reporting.timestep_xy_reports['Cof_max'][-1]
        Cofs.append(Cof_max)
        
        print(use_cpp, repr(bn), bn - ANS1, repr(Cof_max), Cof_max - ANS2)
        assert abs(bn - ANS1) < 1e-15
        assert abs(Cof_max - ANS2) < 1e-15
    
    diff = betas[0].copy(deepcopy=True)
    diff.vector().axpy(-1, betas[1].vector())
    
    assert diff.vector().norm('l2') == 0


def test_gradient_cpp():
    """
    Run the gradient reconstruction with and without C++ to verify the C++
    implementation correctness
    """
    N = 10
    
    schemes = []
    for use_cpp in (False, True):
        convection_inp = {'use_cpp_gradient': use_cpp,
                          'convection_scheme': 'HRIC'}
        scheme = mk_scheme(N, Vname='DG', Vorder=0, 
                           cpp_expr='A + A*sin(B*pi*x[0])*sin(B*pi*x[1])',
                           expr_args={'A': 0.5, 'B': 2},
                           convection_inp=convection_inp)
        scheme.initialize_gradient()
        scheme.gradient_reconstructor.reconstruct()
        schemes.append(scheme)
    
    g0_py = schemes[0].gradient_reconstructor.gradient[0].copy(deepcopy=True)
    g1_py = schemes[0].gradient_reconstructor.gradient[1].copy(deepcopy=True)
    g0_cpp = schemes[1].gradient_reconstructor.gradient[0].copy(deepcopy=True)
    g1_cpp = schemes[1].gradient_reconstructor.gradient[1].copy(deepcopy=True)
    
    for g0, g1 in ((g0_py, g0_cpp), (g1_py, g1_cpp)):
        diff = g0.copy(deepcopy=True)
        diff.vector().axpy(-1, g1.vector())
        g0n = g0.vector().norm('l2')
        assert 20 > g0n > 10
        assert diff.vector().norm('l2') < g0n/1e15
