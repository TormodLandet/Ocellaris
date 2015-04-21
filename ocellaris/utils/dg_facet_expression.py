import dolfin

def FacetExpressionDG0(mesh):
    """
    Create a facet function that can be evaluated on facets in 
    ds and dS integrals
    """
    f = dolfin.Expression(code)
    f.mesh = mesh
    f.facet_data = dolfin.FacetFunctionDouble(mesh)
    f.ocellaris_cpp_expression_type = 'FacetExpressionDG0'
    return f

code = '''
#include <dolfin/log/Event.h>

class FacetExpressionDG0 : public Expression
{
public:

    std::shared_ptr<MeshFunction<double> > facet_data;
    std::shared_ptr<const Mesh> mesh;

    // Constructor
    FacetExpressionDG0()
    :
       Expression(),
       not_on_boundary("*** Warning: evaluating FacetExpressionDG0 on a non-facet domain, returning zero.")
    {
    }
    
    void eval(Array<double>& values, const Array<double>& x,
              const ufc::cell& ufc_cell) const
    {
        int ndim = (*mesh).geometry().dim();
        assert(facet_data);
        assert(ufc_cell.geometric_dimension == ndim);

        if (ufc_cell.local_facet >= 0)
        {
            // We are evaluated on a facet
            // Get the facet index of the current facet 
            const Cell dolfin_cell(*mesh, ufc_cell.index);
            int fidx = dolfin_cell.entities(ndim-1)[ufc_cell.local_facet];
            
            //cout << "at ufc_cell.index " << ufc_cell.index
            //     << " and local facet " << ufc_cell.local_facet
            //     << " and facet idx " << fidx << endl;
            
            // Return the value provided in the facet_data mesh function
            values[0] = (*facet_data)[fidx];
        }
        else
        {
            // We are evaluated in a cell (dx integral), not on a facet
            // (ds or dS integral). This should not happen. Report and
            // return 0 
            not_on_boundary();
            values[0] = 0.0;
        }
  }

private:
    // Warning when evaluating on cells
    mutable Event not_on_boundary;
};
'''

if __name__ == '__main__':
    mesh = dolfin.UnitSquareMesh(3, 3, 'left/right')
    V = dolfin.FunctionSpace(mesh, 'CG', 1)
    
    x = mesh.coordinates()
    x[:,0] = x[:,0]**2
    x[:,1] = x[:,1]**3

    func = FacetExpressionDG0(mesh)
    func2 = dolfin.FacetArea(mesh)

    # Make func into a facet area function
    for cell in dolfin.cells(mesh):
        for i, facet in enumerate(dolfin.facets(cell)):
            func.facet_data[facet.index()] = cell.facet_area(i)

    # Check that the functions are the same
    v = dolfin.TestFunction(V)
    a = dolfin.assemble(func('+')*v('+')*dolfin.dS).array()
    a2 = dolfin.assemble(func2('+')*v('+')*dolfin.dS).array()
    print a
    print a2
    assert (abs(a - a2) < 1e-10).all()
