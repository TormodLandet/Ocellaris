import dolfin


class OcellarisConstant(dolfin.Constant):
    def __init__(self, value):
        """
        This is an overloaded dolfin.Constant which can be asked for its value
        without having to supply coordinates and evaluating it as a function
        (it is constant after all)
        """
        super(OcellarisConstant, self).__init__(value)
        self.py_value = value
    
    @property
    def py_value(self):
        return self._py_value
    
    @py_value.setter
    def py_value(self, value):
        self.assign(dolfin.Constant(value))
        self._py_value = value
