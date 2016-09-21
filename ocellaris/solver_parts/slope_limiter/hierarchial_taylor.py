# enconding: utf8
from __future__ import division
import numpy
import dolfin as df
from ocellaris.cpp import load_module
from ocellaris.utils import ocellaris_error, verify_key
from . import register_slope_limiter, SlopeLimiterBase


@register_slope_limiter('HierarchalTaylor')
class NaiveNodalSlopeLimiter(SlopeLimiterBase):
    description = 'Uses a Taylor DG decomposition to limit derivatives at the vertices in a hierarchal manner'
    
    def __init__(self, phi_name, phi, boundary_condition, filter_method='nofilter', use_cpp=True):
        """
        Limit the slope of the given scalar to obtain boundedness
        """
    
    def run(self):
        """
        Perform slope limiting
        """
