#cython: c_string_encoding=utf8, c_string_type=unicode, language_level=3, nonecheck=True
# distutils: language=c++

## Imports.

# Cython operator workarounds.
from cython.operator cimport address as addr       # &foo
from cython.operator cimport dereference as deref  # *foo
from cython.operator cimport preincrement as inc   # ++foo

## my imports
from libcpp.vector cimport vector
from .pywrapfst cimport *
from .cops cimport *
import numpy as np

def runme():
    print('This works')

cpdef MutableFst align(Fst model, Fst ifst1, Fst ifst2):
  
  cdef unique_ptr[fst.VectorFstClass] tfst
  tfst.reset(new fst.VectorFstClass(ifst1.arc_type()))

  cdef fst.WeightClass distance
  shortest_path(deref(model._fst), deref(ifst1._fst), deref(ifst2._fst), tfst.get())

  return _init_MutableFst(tfst.release())

cpdef Weight score_std(Fst model, Fst ifst1, Fst ifst2):
  distance = shortest_distance_std(deref(model._fst), deref(ifst1._fst), deref(ifst2._fst))
  retval = Weight(model._fst.get().WeightType(), distance.ToString())
  return retval


cpdef Weight score_log(Fst model, Fst ifst1, Fst ifst2):
  distance = shortest_distance_log(deref(model._fst), deref(ifst1._fst), deref(ifst2._fst))
  retval = Weight(model._fst.get().WeightType(), distance.ToString())
  return retval

cpdef Weight kernel_score_std(Fst model, Fst ifst1, Fst ifst2):
  distance = kernel_score_std_impl(deref(model._fst), deref(ifst1._fst), deref(ifst2._fst))
  retval = Weight(model._fst.get().WeightType(), distance.ToString())
  return retval


