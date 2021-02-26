#cython: language_level=3

from libcpp.vector cimport vector

cimport cpywrapfst as fst

cdef extern from "cops.h" nogil:
	cdef void shortest_path(const fst.FstClass &, const fst.FstClass &, const fst.FstClass &, fst.MutableFstClass *)
	cdef fst.WeightClass shortest_distance_std(const fst.FstClass &, const fst.FstClass &, const fst.FstClass &)
	cdef fst.WeightClass shortest_distance_log(const fst.FstClass &, const fst.FstClass &, const fst.FstClass &)
	cdef fst.WeightClass kernel_score_std_impl(const fst.FstClass &, const fst.FstClass &, const fst.FstClass &)