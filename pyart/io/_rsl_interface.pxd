""" Shared extensions in _rsl_interface, used by _fourdd_interface. """


cimport _rsl_h

cdef class _RslVolume:
    cdef _rsl_h.Volume * _Volume
    cdef int _dealloc
    cdef load(self, _rsl_h.Volume * Volume)
    cdef _prtmode(self, _rsl_h.Ray_header h)
