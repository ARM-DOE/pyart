========================================
Radar Corrections (:mod:`pyart.correct`)
========================================

.. currentmodule:: pyart.correct

Py-ART has modules, classes and functions which are able to correct a 
number of common problems with radar data.

Field Corrections
=================

.. autosummary::
    :toctree: generated/

    dealias.dealias
    attenuation.calculate_attenuation
    phase_proc.phase_proc
    
Phase Processing
================
.. autosummary::
    :toctree: generated/

    phase_proc.det_sys_phase
    phase_proc.det_process_range
    phase_proc.fzl_index
    phase_proc.snr
    phase_proc.noise
    phase_proc.sobel
    phase_proc.unwrap_masked
    phase_proc.smooth_and_trim
    phase_proc.get_phidp_unf
    phase_proc.construct_A_matrix
    phase_proc.construct_B_vectors
    phase_proc.LP_solver

Utilities
=========

.. autosummary::
    :toctree: generated/

    dealias.find_time_in_interp_sonde

