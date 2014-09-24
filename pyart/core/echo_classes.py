"""
pyart.core.echo_classes
===============

A Module where enumerated echo classifications are held

.. autosummary::
    :toctree: generated/
    :template: dev_template.rst

    echo_classes

"""

CLASS_NOISE = 0 #returs for which there is no significant scattering
CLASS_CLUTTER = 1 #Stationary persistant clutter.
CLASS_NONMETEO = 2 # Non-stationary but non-weather related, eg birds/bragg etc
CLASS_MULTITRIP = 3
CLASS_LIQUID = 4
CLASS_ICE = 5
CLASS_MIXEDPHASE = 6
CLASS_HAIL = 7


#This string will be helpful when saving files etc.
CLASS_METADATA = """Numerical classifications for dominant scatterer at gate:
0: Returs for which there is no significant scattering
1: Stationary persistant clutter.
2: Non-stationary but non-weather related, eg birds/bragg etc
3: Multi-trip echoes beyond the MUR
4: Liquid precipitation
5: Solid (ice) precipitation
6: Mixed ice and liquiq
7: Dominant hail (ie hail contaminated)
"""
