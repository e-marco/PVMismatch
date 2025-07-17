# -*- coding: utf-8 -*-
"""
This is the PVMismatch Package. It contains :mod:`~pvmismatch.pvmismatch_lib`
and :mod:`~pvmismatch.pvmismatch_tk`.

:mod:`~pvmismatch.pvmismatch_lib`
=================================
This package contains the basic library modules, methods, classes and
attributes to model PV system mismatch.

.. note::
   The main library classes and modules are exposed through this package for
   convenience.

   For example::

       >>> from pvmismatch import PVcell  # imports the PVcell class
       >>> # import pvconstants, pvcell, pvmodule, pvstring and pvsystem
       >>> from pvmismatch import *

:mod:`~pvmismatch.pvmismatch_tk`
================================
This package contains an application that can be run using
:mod:`pvmismatch.pv_tk`.
"""

# Import the version from the package metadata.
try:
    # Python 3.8+
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    # Python <3.8
    from importlib_metadata import version, PackageNotFoundError

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "0+unknown"

from .pvmismatch_lib import (
    pvconstants, pvcell, pvcell_from_data,
    pvmodule, pvstring, pvsystem, pvexceptions,
)
PVconstants = pvconstants.PVconstants
PVcell      = pvcell.PVcell
PVmodule    = pvmodule.PVmodule
PVstring    = pvstring.PVstring
PVsystem    = pvsystem.PVsystem

__all__ = [
    "PVconstants", "PVcell", "PVmodule", "PVstring", "PVsystem",
    "pvconstants", "pvcell", "pvcell_from_data",
    "pvmodule", "pvstring", "pvsystem", "pvexceptions",
    "__version__",
]