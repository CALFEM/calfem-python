# -*- coding: utf-8 -*-
"""
CALFEM Visualisation module

Alias for calfem.vis_mpl, the matplotlib based visualisation module::

    import calfem.vis as cfv     # same as: import calfem.vis_mpl as cfv

calfem.vis was previously based on visvis. That module is deprecated and
available as calfem.vis_visvis.
"""

import sys

import calfem.vis_mpl as _vis_mpl

# Make calfem.vis the very same module object as calfem.vis_mpl, so that
# module state (current figure, colorbar mappable, ...) is shared.
sys.modules[__name__] = _vis_mpl
