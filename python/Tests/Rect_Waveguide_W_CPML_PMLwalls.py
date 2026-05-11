"""
Rectangular Waveguide with 4-cell CFS-CPML strip — PML_8 outer boundary,
PEC waveguide walls as explicit blocks. Production-shape configuration
that the memory note flags as triggering CPML positive-feedback divergence.

Goal: reproduce the divergence so we can capture |Ψ|max trace through the
unstable mode and bisect which buffer / sign / coefficient drives it.
"""

import os, sys, tempfile
import numpy as np
from pylab import *

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

from CSXCAD.CSProperties import ABCtype

Sim_Path = os.path.join(tempfile.gettempdir(), 'Rect_WG_CPML_PMLwalls')

post_proc_only = False
unit = 1e-6

# WR-WideBand-ish (closer to production manifest: a=100mm, b=50mm, 1.5-3 GHz)
a = 100e3       # 100 mm
b = 50e3        # 50 mm
length = 200e3  # 200 mm

f_start = 1.5e9
f_0     = 2.25e9
f_stop  = 3.0e9
lambda0 = C0 / f_0 / unit

TE_mode = 'TE10'
mesh_res = lambda0 / 60  # ~2.2 mm at 2.25 GHz, gets cell count near production

PML_DEPTH = 4  # CPML port absorber depth
PML_OUTER = 8  # outer simulation boundary PML cells

# Production CPML params (R0=1e-3, α=0). Computed to match production
# exactly; in production these come out via the same formula in ports.py.
import math
_p_order = 3
_eta_0 = 376.730313461
_pml_thickness_m = PML_DEPTH * mesh_res * unit
_sigma_max = -(_p_order + 1) * math.log(1e-3) / (2.0 * _eta_0 * _pml_thickness_m)

# Padding around WG: λ/4 on x, y; the WG length itself is the z extent.
pad_xy = lambda0 / 4

FDTD = openEMS(NrTS=int(1e6), EndCriteria=1e-4)  # production: NrTS=1e6, end_criteria=-40 dB
FDTD.SetGaussExcite(0.5 * (f_start + f_stop), 0.5 * (f_stop - f_start))
FDTD.SetBoundaryCond([f'PML_{PML_OUTER}'] * 6)

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

# Mesh extends past WG by pad_xy in x, y and just covers WG in z.
mesh.AddLine('x', [-pad_xy, 0, a, a + pad_xy])
mesh.AddLine('y', [-pad_xy, 0, b, b + pad_xy])
mesh.AddLine('z', [0, length])

# Port plane and absorber/PEC backing geometry — same convention as production.
port_start_z0 = PML_DEPTH * mesh_res
port_start_z1 = length - PML_DEPTH * mesh_res

abs_inner_z0  = port_start_z0 + mesh_res
abs_inner_z1  = port_start_z1 - mesh_res

wg_port_z0    = port_start_z0 + 2 * mesh_res
wg_port_z1    = port_start_z1 - 2 * mesh_res

ports = []
start = [0, 0, wg_port_z0]
stop  = [a, b, wg_port_z0 + mesh_res]
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort(0, start, stop, 'z', a * unit, b * unit, TE_mode, 1))

start = [0, 0, wg_port_z1]
stop  = [a, b, wg_port_z1 - mesh_res]
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort(1, start, stop, 'z', a * unit, b * unit, TE_mode))

# PEC waveguide side walls (4 sheets), one cell thick each, around the WG cavity.
wg_walls = CSX.AddMetal('WG_PEC_walls')
wg_walls.AddBox(priority=4, start=[-pad_xy, -pad_xy, 0],
                            stop=[0,        b + pad_xy, length])
wg_walls.AddBox(priority=4, start=[a,       -pad_xy, 0],
                            stop=[a + pad_xy, b + pad_xy, length])
wg_walls.AddBox(priority=4, start=[-pad_xy, -pad_xy, 0],
                            stop=[a + pad_xy, 0, length])
wg_walls.AddBox(priority=4, start=[-pad_xy, b, 0],
                            stop=[a + pad_xy, b + pad_xy, length])

# PEC backing block at each end.
pecBlocks = CSX.AddMetal('PEC_back')
pecBlocks.AddBox(priority=5, start=[0, 0, 0.0],
                             stop=[a, b, 1 * mesh_res])
pecBlocks.AddBox(priority=5, start=[0, 0, length - 1 * mesh_res],
                             stop=[a, b, length])

# CPML strips.
mesh.AddLine('z', [abs_inner_z0, abs_inner_z1])
abs1 = CSX.AddAbsorbingBC('abs1',
                          NormalSignPositive=True,
                          AbsorbingBoundaryType=ABCtype.CPML,
                          CPMLDepth=PML_DEPTH,
                          CPMLAlphaMax=0.0,
                          CPMLSigmaMax=_sigma_max,
                          CPMLProfileOrder=_p_order)
abs1.AddBox([0, 0, abs_inner_z0], [a, b, abs_inner_z0], priority=6)

abs2 = CSX.AddAbsorbingBC('abs2',
                          NormalSignPositive=False,
                          AbsorbingBoundaryType=ABCtype.CPML,
                          CPMLDepth=PML_DEPTH,
                          CPMLAlphaMax=0.0,
                          CPMLSigmaMax=_sigma_max,
                          CPMLProfileOrder=_p_order)
abs2.AddBox([0, 0, abs_inner_z1], [a, b, abs_inner_z1], priority=6)

mesh.SmoothMeshLines('all', mesh_res, ratio=1.4)

if 1:
    CSX_file = os.path.join(Sim_Path, 'rect_wg_cpml_pmlwalls.xml')
    if not os.path.exists(Sim_Path):
        os.mkdir(Sim_Path)
    CSX.Write2XML(CSX_file)

if not post_proc_only:
    FDTD.Run(Sim_Path, verbose=3, cleanup=True)

freq = linspace(f_start, f_stop, 201)
for port in ports:
    port.CalcPort(Sim_Path, freq)

s11 = ports[0].uf_ref / ports[0].uf_inc
s21 = ports[1].uf_ref / ports[0].uf_inc

if "--no-show" not in sys.argv:
    figure()
    plot(freq * 1e-9, 20 * log10(abs(s11) + 1e-30), 'k-', label='|S11|')
    plot(freq * 1e-9, 20 * log10(abs(s21) + 1e-30), 'r--', label='|S21|')
    grid(); legend(); xlabel('frequency [GHz]'); ylabel('|S| [dB]')
    title('Rect WG, PML_8 outer + 4-cell CPML port termination')
    show()
