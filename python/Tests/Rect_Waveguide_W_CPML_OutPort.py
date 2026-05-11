"""
Rectangular Waveguide with 4-cell CFS-CPML strip — port OUTSIDE the strip.

Mirrors the production ports.py geometry (port plane on the simulation-interior
side of the absorber strip, PEC backing one PML-depth behind the strip).

Goal: reproduce the late-time positive-feedback divergence the standalone
PEC-walled test (port-INSIDE-PML) does not exhibit.
"""

import os, sys, tempfile
from pylab import *

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

from CSXCAD.CSProperties import ABCtype

Sim_Path = os.path.join(tempfile.gettempdir(), 'Rect_WG_CPML_OutPort')

post_proc_only = False
unit = 1e-6

a = 10700
b = 4300
length = 50000

f_start = 20e9
f_0     = 24e9
f_stop  = 26e9
lambda0 = C0 / f_0 / unit

TE_mode = 'TE10'
mesh_res = lambda0 / 50

PML_DEPTH = 4

# Layout (direction = +1, D = PML_DEPTH):
#   z =  0                           sim domain start (PEC boundary)
#   z =  0 .. 1*mesh_res             PEC backing block (1 cell thick)
#   z =  1*mesh_res .. (D+1)*mesh    CPML strip (k=0 high-σ at PEC end, k=D-1 low-σ at sheet)
#   z = (D+1)*mesh_res               PML inner edge / absorber sheet
#   z = (D+2)*mesh_res               wg_port plane (port_start + 2*mesh)
#                                    interior simulation volume
#   ... mirror on the other end ...

# port_start is at (D)*mesh_res (one cell past the PEC backing). The wg_port
# itself sits at port_start + 2*mesh_res, and the absorber sheet at
# port_start + 1*mesh_res. The PML strip extends D-1 more cells back.
port_start_z0 = PML_DEPTH * mesh_res
port_start_z1 = length - PML_DEPTH * mesh_res

abs_inner_z0  = port_start_z0 + mesh_res
abs_inner_z1  = port_start_z1 - mesh_res

wg_port_z0    = port_start_z0 + 2 * mesh_res
wg_port_z1    = port_start_z1 - 2 * mesh_res

FDTD = openEMS(NrTS=2e4, EndCriteria=1e-9)
FDTD.SetGaussExcite(0.5 * (f_start + f_stop), 0.5 * (f_stop - f_start))
# All-PEC outer boundary mirrors standalone test; bug is in the CPML
# kernel itself, so the simplest non-PML outer config still triggers it
# when the port-vs-strip topology matches production.
FDTD.SetBoundaryCond([0, 0, 0, 0, 0, 0])

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

mesh.AddLine('x', [0, a])
mesh.AddLine('y', [0, b])
mesh.AddLine('z', [0, length])

ports = []
# Port 0 (z = +): wg_port plane just past the strip on the interior side.
start = [0, 0, wg_port_z0]
stop  = [a, b, wg_port_z0 + mesh_res]
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort(0, start, stop, 'z', a * unit, b * unit, TE_mode, 1))

# Port 1 (z = -): mirrored.
start = [0, 0, wg_port_z1]
stop  = [a, b, wg_port_z1 - mesh_res]
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort(1, start, stop, 'z', a * unit, b * unit, TE_mode))

# PEC backing blocks (1 cell thick), placed PML_DEPTH cells behind each
# port plane (= immediately at the simulation wall in this layout).
pecBlocks = CSX.AddMetal('PEC_back')
pecBlocks.AddBox(priority=5,
                 start=[0, 0, 0.0],
                 stop=[a, b, 1 * mesh_res])
pecBlocks.AddBox(priority=5,
                 start=[0, 0, length - 1 * mesh_res],
                 stop=[a, b, length])

# CPML strips. NormalSignPositive=True for port 0 means the strip absorbs
# waves traveling in -z (toward the PEC). The sheet sits at abs_inner_z0,
# the strip extends D cells toward the PEC backing.
mesh.AddLine('z', [abs_inner_z0, abs_inner_z1])
abs1 = CSX.AddAbsorbingBC('abs1',
                          NormalSignPositive=True,
                          AbsorbingBoundaryType=ABCtype.CPML,
                          CPMLDepth=PML_DEPTH,
                          CPMLProfileOrder=3)
abs1.AddBox([0, 0, abs_inner_z0], [a, b, abs_inner_z0], priority=6)

abs2 = CSX.AddAbsorbingBC('abs2',
                          NormalSignPositive=False,
                          AbsorbingBoundaryType=ABCtype.CPML,
                          CPMLDepth=PML_DEPTH,
                          CPMLProfileOrder=3)
abs2.AddBox([0, 0, abs_inner_z1], [a, b, abs_inner_z1], priority=6)

mesh.SmoothMeshLines('all', mesh_res, ratio=1.4)

if 1:
    CSX_file = os.path.join(Sim_Path, 'rect_wg_cpml_outport.xml')
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
    title('Rect WG, port-OUTSIDE 4-cell CFS-CPML')
    show()
