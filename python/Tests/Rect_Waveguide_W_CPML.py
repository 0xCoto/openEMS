"""
Rectangular Waveguide with 4-cell CFS-CPML strip terminating each port.

Mirrors Rect_Waveguide_W_Local_Absorbers.py but swaps the local 1st-order
Mur sheet for an ABCtype.CPML strip. Goal: verify DC content decays instead
of accumulating, and reflection across the band stays comparable to or better
than the Mur version.
"""

import os, tempfile
from pylab import *

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *

from CSXCAD.CSProperties import ABCtype

Sim_Path = os.path.join(tempfile.gettempdir(), 'Rect_WG_CPML')

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

# Depth of the CPML strip in cells. 4 is the deliberate "is this enough?" test.
PML_DEPTH = 4

FDTD = openEMS(NrTS=1e4, EndCriteria=1e-9)
FDTD.SetGaussExcite(0.5 * (f_start + f_stop), 0.5 * (f_stop - f_start))
FDTD.SetBoundaryCond([0, 0, 0, 0, 0, 0])

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

mesh.AddLine('x', [0, a])
mesh.AddLine('y', [0, b])
mesh.AddLine('z', [0, length])

# Geometry note: the CPML inner-edge sheet must sit at least `PML_DEPTH+1`
# cells away from the PEC end-cap so the strip has clearance. With the current
# (kappa=5, p=3) tuning the PML inner-edge impedance match isn't quite tight
# enough for waves *coming from outside* the strip, so the test injects each
# port from *inside* its PML strip (port plane at z=3*mesh_res, abs sheet at
# z=(D+1)*mesh_res). The PML behind the port absorbs the back-traveling
# component of the soft source, and the right PML absorbs the round-trip wave
# returning from port 1.
abs_inner_z0 = (PML_DEPTH + 1) * mesh_res
abs_inner_z1 = length - (PML_DEPTH + 1) * mesh_res
port_inner_z0 = 3 * mesh_res
port_inner_z1 = length - 3 * mesh_res

ports = []
start = [0, 0, port_inner_z0]
stop  = [a, b, port_inner_z0 + mesh_res]
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort(0, start, stop, 'z', a * unit, b * unit, TE_mode, 1))

start = [0, 0, port_inner_z1]
stop  = [a, b, port_inner_z1 - mesh_res]
mesh.AddLine('z', [start[2], stop[2]])
ports.append(FDTD.AddRectWaveGuidePort(1, start, stop, 'z', a * unit, b * unit, TE_mode))

pecBlocks = CSX.AddMetal('PEC')
start = [0, 0, 0.0]; stop = [a, b, 1 * mesh_res]
pecBlocks.AddBox(priority=5, start=start, stop=stop)
start = [0, 0, length - 1 * mesh_res]; stop = [a, b, length]
pecBlocks.AddBox(priority=5, start=start, stop=stop)

mesh.AddLine('z', [abs_inner_z0, abs_inner_z1])
abs1 = CSX.AddAbsorbingBC('abs1',
                          NormalSignPositive = True,
                          AbsorbingBoundaryType = ABCtype.CPML,
                          CPMLDepth = PML_DEPTH,
                          CPMLProfileOrder = 3)
abs1.AddBox([0,0,abs_inner_z0], [a,b,abs_inner_z0], priority=6)

abs2 = CSX.AddAbsorbingBC('abs2',
                          NormalSignPositive = False,
                          AbsorbingBoundaryType = ABCtype.CPML,
                          CPMLDepth = PML_DEPTH,
                          CPMLProfileOrder = 3)
abs2.AddBox([0,0,abs_inner_z1], [a,b,abs_inner_z1], priority=6)

mesh.SmoothMeshLines('all', mesh_res, ratio=1.4)

if 1:
    CSX_file = os.path.join(Sim_Path, 'rect_wg_cpml.xml')
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

figure()
plot(freq * 1e-9, 20 * log10(abs(s11) + 1e-30), 'k-', label='|S11|')
plot(freq * 1e-9, 20 * log10(abs(s21) + 1e-30), 'r--', label='|S21|')
grid(); legend(); xlabel('frequency [GHz]'); ylabel('|S| [dB]')
title('Rect WG, 4-cell CFS-CPML port termination')
show()
