# Experimental reactive Voronoi acceleration

The three scalar actions retain their existing mathematical definitions.
CPU evaluation automatically shares compatible geometry, softmax weights and
occupancies within one PLUMED instance. NO_SHARE selects the independent
reference path. No process-global cache is used.

Compatibility requires identical ordered CENTERS/ASSIGNED, PBC and SERIAL
settings, neighbor cutoff/skin/stride and neighbor algorithm. Different KAPPA
values share CPU geometry but keep separate normalized weights. REFERENCE and
all scalar reduction parameters remain consumer-local. Compatible strided
neighbor lists share candidate-list ownership and displacement references.
Exact coordinates, box and step are checked on every evaluation request.
Coordinate changes at a repeated MD step invalidate the cache.

Each action still computes its own analytic derivatives. PLUMED combines bias
sensitivities through its usual force graph. Direct distance/position terms and
box derivatives are not replaced. Sharing the forward pass is not equivalent
to fusing all reverse passes.

## Optional CUDA backend

Build on a Linux CUDA host, selecting the architecture of the intended GPU:

~~~bash
nvcc -O3 -std=c++17 -arch=sm_70 -Xcompiler=-fPIC -shared \
  ReactiveVoronoiCuda.cu -o libReactiveVoronoiCuda.so
~~~

CPU-only Linux installations do not require CUDA headers or the optional library.
GPU selection is explicit on each compatible action:

~~~plumed
UNITS LENGTH=A
q: VORONOI_COORDINATION CENTERS=1,2 ASSIGNED=3,4 KAPPA=5 REFERENCE=1 POWER=2 BACKEND=CUDA CUDA_LIBRARY=./libReactiveVoronoiCuda.so
d: VORONOI_DISTANCE CENTERS=1,2 ASSIGNED=3,4 KAPPA=5 REFERENCE=1 GROUP1=1 GROUP2=2 BACKEND=CUDA CUDA_LIBRARY=./libReactiveVoronoiCuda.so
PRINT ARG=q,d FILE=COLVAR
~~~

The prototype supports complete exponential assignments, a single MPI rank,
orthorhombic PBC and NOPBC. It rejects NLIST, triclinic PBC and multi-rank CUDA
requests explicitly; these inputs remain supported by the CPU backend.
The regional compact-kernel action is a separate implementation, not
implicitly converted to the exponential kernel.

FP64 forward and reverse kernels stream over centers for each assigned atom.
They retain O(Ncenters+Nassigned) device data rather than a dense weight matrix.
One context is reused per compatible KAPPA. Coordinate, occupancy, sensitivity,
force and virial transfers use an owned CUDA stream. OPES, scalar reductions
and direct terms remain on the CPU. The plugin uses the currently selected
CUDA device and respects the host CUDA_VISIBLE_DEVICES mapping.

The current LAMMPS-PLUMED host interface still supplies host coordinates.
This implementation does not claim device-resident MD integration or overlap
with DeePMD execution. GPU kernel throughput and complete MD throughput must
be measured separately. Small systems may be slower on CUDA.

## Validation

The rt-reactive-voronoi-shared test checks simultaneous biases, mixed KAPPA,
reordered selections, repeated-step coordinate and box changes, independent
PLUMED instances, strided-list behavior and finite-difference forces. Its
executable optionally takes a CUDA keyword suffix to compare GPU evaluation
against the independent CPU path.

Use native driver trajectories to compare CVs and every coordinate/box
derivative, then compare fixed-configuration total bias forces. Long chaotic
trajectories are not a pointwise equivalence criterion. Match hardware,
timestep, output and OPES physical update intervals in MD benchmarks.
Keep a nonempty, matching OPES STATE with each MD checkpoint.

A host MD engine must accept the actual PLUMED command API of the build.
For an older host with a strict upper API guard, backport these CV changes to
a qualified PLUMED base and repeat interface/driver/force tests. Do not spoof
the command API or change a shared installation to bypass that guard.
