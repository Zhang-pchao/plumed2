# Upstream preparation

The CPU-sharing change, the triclinic VORONOI_POSITION virial correction and
the optional CUDA backend can be reviewed as separate changes.

## Existing PLUMED GPU interfaces

The official repository provides cudaCoord and openaccPTM plugins. Their
device selection honors the MD-assigned device, otherwise captures the MD's
current device, and restores the caller's context around GPU work.

- https://github.com/plumed/plumed2/tree/master/plugins/cudaCoord
- https://github.com/plumed/plumed2/tree/master/plugins/openaccPTM
- https://github.com/plumed/plumed2/blob/master/.github/CONTRIBUTING.md

The existing COORDINATION CUDA plugin evaluates pair switching functions.
Reactive Voronoi assignments normalize each assigned atom over all candidate
centers and require a coupled reverse derivative. This backend preserves
that normalization and the per-action scalar interfaces.

BACKEND is the current experimental Voronoi interface. A registered plugin
using LOAD and USEGPU, or integration with the task-manager interfaces,
requires an upstream API and maintenance discussion. The current core keeps
CUDA headers and runtime linking confined to the optional library.

## Source and maintenance checks

Core CPU regression coverage is in rt-reactive-voronoi-backends and
rt-reactive-voronoi-shared. The backend test compares the three actions'
values, energies, all forces and virial components, and finite-difference
coordinate and box derivatives. The triclinic position case guards the
full periodic-image vector used by its direct virial.

Optional GPU checks are in regtest/backend_checks.py. They compare driver
values and all coordinate/box derivatives and exercise invalid atom lists,
parameters, zero normalized weight and unsupported CUDA inputs. GPU device
ownership and context lifetime also need device-level tests on allocated GPUs.

Source formatting follows .astyle.options. Source files use the repository's
LGPL license. Hardware-specific performance data, models, trajectories,
absolute execution paths and private runtime records belong in the case
workspace. Reusable test structures in this plugin are synthetic small
systems.

The contribution guide requests new features against master, regression
tests, and preservation of existing tests. A bug fix should target the
oldest maintained branch containing the affected action. Maintainers also
review ownership and maintenance responsibility for new modules.

## Review boundaries

The FP64 CUDA prototype supports one MPI rank, complete exponential
assignments and axis-aligned orthorhombic PBC or NOPBC. AUTO uses CPU for the
other currently supported CPU modes. Scalar/direct reductions and OPES
execute on CPU.

Official GPU CI, the final public API, plugin registration and broader GPU
architectures remain upstream integration items. Public submission should
include the focused test results and document these limits.
