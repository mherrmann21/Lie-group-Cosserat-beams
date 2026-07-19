[![MATLAB](https://img.shields.io/badge/language-MATLAB-blue.svg)](https://www.mathworks.com/products/matlab.html)
![MATLAB version](https://img.shields.io/badge/MATLAB-R2024b%2B-blue)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Paper DOI](https://img.shields.io/badge/paper_DOI-10.1016%2Fj.cma.2024.117367-blue)](https://doi.org/10.1016/j.cma.2024.117367)

# Lie-group Cosserat beams

MATLAB research code for the dynamic and static simulation of geometrically
exact Cosserat beams on $SE(3)$. The repository implements Lie-group
variational integrators (LGVIs) in absolute- and relative-kinematic
coordinates, including reduced beam models that constrain selected
deformation modes.

Moreover, the repository implements the simulation studies in [1] and Chapter 3 of [2].

> **Project status:** For new implementations, we recommend the actively
> developed MATLAB toolbox [ELARA](https://github.com/ELARA-Toolbox/ELARA).
> ELARA covers most of the functionality in this research repository and
> provides additional features.

## Features

- Geometrically exact, three-dimensional Cosserat beam models in absolute- and relative-kinematic spatial descriptions.
- Structure-preserving time integration with LGVIs.
- Full six-mode Simo-Reissner (SR) models (with bending, torsion, shear and extension) and reduced 
  models such as Kirchhoff (KH) beams.
- Implicit solution of the relative-kinematic models with Broyden's good method or a linear-time
  solver based on recursive articulated-body
  inertia algorithms [4,5].
- A continuous-time absolute-kinematic model for use with MATLAB ODE solvers.
- Static and dynamic simulations that include gravity, body-fixed and spatial external forces (wrenches), time-dependent spatial loads, Kelvin-Voigt-type strain damping, quadratic air drag, and additional
  node-attached rigid-body inertias.
- Post-processing, plotting, 3-D
  visualization, and animation functionality.
- MEX code generation with MATLAB Coder of important core functions to achieve high computational speed.

## Available models

| Formulation | Time model / solver | Deformation modes | Main entry point |
| --- | --- | --- | --- |
| Absolute kinematics | Coupled LGVI | Simo-Reissner, all six modes | `beamMdlAbsKinLGVI_general` |
| Absolute kinematics | Per-node LGVI | Simo-Reissner, all six modes | `beamMdlAbsKinLGVI_perNode` |
| Absolute kinematics | Continuous-time model with MATLAB ODE solvers | Simo-Reissner, all six modes | `simulateBeam_absKin_cont` |
| Relative kinematics | LGVI with Broyden's good method | Full or reduced | `beamMdlRelKinVarInt_Broyden` |
| Relative kinematics | Linear-time recursive LGVI | Full or reduced | `beamMdlRelKinVarInt_Recursive` |
| Relative kinematics | Static equilibrium with load continuation and `fsolve` | Full or reduced | `computeStaticEquilibriumRelKin` |

The generated MEX functions use the same names with the suffix `_mex`. The common
model definitions used by the studies are collected in
[`defineSimStudyBeamModels.m`](src/sim-studies/defineSimStudyBeamModels.m).

## Requirements

The project baseline is MATLAB R2024b or newer. Building and running every
included workflow requires the following:

- MATLAB Coder (to generate the MEX functions used by the simulations)
- A [supported C++ compiler](https://www.mathworks.com/support/requirements/supported-compilers.html) (to compile the generated MEX functions)
- Optimization Toolbox (for the static equilibrium solver, `fsolve`)
- Signal Processing Toolbox (for the exponential-map helper `sinc` used by the static solver)
- Curve Fitting Toolbox (for convergence-rate estimation in the study evaluation scripts)
- A toolbox providing `eul2rotm`, for example Robotics System Toolbox (for the rubber-rod initial configuration and out-of-plane-load scripts)


The dynamic Cayley-map models themselves do not call Optimization Toolbox or
Signal Processing Toolbox. However, the supplied build script also builds the
static solver, so those products are needed to run the complete build
unchanged.

In principle, code generation is not required for the simulations in this repository;
however, to run the scripts without generated MEX functions, the suffix `_mex` must be removed manually from various function calls across the code base,
as there is no automatic check for compiled MEX files.

## Installation and quick start

Clone the repository and start MATLAB in its root directory:

```console
git clone https://github.com/mherrmann21/Lie-group-Cosserat-beams.git
cd Lie-group-Cosserat-beams
```

Then run the following in MATLAB:

```matlab
% Add the source, scripts, studies, tests, and build directory to the path.
startup_cosserat_beams

% Select a compiler once, if MATLAB has not already configured one.
mex -setup C++

% Generate all MEX functions in build/.
buildMexFuns

% Run the relative-kinematic Kirchhoff LGVI example.
example_steel_string
```

Generated C/C++ sources and platform-specific MEX binaries are written to
`build/` and intentionally excluded from version control. Rebuild them after
changing the generated entry-point source, changing MATLAB releases, or
moving to another platform.

The startup script must be run once per MATLAB session unless the repository
paths are managed by another mechanism.

## Example scripts

Two end-to-end examples are provided:

- [`example_steel_string.m`](scripts/example_steel_string.m): Dynamic
  simulation of a steel string using the relative-kinematic LGVI and a
  reduced inextensible Kirchhoff model.
- [`example_rubber_rod.m`](scripts/example_rubber_rod.m): Dynamic simulation
  of a rubber rod using a relative-kinematic Kirchhoff LGVI, an
  absolute-kinematic Simo-Reissner LGVI, and an absolute-kinematic
  Simo-Reissner model integrated with `ode15s`.

Before running the examples, run `startup_cosserat_beams` and `buildMexFuns` (if not done previously).
Plotting, animation, and saving are controlled by flags at the top of each
script.

### Test and validation scripts

The files in [`tests/`](tests) are unit tests and interactive validation
scripts:

- `validation_all_beam_models` simulates the available dynamic formulations
  under common conditions for comparison.
- `validation_step_force_function` visualizes all implemented temporal load
  scaling modes for comparison.
- `runCayleyMapTests.m` and `runExponentialMapTests.m` run unit tests for
  selected Lie-group functions under `src/math`.

## Running the simulation studies

The study scripts are in [`simulation-studies/`](simulation-studies):

| Data generation | Evaluation | Main result file |
| --- | --- | --- |
| `sim_study_convergence` | `sim_study_convergence_evaluation` | `convStudyResults.mat` |
| `sim_study_computational_efficiency` | `sim_study_computational_efficiency_evaluation` | `simTimeStudyResults.mat` |
| `static_out_of_plane_load_test` | Results and plots are produced by the same script | `table_end_pos` and figure files |

Before starting a study:

1. Run `startup_cosserat_beams` to add all required folders to the MATLAB path.
2. Run `buildMexFuns` to build the required MEX functions (if the functions have not already been built).
3. Review the settings block at the top of the script.
4. For the convergence study, run the generation script twice with
   `CONV_STUDY_TYPE = "space"` and `CONV_STUDY_TYPE = "time"`.
5. In the corresponding evaluation script, set `SIM_SUBFOLDER` to the
   timestamped folder created by the generation script.

By default, study data are saved under `results/runs/`, and plots are
generated under `results/plots/`. Both locations are excluded
from version control. Full studies use fine spatial and temporal grids and can
require substantial runtime, memory, and disk space.

## Main data structures

The simulation interface is organized around MATLAB value classes:

- `beamSimulation` is the main simulation class that combines a model, beam parameters, simulation parameters,
  solver settings, and results. Its `simulateModel` method executes the
  selected model.
- `beamParams` stores all beam parameters: length, material and cross-section data, distributed inertia, stiffness, and damping.
- `beamSimPars` stores the simulation parameters: reference and initial configurations, the time grid,
  gravity, external loads, and optional attached-body parameters.
- `beamSimModel` defines the beam model (deformation modes) and integrator for the simulation. It stores a model name, function handle, solver configuration,
  and relative-model reduction matrices.

    - The integrator (absolute/relative kinematics, LGVI/ODE) is selected by the function handle stored in `funHandle`; see above for the available model functions.

    - The deformation modes are selected by the `Ba` and `Bc` matrices stored in `reducedParams`.
- `beamSolverConfig` stores nonlinear-solver tolerances, iteration limits, and
  Jacobian update settings.
- `beamSimRes` stores complete simulation results, including time-step
  metadata, aggregate metadata, energies, and a
  `beamSimData` object.
- `beamSimData` stores the numerical results: node rotations, positions, velocities, momenta,
  segment strains, and output times.

See the class definitions in [`src/simulation/`](src/simulation) and the
examples for complete construction and execution patterns.

## Numerical and implementation notes

### Beam deformation modes and solver time step

- The absolute-kinematic models retain all six deformation modes. Use a
  relative-kinematic model to constrain deformation modes;
  excluding stiff deformation modes such as shear or extension is important to achieve fast simulations for slender and stiff beams.

- The discrete LGVI models require an appropriate time step $h$ to work correctly.
This time step depends on the used beam model (included deformation modes, e.g., SR or KH models) and the beam parameters (in particular, slenderness and stiffness).

- For slender beams, SR models usually require much smaller time steps than constrained models such as KH models.

- In practice, a time step that is too large can prevent the implicit solver
  from converging. For simulations executed with the `beamSimulation` class,
  the console reports an `exitCode`: `0` means failure, `1` means that every
  step satisfied the target tolerance, and `2` means that at least one step
  reached the iteration limit but remained within `errorMarginLimit`.

- If a simulation fails for given beam parameters and beam model, decrease the time step $h$ until the simulation runs.

- The required time step can depend on additional simulation parameters, such as the initial configuration or external forces;
generally, high-frequency excitation (e.g., due to step inputs or initial deformations) can cause numerical issues and may require finer time steps.

### General notes

- The recursive relative-kinematic solver is intended primarily for
  conservative or weakly dissipative cases. Use the relative Broyden or an absolute model
  for strongly dissipative simulations.
- The current dynamic models use a cantilever boundary condition with the
  first node fixed and initialize the beam with zero velocity.
- The per-node absolute LGVI does not currently implement external forces.
- Fine time steps produce very large in-memory result arrays. Results can be
  resampled in time before expensive plotting or saving:

  ```matlab
  beamSim.simRes = beamSim.simRes.interpolateSimResTime(hOutput);
  ```

  This replaces the fine time grid in the assigned result object; keep a copy
  first if the original resolution is required.

## Repository structure

```text
build/                 Generated MEX files and C/C++ code (ignored)
results/               Generated study data and plots (ignored)
scripts/               MEX build script and runnable examples
simulation-studies/    Reproducibility and evaluation scripts
src/
  beam-definition/     Continuous beam parameter classes
  beam-parameters/     Ready-to-use parameter sets
  integration/         Dynamic and static integrators
  kinematics/          Relative kinematics and beam Jacobians
  math/                SO(3) and SE(3) operations and retractions
  plotting/            Result plotting and figure export
  sim-studies/          Shared study helpers
  simulation/           Simulation, solver, and result classes
  visualization/        3-D beam rendering and animation
tests/                  Interactive validation scripts
startup_cosserat_beams.m
```

## Citation

If you use the relative-kinematic formulation or this implementation in
academic work, please cite [1]:

```bibtex
@article{Herrmann2024RelativeKinematic,
  author  = {Herrmann, Maximilian and Kotyczka, Paul},
  title   = {Relative-kinematic formulation of geometrically exact beam
             dynamics based on {Lie} group variational integrators},
  journal = {Computer Methods in Applied Mechanics and Engineering},
  volume  = {432},
  pages   = {117367},
  year    = {2024},
  doi     = {10.1016/j.cma.2024.117367}
}
```

For work specifically using the absolute-kinematic LGVI or the recursive
solver, also cite the corresponding methodological references [3] or [4,5].

## License

This repository is available under the MIT License. See [`LICENSE`](LICENSE).

## References

[1] M. Herrmann and P. Kotyczka, “Relative-kinematic formulation of geometrically exact beam dynamics based on Lie group variational integrators,” *Computer Methods in Applied Mechanics and Engineering*, vol. 432, art. 117367, 2024. [doi:10.1016/j.cma.2024.117367](https://doi.org/10.1016/j.cma.2024.117367)

[2] M. Herrmann, *Geometric Modeling and Optimal Control of Rigid-Flexible
   Robot Manipulators*, Ph.D. thesis, Technical University of Munich, 2026
   (in preparation).

[3] F. Demoures, F. Gay-Balmaz, S. Leyendecker, S. Ober-Blöbaum, T. S. Ratiu,
   and Y. Weinand, “Discrete variational Lie group formulation of
   geometrically exact beam dynamics,” *Numerische Mathematik*, vol. 130,
   no. 1, pp. 73–123, 2015.
   [doi:10.1007/s00211-014-0659-4](https://doi.org/10.1007/s00211-014-0659-4)

[4] J. Lee, C. K. Liu, F. C. Park, and S. S. Srinivasa, “A linear-time
   variational integrator for multibody systems,” in *Algorithmic Foundations
   of Robotics XII*, Springer Proceedings in Advanced Robotics, vol. 13,
   pp. 352–367, 2020.
   [doi:10.1007/978-3-030-43089-4_23](https://doi.org/10.1007/978-3-030-43089-4_23)

[5] J. Kim, *Lie Group Formulation of Articulated Rigid Body Dynamics*, 2012.
   [Project page](https://www.cs.cmu.edu/~junggon/tools/gear.html)
