# Transitioning from Particle Dynamics to Continuum Dynamics in Yeast Metabolism

## Introduction
Yeast metabolism offers an intriguing opportunity to investigate mean field theory in the context of many-particle systems with non-local interactions. In this context, we wish to understand how solutions to the $$N$$-particle system of ordinary differential equations (NODE) relate to solutions of the proposed mean field limit partial differential equation (MFPDE). This repository provides means to simulate both contexts numerically.

## Usage
I recommend starting any MATLAB script with a boiler plate segment to
1. establish the parameter structure
2. adjust the MATLAB path to include subdirectories

An example of such a code snippet is available at the top of `testing.m`.

## Functionality
* Simulation of NODE 
  * Forward Euler method (`simulation_methods/forward_euler.m`)
  * Forward Euler-Maruyama method (`simulation_methods/forward_euler_noise.m`)
* Simulation of MFPDE
  * Pseudospectral in space, backwards Euler in time (`simulation_methods/pseudospectral.m`)
* Sampling (`simulation_methods/simulation_utilities/construct_em.m`)
* Graphics
  * Still frame (`graphics/frame.m`)
  * Animation (`graphics/animate.m`)
* Metrics
  * $$L^1$$, accounting for thickness and translation (`metrics/metric_lp_1.m`)
  * Wasserstein (*in development*) (`metrics/metric_wasserstein.m`)
  * Maximum Mean Discrepancy (MMD) (*in development*)

## Credits
This material was used to support the doctoral thesis research of Adam Petrucci under the supervision of [Keith Promislow](https://users.math.msu.edu/users/promislo/) at [Michigan State University](https://math.msu.edu/). The original "yeast problem" takes inspiration from research led by [Todd Young](http://www.ohiouniversityfaculty.com/youngt/).
