# BenchmarkFlatnessbasedControl

<!--[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://stephans3.github.io/BenchmarkFlatnessbasedControl.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://stephans3.github.io/BenchmarkFlatnessbasedControl.jl/dev/)
[![Build Status](https://github.com/stephans3/BenchmarkFlatnessbasedControl.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/stephans3/BenchmarkFlatnessbasedControl.jl/actions/workflows/CI.yml?query=branch%3Amain)-->

This repository contains simulation files for our article about the benchmarking of flatness-based control of the heat equation. Flatness-based control design is a method to compute open-loop control signals. 
The code is used to compute an input signal for a one-dimensional heat equation with Neumann-type boundary actuation. The simulation results are compared for pure aluminum and steel 38Si7. 

## Article

Please find here the article as PDF:

[Benchmarking of Flatness-based Control of the Heat Equation [Open Access]](https://www.asim-gi.org/fileadmin/user_upload_asim/ASIM_Publikationen_OA/AM190/a4703.arep.47_OA.pdf)

## Citation

```bibtex
@article{scholz2024benchmarking,
  title={Benchmarking of Flatness-based Control of the Heat Equation},
  author={Scholz, Stephan and Berger, Lothar and Lebiedz, Dirk},
  journal={Tagungsband Langbeitr{\"a}ge ASIM SST 2024},
  number={47},
  pages={103},
  year={2024},
  publisher={ARGESIM Publisher, Vienna}
}
```

## How to use the code

The source code is stored in the [source folder](https://github.com/stephans3/BenchmarkFlatnessbasedControl.jl/tree/main/src).

The relation between the figures in the article and the source code files is listed in the table below.

| Figure  | Explanation | Source Code File |
|---------|--|----------|
| 2  | Sequences $\eta_{i}$ and ratio $\frac{\eta_{i+1}}{\eta_{i}}$ | *sequence_eta.jl*  |
| 3  | Transition and 1. derivative | *transition_trajectory_vary_w.jl*  |
| 4  | Norm of bump function | *norms_derivative_omega.jl*  |
| 5  | Sequence $\mu_{i}$ and ratio  $\frac{\mu_{i}}{\max_{j\in \{1,...,i\}} \mu_{j} }$ | *input_signal_progress.jl*  |
| 6  | Approximated input signals | *input_signal_progress.jl*  |
| 7  | Input signals (final simulation) | *input_signal_final.jl*  |
| 7  | Solution of heat equation for aluminum | *heat_eq_aluminum.jl*  |
| 7  | Solution of heat equation for steel | *heat_eq_steel.jl*  |


The data of these plots is also available for [pgf | TikZ | LaTeX export](https://github.com/stephans3/BenchmarkFlatnessbasedControl.jl/tree/main/src/pgf_tikz_export).

## Computation of the derivatives

Please take a look on the Julia library [BellBruno.jl](https://github.com/stephans3/BellBruno.jl/)