# Error Estimates for Golub–Kahan Bidiagonalization with Tikhonov Regularization

This repository contains MATLAB code supporting the paper: "Error Estimates for Golub–Kahan Bidiagonalization with Tikhonov Regularization for Ill–posed Operator Equations" by Abdulaziz Alqahtani, Ronny Ramlau, and Lothar Reichel, published in *Inverse Problems*, Volume 39, 2023. DOI: https://doi.org/10.1088/1361-6420/aca754

This project provides a practical and theoretical framework for solving linear ill-posed operator equations \( Ax = y \) where both A and y are contaminated with noise. The method uses a continuous variant of the Golub–Kahan bidiagonalization (GKB) algorithm to construct a low-rank approximation of the operator, followed by Tikhonov regularization in the reduced subspace. A nonlinear equation is solved to determine an optimal regularization parameter. The approach is supported by rigorous error analysis and illustrated using several benchmark test problems from inverse problems and integral equations.

This repository includes implementations for continuous and discrete variants of the method, as well as visualizations and parameter selection strategies. The following MATLAB files are included, with their original filenames preserved:
- `Baart_GKB.m`: Applies continuous GKB and Tikhonov regularization to the Baart problem.
- `Baart_GKB_discrete.m`: Discrete version of the Baart example using matrix formulation.
- `Foxgood_GKB.m` and `Foxgood_GKB_discrete.m`: Continuous and discrete versions of the Foxgood problem.
- `Gravity_GKB.m` and `Gravity_GKB_discrete.m`: Continuous and discrete implementations for the Gravity test problem.
- `Shaw_GKB.m` and `Shaw_GKB_discrete.m`: Implements the method for the Shaw example in both continuous and discrete forms.
- `continuous_GKB.m`: Core function implementing Golub–Kahan bidiagonalization for continuous operators.
- `discrete_GKB.m`: Core function for discrete GKB.
- `fofalpha.m`: Adds noise to the kernel/operator to simulate a perturbed problem instance.
- `fig_code.m`: Generates plots corresponding to figures in the published paper.
- `newton_project3.m`: Implements Newton's method to select the regularization parameter based on a nonlinear equation.
- `numerical_examples_project_3.m`: Master script to run all examples, generate reconstructions, and test convergence.

Requirements:
- MATLAB R2017a or newer
- Chebfun toolbox (https://www.chebfun.org/)
- Regularization Tools by Hansen (https://www.mathworks.com/matlabcentral/fileexchange/52-regularization-tools)
- Optionally, IR Tools for 2D inverse diffusion problem examples (not bundled in this repo)

To get started:
1. Clone the repository: `git clone https://github.com/your-username/gkb-tikhonov-estimates.git` and open it in MATLAB.
2. Run the script `numerical_examples_project_3.m` to execute all included examples automatically and verify correctness.
3. You can also run any individual problem script manually, such as `Baart_GKB.m`, `Foxgood_GKB_discrete.m`, or `Gravity_GKB.m`, to experiment with different configurations.
4. To reproduce figures used in the paper (e.g., comparing error bounds and solution accuracy), run `fig_code.m`.

The test problems included are:
- Baart: A mildly ill-posed problem with an exponential kernel, useful for evaluating stabilization in smooth cases.
- Foxgood: A smooth symmetric problem with known exact solution.
- Gravity: A model derived from inverse gravity problems, highlighting performance under stronger ill-posedness.
- Shaw: A challenging example with oscillatory and localized structure in the solution.
- 2D Inverse Diffusion: Based on IR Tools (not included in this archive), showcasing scalability and applicability to PDE-based inverse problems.

If you use this code or build on this work, please cite the following publication:

@article{Alqahtani2023GKB,
author = {Abdulaziz Alqahtani and Ronny Ramlau and Lothar Reichel},
title = {Error Estimates for Golub--Kahan Bidiagonalization with Tikhonov Regularization for Ill--posed Operator Equations},
journal = {Inverse Problems},
volume = {39},
number = {2},
pages = {025002},
year = {2023},
doi = {10.1088/1361-6420/aca754}
}

For questions, suggestions, or collaboration inquiries, feel free to reach out via email: abdulaziz.alqahtani@protonmail.com

This project is licensed under the Creative Commons Attribution 4.0 International License. You are free to share and adapt the code with appropriate credit. License details: https://creativecommons.org/licenses/by/4.0/
