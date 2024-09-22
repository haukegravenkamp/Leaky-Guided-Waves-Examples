# Leaky Guided Waves Examples

## purpose

This repository includes Matlab codes to compute dispersion curves and mode shapes of elastic waves in plate structures that are coupled to unbounded fluid or solid media on either or both surfaces.
Specifically, the codes reproduce the numerical examples published in

> [1] Gravenkamp, H., Plestenjak, B., Kiefer, D.A., and Jarlebring, E.. “Computation of Leaky Waves in Layered Structures Coupled to Unbounded Media by Exploiting Multiparameter Eigenvalue Problems.” Journal of Sound And Vibration, 2024. <https://doi.org/10.1016/j.jsv.2024.118716>.

**Note:** This repository includes only the previously computed finite element matrices for precisely the examples discussed in the paper. That is to say, this minimalistic repository is intended for reproducing the results in the paper and helping interested readers understand the formulation and code structure. If you want to create your own examples and apply this method to problems of practical interest, you can use the larger framework **SAMWISE**

> [2] Gravenkamp, H. “SAMWISE - Semi-Analytical Modeling of Waves in Structural Elements,” 2024. <https://github.com/haukegravenkamp/SAMWISE>.

SAMWISE lets you generate the (semi-analytical) finite element matrices for arbitrary layered structures as well as generic two-dimensional cross-sections. It incorporates the solvers for leaky waves from this repository.

## content

- four main files 'example*.m', one for each example published in the paper
- 'pathCode.m' adds the subfolders to the user's Matlab path for your convenience
- .\matrices contains the finite element matrices required for the four examples. They have been computed using SAMWISE.
- .\referenceSolutions contains previously computed solutions of the same example problems using alternative methods
- .\solvers contains the newly developed methods for solving the leaky waveguide problem
- .\MultiParEig is a toolbox for the solution of multiparameter eigenvalue problems by Bor Plestenjak. It is included here for convenience but can be obtained separately.

## usage

- download or clone the repository
- add the four subfolders to the Matlab path; you can do this manually or just run the script 'pathCode.m' in the parent directory.
- run any of the 'example*.m' files.

You should obtain the dispersion curves as published in the paper - if not, please let us know!
The simplest example 'example_brassWater.m' should take few seconds to run; the most expensive 'example_titaniumTeflonBrass.m' will take a few minutes.

## other related work

- The approaches presented here rely heavily on the toolbox MultiParEig [3].
- Some of the reference solutions are computed using the Global Matrix Method implemented in the software *disperse* [4].
- The reference solution of the brass plate immersed in water was obtained using the approach presented in [5].
- The finite element matrices were computed using high-order spectral elements as discussed in [6].
- Many papers discuss similar semi-analytical approaches to computing dispersion curves of plates and other prismatic structures in different contexts, see, in particular, the publications on the Thin Layer Method (TLM) or the Semi-Analytical Finite Element Method (SAFE), e.g. [7-10].

> [3] Plestenjak, B., “MultiParEig,” 2023. <https://www.mathworks.com/matlabcentral/fileexchange/47844-multipareig>.
>
> [4] Pavlakovic, B., Lowe, M.J.S., and Alleyne, D.N., “Disperse: A General Purpose Program for Creating Dispersion Curves.” In Review of Progress in Quantitative NDE, 185–92. Plenum Press, 1997. <https://doi.org/10.1007/978-1-4615-5947-4_24>.
>
> [5] Kiefer, D.A., Ponschab, M., Rupitsch, S.J., and Mayle, M., “Calculating the Full Leaky Lamb Wave Spectrum with Exact Fluid Interaction.” The Journal of the Acoustical Society of America 145, no. 6 (2019): 3341–50. <https://doi.org/10.1121/1.5109399>.
>
> [6] Gravenkamp, H., Song, C., and Prager, J., “A Numerical Approach for the Computation of Dispersion Relations for Plate Structures Using the Scaled Boundary Finite Element Method.” Journal of Sound and Vibration 331 (2012): 2543–57. <https://doi.org/10.1016/j.jsv.2012.01.029>.
>
> [7] Waas, G. “Linear Two-Dimensional Analysis of Soil Dynamics Problems in Semi-Infinite Layered Media.” PhD Thesis, University of California, Berkeley, 1972.
>
> [8] Aalami, B. “Waves in Prismatic Guides of Arbitrary Cross Section.” Journal of Applied Mechanics 40, no. 4 (1973): 1067–72.
>
> [9] Gavrić, L. “Finite Element Computation of Dispersion Properties of Thin-Walled Waveguides.” Journal of Sound and Vibration 173, no. 1 (1994): 113–24. <https://doi.org/10.1006/jsvi.1994.1221>.
>
> [10] Kausel, E. “Accurate Stresses in the Thin-Layer Method.” International Journal for Numerical Methods in Engineering 61 (2004): 360–79. <https://doi.org/10.1002/nme.1067>.

## authors

In addition to the paper [1] describing the theory behind this code, the implementation itself may be cited as

> Gravenkamp, H., Plestenjak, B., Kiefer, D. A., & Jarlebring, E. (2024), Leaky Guided Waves Examples (v1.0.0), doi:10.5281/zenodo.13825263, <https://github.com/haukegravenkamp/Leaky-Guided-Waves-Examples>

[![DOI](https://zenodo.org/badge/859846611.svg)](https://zenodo.org/doi/10.5281/zenodo.13825262)

