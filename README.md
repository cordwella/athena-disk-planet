athena-disk-planet
======
Modified version of Athena++ for planet-disc interactions. Pairs nicely with [this disc\_planet\_analysis](https://github.com/cordwella/disc_planet_analysis) code. 

Please see the [Athena++ main GitHub](https://github.com/PrincetonUniversity/athena) for 
up to date changes on the underlying hydrodynamics code.

## New feature overview
This fork adds two new problem files, `diskplanet_3d_sph.cpp` and `diskplanet_2d.cpp`, and two additional utilities, one for the calculation of Bessel functions based off the Cephes library and the other for a custom radial output type (this is used in `diskplanet_2d.cpp` to output $F\_wave$ and $dT/dR$, and originally written for Cimerman & Rafikov 2021). In addition, compile scripts, slurm scripts, and input files are provided. These are designed for the University of Cambridge Faculty of Mathematics's [Swirles HPC Cluster](https://www.maths.cam.ac.uk/computing/faculty-hpc-system-swirles), but more generally I hope they can be used as a starting point for researchers begining to use Athena++ for planet-disc interaction problems. Those scripts can be found in the `cambridge_scripts` directory.

It also adds a 1D output for 2D simulations. This output code is based on that of Nicolas Cimerman, developed for Cimerman & Rafikov (2021), "Planet-driven density waves in protoplanetary discs: Numerical verification of non-linear evolution theory"

### Setting up 2D planet-disc simulations

The potential $\Phi_{B, H\_p}$ is reccomended. Set this using
```
potential_order = -2 
```
if instead you must the traditional second order smoothed potential use
```
potential_order = 2
eps = 0.65
```
where `eps` is set to your smoothing length (0.65 is reccomended for the best match to the one-sided Lindblad torque).

## Citing this repository

Please first cite the main Athena++ repository and paper (details below) and then cite that you have used problem files generated for Cordwell, Ziampras, Brown & Rafikov (2025). Cordwell et al (2025) is also the correct citation for disc\_planet\_analysis.

# The original Athena++ README.md is as follows
<p align="center">
	  <img width="345" height="345" src="https://user-images.githubusercontent.com/1410981/115276281-759d8580-a108-11eb-9fc9-833480b97f95.png">
</p>

Athena++ radiation GRMHD code and adaptive mesh refinement (AMR) framework

Please read [our contributing guidelines](./CONTRIBUTING.md) for details on how to participate.

## Citation
To cite Athena++ in your publication, please use the following BibTeX to refer to the code's [method paper](https://ui.adsabs.harvard.edu/abs/2020ApJS..249....4S/abstract):
```
@article{Stone2020,
	doi = {10.3847/1538-4365/ab929b},
	url = {https://doi.org/10.3847%2F1538-4365%2Fab929b},
	year = 2020,
	month = jun,
	publisher = {American Astronomical Society},
	volume = {249},
	number = {1},
	pages = {4},
	author = {James M. Stone and Kengo Tomida and Christopher J. White and Kyle G. Felker},
	title = {The Athena$\mathplus$$\mathplus$ Adaptive Mesh Refinement Framework: Design and Magnetohydrodynamic Solvers},
	journal = {The Astrophysical Journal Supplement Series},
}
```
Additionally, you can add a reference to `https://github.com/PrincetonUniversity/athena` in a footnote.

Finally, we have minted DOIs for each released version of Athena++ on Zenodo. This practice encourages computational reproducibility, since you can specify exactly which version of the code was used to produce the results in your publication. `10.5281/zenodo.4455879` is the DOI which cites _all_ versions of the code; it will always resolve to the latest release. Click on the Zenodo badge above to get access to BibTeX, etc. info related to these DOIs, e.g.:

```
@software{athena,
  author       = {Athena++ development team},
  title        = {PrincetonUniversity/athena: Athena++ v24.0},
  month        = jun,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {24.0},
  doi          = {10.5281/zenodo.11660592},
  url          = {https://doi.org/10.5281/zenodo.11660592}
}
```
