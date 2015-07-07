## ABOUT ##

The GOPT toolbox contains functions for optimization over the manifold of PSD (HPD) matrices, as well as a collection of functions for computing maximum-likelihood estimates for Elliptical Gamma Distributions (EGDs) and calculating the KL-divergence between EGDs.

Please see the end of this README for a list of the algorithms available.

This toolbox is copyright (C) 2014 by Reshad Hosseini & Suvrit Sra and is distributed under the terms of the GNU General Public License (GPL) version 3 (or later).

Contact: [Reshad Hosseini](mailto:reshad.hosseini@ut.ac.ir) or  [Suvrit Sra](mailto:suvrit@mit.edu)


## Quick installation guide ##

* Unzip and copy the whole gopt directory you just downloaded in a location of your choice on disk, say, in /my/directory/.

* Go to /my/directory/eg/ at the Matlab command prompt and execute 'install' (or 'install_eg'). You may save this path for your next Matlab sessions: follow the menu File . Set Path... and save.


##Directory structure##

<pre>
./ The top directory with README.md
 |manoptAuxiliary   - Added solver and manifold (use with manopt)
   |--- lbfgsSolver/       - Manifold BFGS
   |--- pdManifold/        - psd and hpd classes 
   |--- wolfeLinesearch/   - Line search algorithm satisfy Wolfe conditions
 | means-medians/   - algos for means and medians of PSD matrices
 | examples/        - some examples showing how to use GOPT
 | eg/              - algorithms specific to EGDs
 | kotz/            - algorithms specific to Kotz type distributions
 | ThirdParty/      - Third party tools
   |--- randraw/          - m-function written by Alex Bar-Guy
   |--- manopt/           - the manopt toolbox
   |--- prettyPlot/       - our version of Mark Schmidt's plot tools

</pre>
## HOW TO CITE ##

If you find manifold LBFGS or the PSD (HPD) manifold useful in your work, please cite the following two papers:

<pre>

@Article{gopt,
  title={Conic Geometric Optimization on the Manifold of Positive Definite Matrices},
  author={Sra, Suvrit and Hosseini, Reshad},
  journal={SIAM Journal on Optimization},
  volume={25},
  number={1},
  pages={713--739},
  year={2015},
  publisher={SIAM}
}

@Article{manopt,
  author  = {Nicolas Boumal and Bamdev Mishra and P.-A. Absil and Rodolphe Sepulchre},
  title   = {{M}anopt, a {M}atlab Toolbox for Optimization on Manifolds},
  journal = {Journal of Machine Learning Research},
  year    = {2014},
  volume  = {15},
  pages   = {1455--1459},
  url     = {http://www.manopt.org}
}

</pre>

If you are using functions for estimating KL-divergence between EGDs or ML-estimation of parameters of EGD, please cite the following paper:

<pre>
@Article{egd,
  author  = {Reshad Hosseini and Suvrit Sra and Lucas Theis and Matthias Bethge},
  title   = {Statistical Inference with the Elliptical Gamma Distribution},
  journal = {arXiv:1410.4812},
  year    = {2014},
}

</pre>


## ALGORITHMS ##

- Fixed point iteration for S-Divergence means
- Fixed point iteration for S-Divergence medians
- Sampling based algorithm (of M. Bacak) for the Riemannian / Karcher Mean
- Limited memory Riemannian BFGS algorithm for optimization over smooth manifolds using [manopt](http://manopt.org)
- Fixed point iterations for ML estimation for EGDs
- Numerical computation of KL-Divergence between EGDs
- Fixed point iterations for ML estimation for Kotz type distributions


## COMING UP ##

- Manifold optimization algorithms for Riemannian Dictionary Learning (based on the paper: A. Cherian, S. Sra (Dictionary learning, sparse coding with PSD matrices)
- Faster computation of matrix geometric means (S. Sra (2013))
- Implementation of methods for computing matrix means and medians


-------------------------------------------------------------------------------

Last updated *Thu Jul 02 14:34:20 EST 2015*

