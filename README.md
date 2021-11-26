# Coad.jl

This library is a Julia-based re-implementation of the flux-based numerical
solution to the stochastic collection equation developed by Andreas Bott in
his [well-known 1998 JAS paper (B98)][B98]. Bott provides legacy FORTRAN
implementations on his [academic website][bott_fortran], and while these are a
useful reference, they have a few idiosyncrasies that make them challenging
for new researchers in cloud microphysicists: 

- They are implemented in legacy FORTRAN 77 conventions with `common` and `data`
  blocks rather than use proper subroutines with clearly declared argument
  intent.
- The variable names and conventions leverage legacy implicit conventions for
  typing with which contemporary coders may not be familiar.
- The chosen variable names and overall structure of the code can be difficult
  to relate to the actual equations in [B98][].

## Stochastic Collection Equation

The SCE describes the evolution of a droplets size distribution over time due to
random collision events sampled from elements of droplet population:

![sce_eq](https://latex.codecogs.com/svg.image?\begin{align*}\frac{\partial&space;n(x,t)}{\partial&space;t}&space;&=&space;\int\limits_{x_0}^{x_1}&space;n(x_c,&space;t)K(x_c,x')n(x',&space;t)\,dx'&space;\\&space;&-&space;\int\limits_{x_0}^\infty&space;n(x,t)K(x,&space;x')n(x',t)\,&space;dx'\end{align*}&space;)

## Julia Implementation

The implementation presented here is a 50/50 re-write of the original FORTRAN.
The core algorithm (inner loop) is more-or-less unchanged from the original
version. However, we liberally leverage the core features of the Julia language
to enable future development and extensibility while retaining very fast
computational performance. The end result is an implementation which is about
as fast as the original FORTRAN code but far more practically useful.

We expect that the code here could easily be implemented in cloud microphysical
or other process models boot-strapped in Julia.

## Example Application

We've included a very basic demo program reproducing Figure (3) of [B98][] in
the script [`examples/simple.jl`](). In this example we initialize a basic
exponential droplet size distribution and configure the model to use the Long
(1974) collision kernel. This results in a droplet autoconversion simulation
where a secondary raindrop mode "appears" out of nowhere after about 15 minutes
in the simulation.

The result of this simulation can be visualized (e.g. in the demo Python Jupyter
Notebook at [`examples/plot_brm_ridgeline.ipynb`]() to produce a "ridgeline"
style plot showing the evolution of the multi-modal droplet population over
timme:

![BR74 reproduction](/examples/br74_example.pngexamples/br74_example.png)

## Acknowledgments

The core algorithm and basis implementation that inspired this project are
originally by [Andreas Bott](mailto:a.bott@uni-bonn.de).

[B98]: https://doi.org/10.1175/1520-0469(1998)055<2284:AFMFTN>2.0.CO;2
[bott_fortran]: https://www2.meteo.uni-bonn.de/forschung/gruppen/tgwww/people/abott/fortran/fortran_english.html