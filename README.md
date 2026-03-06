# Discrete Sampling

Routines for sampling from discrete probability distributions[^1]

Currently supporting:
  * [CDF Inversion](https://www.pbr-book.org/3ed-2018/Monte_Carlo_Integration/Sampling_Random_Variables#x1-Example:Piecewise-Constant1DFunctions)
  * [Integer CDF Inversion](https://peterwonka.net/Publications/pdfs/2008.CGF.Cline.AComparisonOfTabularPDFInversionMethods.pdf)
  * [Alias Method](https://www.keithschwarz.com/darts-dice-coins/)
  * [Hierarchical Warping](https://cs.dartmouth.edu/~wjarosz/publications/clarberg05wavelet.html)

Other references:
  * [A Comparison of Tabular PDF Inversion Methods](https://peterwonka.net/Publications/pdfs/2008.CGF.Cline.AComparisonOfTabularPDFInversionMethods.pdf), *Cline et al.*
  * [Visualizing Warping Strategies for Sampling Environment Map Lights](https://pharr.org/matt/blog/2019/06/05/visualizing-env-light-warpings), *Matt Pharr's blog*

### License

Your choice of MIT or BSD0.

[^1]: Technically, step functions rather than arbitrary discrete distributions
