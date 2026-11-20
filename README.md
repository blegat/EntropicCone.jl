# Entropic Cone

| **Documentation** | **Build Status** |
|:-----------------:|:----------------:|
| [![][docs-stable-img]][docs-stable-url] | [![Build Status][build-img]][build-url] |
| [![][docs-dev-img]][docs-dev-url]       | [![Codecov branch][codecov-img]][codecov-url] |

This branch uses JuMP 1 and the unreleased StructDualDynProg revision pinned in `Project.toml`. With Julia 1.12 or later, run `julia --project -e 'using Pkg; Pkg.instantiate()'` from the repository root to set up the package, tests, and documentation.

Package for approximating the Entropic Cone and solving optimization problems on it.

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **most recently tagged version of the documentation.**
- [**LATEST**][docs-dev-url] &mdash; *in-development version of the documentation.*

## Examples

Example notebooks are available in the [`examples` folder](examples).
We link them below with the literature.

### Reproducing

The linked notebooks reproduce the results of the following papers:

* [BJ16] [B. Legat](https://perso.uclouvain.be/benoit.legat), [R. M. Jungers](https://perso.uclouvain.be/raphael.jungers/content/home).
[*Parallel optimization on the Entropic Cone*](http://sites.uclouvain.be/sitb2016/Proceedings_SITB2016_preliminary.pdf).
[37rd Symposium on Information Theory in the Benelux](http://sites.uclouvain.be/sitb2016), **2016**:
[Ingleton Score result](examples/Parallel%20optimization%20on%20the%20Entropic%20Cone.ipynb).

### Exploring

The linked notebooks explores the examples of the following papers using this
package:

* [ZY97] Z. Zhang, R. W. Yeung.
*A non-Shannon-type conditional inequality of information quantities.*
IEEE Transactions on Information Theory, **1997**:
[Zhang-Yeung inequality](examples/Zhang-Yeung_inequality.ipynb).

## How to cite

To cite this package, use
```
@InProceedings{legat2016parallel,
  author    = {Legat, B. and Jungers, R. M.},
  title     = {Parallel optimization on the Entropic Cone},
  booktitle = {Proceedings of the 37rd Symposium on Information Theory in the Benelux},
  year      = {2016},
  series    = {SITB '16},
  pages     = {206--211},
}
```

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-stable-url]: https://blegat.github.io/EntropicCone.jl/stable
[docs-dev-url]: https://blegat.github.io/EntropicCone.jl/dev

[build-img]: https://github.com/blegat/EntropicCone.jl/actions/workflows/ci.yml/badge.svg
[build-url]: https://github.com/blegat/EntropicCone.jl/actions/workflows/ci.yml
[codecov-img]: http://codecov.io/github/blegat/EntropicCone.jl/coverage.svg?branch=main
[codecov-url]: http://codecov.io/github/blegat/EntropicCone.jl?branch=main
