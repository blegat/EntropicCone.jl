# Entropic Cone

| **Documentation** | **Build Status** |
|:-----------------:|:----------------:|
| [![][docs-stable-img]][docs-stable-url] | [![Build Status][build-img]][build-url] |
| [![][docs-dev-img]][docs-dev-url]       | [![Codecov branch][codecov-img]][codecov-url] |

With Julia 1.12 or later, run `julia --project -e 'using Pkg; Pkg.instantiate()'` from the repository root to set up the package, tests, and documentation.

Package for approximating the Entropic Cone and solving optimization problems on it.

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **most recently tagged version of the documentation.**
- [**LATEST**][docs-dev-url] &mdash; *in-development version of the documentation.*

## Examples

Example notebooks are available in the [`examples` folder](examples).
We link them below with the literature.

### Paper example

The Ingleton notebook implements the hierarchy studied in the following paper.

* [BJ16] [B. Legat](https://perso.uclouvain.be/benoit.legat), [R. M. Jungers](https://perso.uclouvain.be/raphael.jungers/content/home).
[*Parallel optimization on the Entropic Cone*](http://sites.uclouvain.be/sitb2016/Proceedings_SITB2016_preliminary.pdf).
[37rd Symposium on Information Theory in the Benelux](http://sites.uclouvain.be/sitb2016), **2016**:
[Ingleton Score result](examples/Parallel%20optimization%20on%20the%20Entropic%20Cone.ipynb).

### Ingleton bounds with Gurobi

The experiment minimizes the four-variable Ingleton score with joint entropy
at most one. Each run starts with a fresh model, visits every scenario and
keeps all feasibility cuts. The table gives the number of passes through the
first feasible forward pass and the corresponding stage-solve calls in
StructDualDynProg 0.3.0 (Julia 1.13.4, Gurobi 13.0.3, one solver thread).

| Model | Lower bound | Forward passes (calls) | Backward passes (calls) |
| --- | ---: | ---: | ---: |
| Notebook hierarchy, 7 variables, 3 stages | −0.1585365854 | 11 (16,017) | 13 (32,140) |
| DFZ and Csirmaz copy chains, one stage per copy | −5/32 | 134 (576) | 25 (211) |
| Chang copy chain, stages of 4, 6 and 10 variables | −19/124 | 68 (158) | 34 (116) |

The 2016 paper reports −3/19. The targeted experiments use two catalogue copy
programs and a four-step, ten-variable copy program. [Chang (2026), Eqs. (11)
and (20)](https://arxiv.org/html/2610.05168v1).
The root can be averaged over the swaps A↔B and C↔D because the normalized
almost-entropic region is convex and the swaps preserve the objective.
Each copy stage uses Shannon inequalities and ordinary self-adhesivity.

From the repository root, with a Gurobi license, run

```sh
julia --project=examples/ingleton -e 'import Pkg; Pkg.instantiate()'
julia --project=examples/ingleton examples/ingleton/hierarchy.jl 7 3 30
julia --project=examples/ingleton examples/ingleton/copy_programs.jl catalogue multistage backward 200
julia --project=examples/ingleton examples/ingleton/copy_programs.jl chang split backward 100
```

Use `forward` in either copy experiment for the matching control. The
`grouped` option combines all copies in each chain into one recourse LP.
The split and grouped runs use primal simplex and a feasibility tolerance of
1e−9 to improve numerical separation by the feasibility cuts.
The CSV reports the root bound before and after each pass. A finite
`forward_upper` confirms feasibility for the finite copy program.
An independent LP solve supplied a rational certificate consisting of 245
Shannon inequalities and 72 equalities. To verify −19/124 exactly, run
`julia examples/ingleton/verify.jl`. This uses only the Julia standard library.

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
