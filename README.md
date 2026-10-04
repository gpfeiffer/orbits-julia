# orbits-julia
A collection of Jupyter notebooks and Julia programs for orbit calculations

[![Open in Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gpfeiffer/orbits-julia/main)

## Contents

* `notebooks/`: the book *Computational Aspects of Complex Reflection Groups*,
  one notebook per chapter (`orbits`, `coxeter`, `enumerate`, `linear`),
  with a `preface`.  Other notebooks there, such as `core-topics`, are drafts
  and not part of the book.
* `gap3/`, `gap4/`: the GAP originals of many of the algorithms.
* `nb2book.py`, `templates/book/`, `Makefile`: the production line for the
  book PDF.

## Running the notebooks

The notebooks use Julia 1.13 and the packages in `Project.toml`, among them
[OrbitAl.jl](https://github.com/gpfeiffer/OrbitAl.jl).  To install them:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This also installs IJulia, whose build step installs the Jupyter kernel
`julia-1.13` that the notebooks ask for.  The kernel starts Julia with
`--project=@.`, so it finds this repository's `Project.toml`.  If the kernel
is missing, `julia --project=. -e 'using Pkg; Pkg.build("IJulia")'` installs
it.  Then start Jupyter in `notebooks/`.  Or run the notebooks online, through
Binder.

The notebooks are kept **without outputs**: please clear them before committing.

## Building the book

```sh
make            # executes the notebooks and builds exports/book.pdf
make clean      # removes _book_build/ and exports/book.pdf
```

A full build takes a few minutes, mostly for the conjugacy classes of $E_7$.
It needs, besides Julia as above:

* Python 3 with `nbconvert` and `nbformat`,
* `pandoc` (3.x),
* `xelatex`, with the LaTeX packages `tcolorbox` and `adjustbox`,
* the font DejaVu Sans Mono (or Liberation Mono, or Noto Sans Mono), for the
  Unicode symbols in the Julia code.

`nb2book.py` executes each notebook in memory, so the notebooks on disk are
never modified.  It moves the exercise solutions into an appendix, turns the
notebooks' coloured boxes into LaTeX boxes, exports everything to LaTeX
through `templates/book/`, and compiles with `xelatex`; its docstring has the
details.  Intermediate files go to `_book_build/`.

To add a chapter, list its notebook in `CHAPTERS` in `nb2book.py` **and** in
`NOTEBOOKS` in the `Makefile`.
