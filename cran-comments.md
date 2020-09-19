## Test environments
* Local install x86_64-apple-darwin17.0 (64-bit), R 4.0.2
* Winbuilder
  * devtools::check_win_devel() - x86_64-w64-mingw32 (64-bit), R Under development (unstable) (2020-09-09 r79174)
  * devtools::check_win_release() - x86_64-w64-mingw32 (64-bit), R 4.0.2
* rhub::check_for_cran(platforms = c('ubuntu-gcc-release', 'debian-gcc-devel', 'fedora-clang-devel'))
  * Ubuntu Linux 16.04 LTS, R-release, GCC - R 3.6.1
  * Debian Linux, R-devel, GCC - R Under development (unstable) (2020-09-12 r79193)
  * Fedora Linux, R-devel, clang, gfortran - R Under development (unstable) (2020-09-13 r79194)

## R CMD check results

0 ERRORs, 0 WARNINGs, 2 NOTEs

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Christoph Hafemeister <christoph.hafemeister@nyu.edu>'

New maintainer:
  Christoph Hafemeister <christoph.hafemeister@nyu.edu>
Old maintainer(s):
  Christoph Hafemeister <chafemeister@nygenome.org>
```

My New York Genome Center email address is no longer functioning. I have switched to my New York University one. This is also the email address listed in the publication associated with this package's description (doi:10.1186/s13059-019-1874-1).

```
* checking package dependencies ... NOTE
Packages which this enhances but not available for checking: 'speedglm', 'glmGamPoi'
```
`speedglm` is an optional package not required for core functionality. It has not been updated in more than three years, hence the choice to put it in Enhance rather than Suggest.
`glmGamPoi` is an entirely optional package that is not required for core functionality, but only needed for alternative/faster implementations of the methods in this package. It is only available on Bioconductor. 

## Reverse dependencies

There are currently two packages importing sctransform: muscat, Seurat. This update does not impact their functionality.
