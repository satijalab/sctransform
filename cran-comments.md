## Test environments
* Local install x86_64-apple-darwin17.0 (64-bit), R 4.0.2
* Winbuilder
  * devtools::check_win_devel() - x86_64-w64-mingw32 (64-bit), R Under development (unstable) (2020-12-13 r79623)
  * devtools::check_win_release() - x86_64-w64-mingw32 (64-bit), R 4.0.3
* rhub::check_for_cran(platforms = c('ubuntu-gcc-release', 'debian-gcc-devel', 'fedora-clang-devel'))
  * Ubuntu Linux 16.04 LTS, R-release, GCC - R 3.6.1
  * Debian Linux, R-devel, GCC - R Under development (unstable) (2020-12-12 r79622)
  * Fedora Linux, R-devel, clang, gfortran - R Under development (unstable) (2020-10-24 r79367)

## R CMD check results

0 ERRORs, 0 WARNINGs, 1 NOTE

```
* checking package dependencies ... NOTE
Packages which this enhances but not available for checking: 'glmGamPoi'
```

`glmGamPoi` is an entirely optional package that is not required for core functionality, but only needed for alternative/faster implementations of the methods in this package. It is only available on Bioconductor. 

## Reverse dependencies

Tested using `revdepcheck::revdep_check`:
We checked 2 reverse dependencies (1 from CRAN + 1 from BioConductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
