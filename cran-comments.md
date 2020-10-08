## Test environments
* Local install x86_64-apple-darwin17.0 (64-bit), R 4.0.2
* Winbuilder
  * devtools::check_win_devel() - x86_64-w64-mingw32 (64-bit), R Under development (unstable) (2020-10-05 r79298)
  * devtools::check_win_release() - x86_64-w64-mingw32 (64-bit), R 4.0.2
* rhub::check_for_cran(platforms = c('ubuntu-gcc-release', 'debian-gcc-devel', 'fedora-clang-devel'))
  * Ubuntu Linux 16.04 LTS, R-release, GCC - R 3.6.1
  * Debian Linux, R-devel, GCC - R Under development (unstable) (2020-10-02 r79291)
  * Fedora Linux, R-devel, clang, gfortran - R Under development (unstable) (2020-10-02 r79291)

## R CMD check results

0 ERRORs, 0 WARNINGs, 1 NOTE

```
* checking package dependencies ... NOTE
Packages which this enhances but not available for checking: 'speedglm', 'glmGamPoi'
```

`glmGamPoi` is an entirely optional package that is not required for core functionality, but only needed for alternative/faster implementations of the methods in this package. It is only available on Bioconductor. 

## Reverse dependencies

There are 2 reverse dependencies (muscat, Seurat). We checked them comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems


