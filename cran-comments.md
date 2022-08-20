## Test environments

* Local install x86_64 (Linux 5.4.0-81-generic 64-bit), R 4.1.2
* Winbuilder
  * devtools::check_win_devel() - x86_64-w64-mingw32 (64-bit), R Under development (unstable) (2020-12-13 r79623)
  * devtools::check_win_release() - x86_64-w64-mingw32 (64-bit), R 4.0.3
* rhub::check_for_cran(platforms = c('ubuntu-gcc-release', 'debian-gcc-devel', 'fedora-clang-devel'))
  * Ubuntu Linux 16.04 LTS, R-release, GCC - R 3.6.1
  * Debian Linux, R-devel, GCC - R Under development (unstable) (2020-12-12 r79622)
  * Fedora Linux, R-devel, clang, gfortran - R Under development (unstable) (2020-10-24 r79367)

## R CMD check results

0 errors | 0 warnings | 0 notes

