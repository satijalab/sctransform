## Test environments
* x86_64-apple-darwin15.6.0 (64-bit), R 3.5.0
* x86_64-pc-linux-gnu (64-bit), R 3.5.3
* devtools::build_win(version = 'R-devel') - R 3.6.0 alpha (2019-04-09 r76362)
* devtools::build_win(version = 'R-release') - R 3.5.3 (2019-03-11)

## R CMD check results

There were no ERRORs or WARNINGs. 

There were 3 NOTEs:

    * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Christoph Hafemeister <chafemeister@nygenome.org>'
    Possibly mis-spelled words in DESCRIPTION:
      Hafemeister (11:29)
      Satija (11:45)

No, these names are not mis-spelled.

    ** running examples for arch 'i386' ... [49s] NOTE
    Examples with CPU or elapsed time     10s
                    user system elapsed
    smooth_via_pca 12.04   0.07   12.21
    correct        10.53   0.06   10.61
    vst            10.43   0.00   10.47
    correct_counts 10.34   0.02   10.42
    ** running examples for arch 'x64' ... [51s] NOTE
    Examples with CPU or elapsed time     10s
                    user system elapsed
    smooth_via_pca 12.36   0.00   12.48
    correct        11.00   0.02   11.10
    correct_counts 10.92   0.03   10.98
    vst            10.42   0.02   10.48

Yes, these examples take slightly longer than they should, but I think it 
is more important to provide meaningful example data here than to shave 
off a second or two.
