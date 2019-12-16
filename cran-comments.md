## Test environments
* x86_64-apple-darwin15.6.0 (64-bit), R 3.5.0
* x86_64-pc-linux-gnu (64-bit), R 3.6.1
* devtools::build_win(version = 'R-devel') - (unstable) (2019-12-09 r77545)

## R CMD check results

There were no ERRORs or WARNINGs. 

There was one NOTE:

    ** running examples for arch 'x64' ... [57s] NOTE
    Examples with CPU (user + system) or elapsed time > 10s
                    user system elapsed
    smooth_via_pca 11.34   0.02   11.37
    vst            10.53   0.01   10.54
    correct_counts 10.24   0.01   10.25
    correct        10.11   0.07   10.18

Yes, these examples take slightly longer than they should, but I think it 
is more important to provide meaningful example data here than to shave 
off a second or two.
