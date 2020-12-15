# muscat

<details>

* Version: 1.4.0
* GitHub: https://github.com/HelenaLC/muscat
* Source code: https://github.com/cran/muscat
* Date/Publication: 2020-10-27
* Number of recursive dependencies: 208

Run `revdep_details(, "muscat")` for more info

</details>

## Newly broken

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      Warning (test-aggregateData.R:2:1): (code run outside of `test_that()`)
      Warning (test-aggregateData.R:2:1): (code run outside of `test_that()`)
      Warning (test-mmDS.R:45:5): mmDS() - filtering
      Warning (test-mmDS.R:45:5): mmDS() - filtering
      Warning (test-mmDS.R:68:9): mmDS-utils
      Warning (test-mmDS.R:68:9): mmDS-utils
      FAILURE (test-pbDS.R:34:9): defaults - pbDS.limma-trend
      FAILURE (test-pbDS.R:65:13): pbDS.limma-trend
      FAILURE (test-pbDS.R:65:13): pbDS.limma-trend
      Warning (test-pbHeatmap.R:2:1): (code run outside of `test_that()`)
      Warning (test-pbHeatmap.R:35:5): pbHeatmap() - input arguments
      
      [ FAIL 3 | WARN 15 | SKIP 0 | PASS 433 ]
      Error: Test failures
      Execution halted
    ```

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  9.1Mb
      sub-directories of 1Mb or more:
        data   2.3Mb
        doc    6.0Mb
    ```

