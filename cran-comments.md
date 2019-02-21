Bioconductor
------------

MassSpecWavelet and impute are Bioconductor packages.

Test environments
-----------------

-   local OS X install, R 3.4.3
-   devtools win-builder, R release (3.5.2)
-   devtools win-builder, R devel (unstable) (2019-02-19 r76133)
-   travis.ci Linux, x64, R oldrel (3.4.4)
-   travis.ci Linux, x64, R release (3.5.2)
-   travis.ci Linux, x64, R devel (unstable) (2019-02-21 r76143)

R CMD check results
-------------------

There were no ERRORs or WARNINGs.

Local OS X build generated no NOTEs

Travis.ci generated no NOTEs

Win-builder generated 1 NOTE

-   Uses the superseded package: 'doSNOW': this is to provide
    continuously updated progress bars when performing parallel
    computations (not available in doParallel so far as doParallel calls
    the combine function only when all results have been accumulated).
    When this option is avalaible in doParallel the switch to this
    package will be made.

-   possibly mis-spelled words in DESCRIPTION: These words are
    not mis-spelled.

Downstream dependencies
-----------------------

I have also run R CMD check on downstream dependencies of speaq. All
packages passed except specmine for which certain dependencies could not
be installed.
