Bioconductor
------------

MassSpecWavelet and impute are Bioconductor packages.

Test environments
-----------------

-   local OS X install, R 3.6.1
-   devtools win-builder, R oldrel
-   devtools win-builder, R release (3.5.2)
-   devtools win-builder, R devel (unstable) (2019-12-02 r77499)
-   travis.ci Linux, x64, R release 3.6.1 (2017-01-27)
-   travis.ci Linux, x64, R devel (unstable) (2019-12-10 r77548)

R CMD check results
-------------------

There were no ERRORs or WARNINGs.

Local OS X build generated 1 NOTE

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
