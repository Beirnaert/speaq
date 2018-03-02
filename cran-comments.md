Bioconductor
------------

MassSpecWavelet is a Bioconductor package.

Test environments
-----------------

-   local OS X install, R 3.4.3
-   devtools win-builder, R-release (3.4.3)
-   devtools win-builder, R devel (unstable) (2018-03-01 r74337)
-   travis.ci Linux, x64, R oldrel (3.3.3)
-   travis.ci Linux, x64, R release (3.4.2)
-   travis.ci Linux, x64, R devel (unstable) (2018-03-01 r74329)

R CMD check results
-------------------

There were no ERRORs or WARNINGs.

Local OS X build generated no NOTEs

Travis.ci generated no NOTEs

Win-builder generated NOTEs

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
packages passed.
