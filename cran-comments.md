Bioconductor
------------

MassSpecWavelet is a Bioconductor package.

Test environments
-----------------

-   local OS X install, R 3.3.2
-   devtools win-builder, R-release (3.4.1)
-   devtools win-builder, R devel (unstable) (2017-09-12 r73242)
-   travis.ci Linux, x64, R release (3.4.1)
-   travis.ci Linux, x64, R devel (unstable)

R CMD check results
-------------------

There were no ERRORs or WARNINGs.

Local OS X build generated NOTEs

Travis.ci generated 1 NOTE

Win-builder generated NOTEs

-   Uses the superseded package: 'doSNOW': this is to provide progress
    bars when performing parallel computations (not available in
    doParallel so far). When this option is avalaible in doParallel the
    switch to this package will be made.

-   possibly mis-spelled words in DESCRIPTION: These words are
    not mis-spelled.

-   Files ‘README.md’ or ‘NEWS.md’ cannot be checked without ‘pandoc’
    being installed: These are common .md files

Downstream dependencies
-----------------------

I have also run R CMD check on downstream dependencies of speaq. All
packages passed.
