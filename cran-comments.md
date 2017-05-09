speaq2 vs speaq:
----------------

Both speaq and speaq2 was developed in our research group but speaq2 is
so fundamentally different than speaq that we deemed it belonged in a
different package. Speaq however will not be discontinued and is still
in an active state (both on the user as the developper end).

Bioconductor
------------

MassSpecWavelet is a Bioconductor package.

Test environments
-----------------

-   local OS X install, R 3.3.2
-   devtools win-builder, R-release (3.3.3)
-   devtools win-builder, R devel (3.4.0 alpha 72482)
-   travis.ci Linux, x64, R oldrel (3.2.5)
-   travis.ci Linux, x64, R release (3.3.3)
-   travis.ci Linux, x64, R devel (unstable r72488)

R CMD check results
-------------------

There were no ERRORs or WARNINGs.

Local OS X build generated no NOTEs

Travis.ci generated no NOTEs

Win-builder generated NOTE

-   Uses the superseded package: 'doSNOW': this is to provide progress
    bars when performing parallel computations (not available in
    doParallel so far). When this option is avalaible in doParallel the
    switch to this package will be made.

CRAN maintainer comments
------------------------

About the preference of not having version forks of CRAN packages. I try
to follow the package guidelines of Hadley Wickham, to tag every
accepted CRAN version with a release on GitHub and keeping the developer
version also there.

Downstream dependencies
-----------------------

There are currently no downstream dependencies for this package
