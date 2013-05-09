hivclust
========

HIV clustering tools

to clone this R package from the git repository:
* if there are issues like "Permission denied (publickey)" , try git clone http://github.com/olli0601/hivclust.git

to build this R package:
* R CMD build pkg: builds the package (generates an archive called foo.tar.gz).
* R CMD build pkg --compact-vignettes --resave-data: same, making the package as small as possible.
* R CMD check hivclust_1.0-0.gz --as-cran: runs package quality checks similar to the checks made by CRAN; must be passed without error/warning before.

to install this R package:
* R CMD INSTALL hivclust_1.0-0.tar.gz

to re-build the R package:
* R CMD build pkg

to re-install the R package:
* R CMD INSTALL  hivclust_1.0-0.tar.gz
