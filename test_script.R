
Sys.getenv('TMPDIR')
.libLoc(new = tempdir())
install.packages('rstan', repos = 'http://cran.studio.com')
library(rstan)

