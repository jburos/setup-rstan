

# Prior to the tutorial make sure that the script below runs without error on your R installation.
# What you need is a working installation of Stan: http://mc-stan.org/ .
# For installation instructions, see here: 
# https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

# After installation you should be able to run this script which should output
# some summary statistics and some pretty plots, :)

# Generating some fake data
set.seed(123)
y <- rbinom(30, size = 1, prob = 0.2016)

# Fitting a simple binomial model using Stan
if (!require("rstan"))
  install.packages("rstan", repos = 'http://cran.rstudio.com")

if (!require("loo"))
  install.packages("loo", repos = "http://cran.rstudio.com")

model_string <- "
data {
  int n;
  int y[n];
}

parameters {
  real<lower=0, upper=1> theta;
}

model {
  y ~ bernoulli(theta);
}

generated quantities {
  int yrep[n];
  real log_lik[n];
  
  for (i in 1:n) {
    yrep[i] <- bernoulli_rng(theta);
    log_lik[i] <- bernoulli_log(y[i], theta);
  }
}
"

stan_samples <- stan(model_code = model_string, data = list(y = y, n = length(y)) )

summary(stan_samples)

## produce output MD file
#if (!require("knitr")) 
#  install.packages("knitr", repos="http://cran.rstudio.org")
#
#knit('./test_rstan_install.Rmd')

if (!require("rmarkdown"))
  install.packages("rmarkdown", repos = "http://cran.rstudio.org")

rmarkdown::render('../projects/setup-rstan/test_rstan_install.Rmd'
                 , output_format = 'all')


