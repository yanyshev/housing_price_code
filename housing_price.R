housing_data = readxl::read_xlsx('~/housing_data.xlsx')
data_to_fit = t(housing_data$real_index)

# Initial parameters
par_init = list(
  alpha1 = 0.8,
  alpha2 = 0.3,
  
  sigma_mu = 0.0035,
  sigma_p_gap = 1.2
)
par_names = names(par_init)

# State Space Representation
ssm <- function(p, yt) {
  # Cov matrix
  P0 <- diag(0.0015, nrow=4, ncol=4)
  
  # State Error matrix
  E <- diag(c(par_init$sigma_mu, par_init$sigma_p_gap),
            nrow=2, ncol=2)
  R <- matrix(c(1, 0,
                1, 0,
                0, 1,
                0, 0),
              nrow=4, ncol=2, byrow=TRUE)
  Qm <- R %*% (E^2) %*% t(R)
  
  # Obs Error matrix
  Rm <- matrix(c(0), nrow=1, ncol=1, byrow=TRUE)
  
  # Initial guess of the unobserved components
  B0 <- matrix(c(0.4, 120, -20, -25), nrow=4, ncol=1, byrow=TRUE)
  
  # Constant in the state eq
  Dm <- matrix(rep(0, 4), nrow=4, ncol=1, byrow=TRUE)
  
  # Constant in the obs eq
  Am <- matrix(rep(0, 1), nrow=1, ncol=1, byrow=TRUE)
  
  # State eq transition matrox
  Fm <- matrix(c(1, 0, 0, 0,
                 1, 1, 0, 0,
                 0, 0, p$alpha1, p$alpha2,
                 0, 0, 1, 0), nrow=4, ncol=4, byrow=TRUE)
  
  # Obs eq matrix
  Hm <- matrix(c(0, 1, 1, 0), nrow=1, ncol=4, byrow=TRUE)
  
  
  return(list(B0 = B0, P0 = P0, Am = Am, Dm = Dm, Hm = Hm, Fm = Fm, Qm = Qm, Rm = Rm))
}

param_list_to_vector <- function(k) {
  # Converts named list to vector, applicable for MaxLik()
  c <- list(
    alpha1 = k[1], alpha2 = k[2], sigma_mu = k[3], sigma_p_gap = k[4])
  return(c)
}

library(kalmanfilter) # source: https://cran.r-project.org/web/packages/kalmanfilter/vignettes/kalmanfilter_vignette.html
library(maxLik) # source: https://cran.csiro.au/web/packages/maxLik/vignettes/using-maxlik.pdf
# Constraint matrices
ineqA <- matrix(c(1, 0, 0, 0,
                  -1, 0, 0, 0,
                  0, 1, 0, 0,
                  0, -1, 0, 0), nrow=4, ncol=4, byrow=TRUE)
ineqB <- matrix(c(0, 1, 0, 1), nrow=4, ncol=1, byrow=TRUE)

# Estimation
objective <- function(p) {
  ssm <- ssm(param_list_to_vector(p), data_to_fit)
  return(kalman_filter(ssm, data_to_fit, smooth = FALSE)$lnl)
}

init_to_fit <- as.numeric(par_init)
solve <- maxBFGS(objective, start = init_to_fit, constraints = list(ineqA = ineqA, ineqB = ineqB))

new_params <- as.list(lastFuncParam)
names(new_params) <- par_names

kf <- kalman_filter(ssm(new_params, data_to_fit), data_to_fit, smooth = FALSE)

# Save data
library(data.table)
unobserved <- data.table(t(kf[["B_tt"]]))[, "date" := housing_data$quarter]
colnames(unobserved) <- c("mu", "price_trend", "price_gap", "price_gap_lag", "date")

estimated <- data.table(t(kf[["y_tt"]]))[, "date" := housing_data$quarter]
colnames(estimated) <- c("price", "date")

# Plot the results
library(zoo)
estimated$date <- as.yearqtr(estimated$date, format = "%YQ%q")
unobserved$date <- as.yearqtr(unobserved$date, format = "%YQ%q")

library(ggplot2)
ggplot(estimated, aes(x=date, y=price)) + geom_line(group = 1) + geom_line(data=unobserved, aes(x=date, y=price_trend), linewidth = 1, group = 1)

write.csv(unobserved, '/Users/yanyshev_dima/Desktop/unobserved.csv')
write.csv(estimated, '/Users/yanyshev_dima/Desktop/estimated.csv')