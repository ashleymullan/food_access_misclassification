library(magrittr)

log_likelihood <- function(theta, Y, Offset, X)
{
  p = exp(theta[1] + theta[2] * X[, 2])
  ll = sum(Y * log(p) + (Offset - Y) * log(1 - p))
  return(ll)
}

neg_log_likelihood <- function(theta, Y, Offset, X)
{
  ll = log_likelihood(theta = theta, Y = Y, Offset = Offset, X = X)
  return(-ll)
}

#calculate gradient of l(alpha,beta) at previous iteration's alpha_n/beta_n
calc_gradient <- function(theta, Y, Offset, X)
{
  p = exp(theta[1] + theta[2] * X[, 2])
  
  grad_mat = matrix(data = 0, 
                    nrow = 2, 
                    ncol = 1)
  
  ## d/dbeta0 
  grad_mat[1] = sum((Y - Offset * p) / (1 - p))
  
  ## d/dbeta1 
  grad_mat[2] = sum(X[, 2] * (Y - Offset * p) / (1 - p))

  return(grad_mat)
}

#calculate Hessian of l(alpha,beta) at previous iteration's alpha_n/beta_n
calc_hessian <- function(theta, Y, Offset, X)
{
  p = exp(theta[1] + theta[2] * X[, 2])
  
  hessian = matrix(data = 0, 
                   nrow = 2, 
                   ncol = 2)
  
  ## d2/dbeta02
  hessian[1, 1] = sum((p / (1 - p) * ((Y - Offset * p) / (1 - p) - Offset)))
  
  ## d2/dbeta0dbeta1
  hessian[1, 2] = sum(X[, 2] * (p / (1 - p) * ((Y - Offset * p) / (1 - p) - Offset)))
  hessian[2, 1] = sum(X[, 2] * (p / (1 - p) * ((Y - Offset * p) / (1 - p) - Offset)))

  ## d2/dbeta12
  hessian[2, 2] = sum(X[, 2] ^ 2 * (p / (1 - p) * ((Y - Offset * p) / (1 - p) - Offset)))
  
  return(hessian)
}

# Load data from GitHub
food_access = read.csv(
  "https://raw.githubusercontent.com/sarahlotspeich/food/main/piedmont-triad-data/analysis_data.csv"
) |> 
  dplyr::select(-CountyName, -Y_BPHIGH, -Y_CHD, -Y_OBESITY) 

## Create binary indicators from X and Xstar
food_access = food_access |> 
  dplyr::mutate(
    Xstar_0.5 = as.numeric(Xstar <= 0.5), ### 0.5 miles
    Xstar_1 = as.numeric(Xstar <= 1), ### 1.0 miles
    Xstar_5 = as.numeric(Xstar <= 5), ### 5.0 miles
    Xstar_10 = as.numeric(Xstar <= 10) ### 10.0 miles
  ) |> 
  tidyr::gather(key = "THRESH", value = "binXstar", -c(1:6)) |> 
  dplyr::mutate(
    THRESH = as.numeric(sub(pattern = ".*_", 
                            replacement = "", 
                            x = THRESH)), 
    binX_full = as.numeric(X_full <= THRESH), 
    binX_partial = as.numeric(X_partial <= THRESH)
  ) |> 
  dplyr::filter(THRESH == 1)
nlm_fit <- nlm(f = neg_log_likelihood, 
               p = c(log(0.5), 0),
               Y = food_access$Y_DIABETES, 
               Offset = food_access$O_POP, 
               X = data.matrix(cbind(int = 1, food_access$binX_full)))

glm_fit <- glm(cbind(Y_DIABETES, O_POP - Y_DIABETES) ~ binX_full, 
               family = binomial(link = log),
               data = food_access)

X_mat <- data.frame(Int = 1, binX_full = food_access$binX_full) %>% data.matrix()
Y = data.matrix(food_access$Y_DIABETES) 
Offset = data.matrix(food_access$O_POP) 
X = X_mat
newton_raphson(Y = data.matrix(food_access$Y_DIABETES), 
               Offset = data.matrix(food_access$O_POP), 
               X = X_mat)

newton_raphson <- function(Y, Offset, X, TOL = 1e-4, MAX_ITER = 200)
{
  theta <- matrix(data = c(log(0.5), 0), 
                  nrow = 2, 
                  ncol = 1)
  Delta_l <- Inf
  l <- log_likelihood(theta = theta, Y = Y, Offset = Offset, X = X)
  iter <- 0
  while(abs(Delta_l) > TOL & iter < MAX_ITER)
  {
    iter <- iter + 1
    g <- calc_gradient(theta = theta, Y = Y, Offset = Offset, X = X)
    hess <- calc_hessian(theta = theta, Y = Y, Offset = Offset, X = X)
    H_inv <- solve(hess)
    
    Delta <- H_inv %*% g
    
    theta <- theta + Delta
    print(theta)
    
    l_new <- log_likelihood(theta = theta, Y = Y, Offset = Offset, X = X)
    Delta_l <- l - l_new
    l <- l_new
  }
  return(theta)
}


X_mat <- data.frame(Int = 1, binX_full = food_access$binX_full) %>% data.matrix()
theta_0 <- matrix(1/2, nrow = 2, ncol = 2)
newton_raphson(X = X_mat, y = food_access$Y_DIABETES)

# Compare to existing
fit <- glm(cbind(Y_DIABETES, O_POP - Y_DIABETES) ~ binX_full, 
           family = binomial(link = log),
           data = food_access)
fit$coefficients

# Simulate data 
N <- 10000
alpha <- 1
beta1 <- 2
beta2 <- 3
X1 <- rnorm(N)
X2 <- rnorm(N)
p <- 1/(1+exp(-(alpha + beta1*X1 + beta2*X2)))
Y <- rbinom(N, 1, p)

# Try GLM
glm(formula = Y ~ X1 + X2, family = binomial)

# Try my code
X_mat <- data.frame(Int = 1, X1, X2) %>% data.matrix()
theta_0 <- matrix(0, nrow = 3, ncol = 1)
newton_raphson(X = X_mat, y = Y)
