#-------------------------------------------------------------------#
# A General Approach to Modeling Latent Variable Interactions
#  and Nonlinear Effects
#
# Algorithmically Derived Conditional Distributions
#
# Distribution-Solver.R
#
# Copyright Brian T. Keller 2024 all rights reserved
#-------------------------------------------------------------------#


## Compute Symbolic Distribution for x | y if
#   y ~ normal(a + b * x, var_y);
#   x ~ normal(mu_x, var_x);
#
compute_dist <- function(y, a, b, mu_x, var_y, var_x) {
    # Obtain substitutions
    y <- substitute(y)
    a <- substitute(a)
    b <- substitute(b)
    mu_x <- substitute(mu_x)
    var_y <- substitute(var_y)
    var_x <- substitute(var_x)
    
    # Check if they exist and obtain values, otherwise use name
       y <- tryCatch(eval(y), error = \(.) y)
    a <- tryCatch(eval(a), error = \(.) a)
    b <- tryCatch(eval(b), error = \(.) b)
    mu_x <- tryCatch(eval(mu_x), error = \(.) mu_x)
    var_y <- tryCatch(eval(var_y), error = \(.) var_y)
    var_x <- tryCatch(eval(var_x), error = \(.) var_x)
    
    # Equation (A3)
    list(
        mean = paste(
            '(', var_x, '*', b, '*', '(', y, '-', a, ')', '+',
            var_y, '*', mu_x, ')', '/',
            '(', var_x, '*', b, '*', b, '+', var_y, ')'
        ),
        variance = paste(
            '(', var_y, '*', var_x, ')', '/', '(', 
            var_x, '*', b, '*', b, '+', var_y, ')'
        )
    )
}


## Factorization:
#
# f(Y_1 | eta) * f(Y_2 | eta) * f(Y_3 | eta) * f(eta)
#
## Models:
#
# f(Y_1 | eta) -> Y_1 ~ normal(nu_1 + l_1 * eta, s2_e1)
# f(Y_2 | eta) -> Y_2 ~ normal(nu_2 + l_2 * eta, s2_e2)
# f(Y_3 | eta) -> Y_3 ~ normal(nu_3 + l_3 * eta, s2_e3)
# f(eta)       -> eta ~ normal(alpha, psi2)
#
## Calculate distribution for f(eta | Y_1, Y_2, Y_3)
#
# Step 1: f(eta | Y_1)           \propto f(Y_1 | eta) * f(eta)
# Step 2: f(eta | Y_1, Y_2)      \propto f(Y_2 | eta) * f(eta | Y_1)
# Step 3: f(eta | Y_1, Y_2, Y_3) \propto f(Y_3 | eta) * f(eta | Y_1, Y_2)
#

# Step 1: f(eta | Y_1) propto f(Y_1 | eta) * f(eta)
dist <- compute_dist(
    y = Y_1,
    a = nu_1,
    b = l_1,
    mu_x = alpha,
    var_y  = s2_e1,
    var_x  = psi2
)

# Step 2:  f(Y_2 | eta) * dist
dist <- compute_dist(
    y = Y_2,
    a = nu_2,
    b = l_2,
    mu_x = dist$mean,
    var_y  = s2_e2,
    var_x  = dist$variance
)

# Step 3: f(Y_3 | eta) * dist
dist <- compute_dist(
    y = Y_3,
    a = nu_3,
    b = l_3,
    mu_x = dist$mean,
    var_y  = s2_e3,
    var_x  = dist$variance
)

# Print as language
lapply(dist, str2lang)

# $mean
# ((s2_e2 * (s2_e1 * psi2)/(psi2 * l_1 * l_1 + s2_e1))/((s2_e1 * 
#     psi2)/(psi2 * l_1 * l_1 + s2_e1) * l_2 * l_2 + s2_e2) * l_3 * 
#     (Y_3 - nu_3) + s2_e3 * ((s2_e1 * psi2)/(psi2 * l_1 * l_1 + 
#     s2_e1) * l_2 * (Y_2 - nu_2) + s2_e2 * (psi2 * l_1 * (Y_1 - 
#     nu_1) + s2_e1 * alpha)/(psi2 * l_1 * l_1 + s2_e1))/((s2_e1 * 
#     psi2)/(psi2 * l_1 * l_1 + s2_e1) * l_2 * l_2 + s2_e2))/((s2_e2 * 
#     (s2_e1 * psi2)/(psi2 * l_1 * l_1 + s2_e1))/((s2_e1 * psi2)/(psi2 * 
#     l_1 * l_1 + s2_e1) * l_2 * l_2 + s2_e2) * l_3 * l_3 + s2_e3)
# 
# $variance
# (s2_e3 * (s2_e2 * (s2_e1 * psi2)/(psi2 * l_1 * l_1 + s2_e1))/((s2_e1 * 
#     psi2)/(psi2 * l_1 * l_1 + s2_e1) * l_2 * l_2 + s2_e2))/((s2_e2 * 
#     (s2_e1 * psi2)/(psi2 * l_1 * l_1 + s2_e1))/((s2_e1 * psi2)/(psi2 * 
#     l_1 * l_1 + s2_e1) * l_2 * l_2 + s2_e2) * l_3 * l_3 + s2_e3)