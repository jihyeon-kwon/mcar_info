set.seed(1234)

#####################################################
# unstructured with cov matrix & indep normal prior #
#####################################################


n <- 8 # number of variances
q <- (n - 1) * n / 2 # number of correlations
# s_lower <- unlist(lapply(1:n, function(i) 0.1^((n - i + 1):1))) # step size
# s_diag <- 0.1^(n:1) # some step size (could be a vector)
s_lower <- 0.1
s_diag <- 0.1


# get elements of Chol lower Tri from previous run
G0 <- solve(rWishart(1, n+1, diag(7, n))[,,1]) # pretend this was our previous cov mat
G0
L0 <- t(chol(G0))
prev_lower <- L0[lower.tri(L0)]
prev_diag <- diag(L0)
G0 - L0 %*% t(L0)

# get new samples of Chol lower Tri
lower_el <- rnorm(q, prev_lower, s_lower) # sample lower triangle of L* (G = LL')
diag_el <- rnorm(n, prev_diag, s_diag)
L <- diag(diag_el)
L[lower.tri(L)] <- lower_el
L # this is actually L*
G <- L %*% t(L)

eigen(G0)$value
eigen(G)$value
cov2cor(G) # legal covariance matrix (always)

diag(G0)
diag(G)

################
# unstructured #
################


n <- 8
q <- (n - 1) * n / 2
lower_el <- rbeta(q, 1, 1)
L <- diag(1, n)
L[lower.tri(L)] <- lower_el
L


# my first attempt : just sample the lower triangle of the covariance matrix
G1 <- L + t(L) - diag(1, n)
isSymmetric(G1)
eigen(G1)$value # note negative eigen values
chol(G1) # doesn't work (not PD)

# my second attempt: sample the L matrix from the G = LL' Cholesky
G2 <- L %*% t(L)
isSymmetric(G2)
eigen(G2)$value # all positive
t(chol(G2)) # works (PD)
t(chol(G2)) - L # recovered (obvisouly)
# but how can I implement 1) sampling a negative correlation 2) evaluating the ration in MH?



#########
# AR(1) #
#########

# AR(1) structure guarantees PD so it's fine
rho <- rbeta(1, 1, 1)
first_row <- rho^(0:(n-1))
G3 <- toeplitz(first_row)
eigen(G3)$value # positive
chol(G3) # works (PD)


