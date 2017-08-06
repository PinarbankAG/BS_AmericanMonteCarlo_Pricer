# Option + Market parameters
T <- 1.0        # Maturity 
r <- -0.35/100  # Risk-free p.a
q <- 3.40/100   # Div yield p.a.
S0 <- 1.0       # Initial stock lvl
Vol <- 0.151    # Volatility
K <- 1.00       # Strike
#Payoff function, given as a function of S (unique Udl)
Payoff <- function(S) {
  return (max(S - K, 0))
}

# Longstaff-Schwartz Parameters
RegDeg <- 5     # Regression Degree
nIter <- 10000  # MC drawings
nTime <- 252    # Time discretization

# Memory allocations
deltaT <- T / nTime
W <- matrix(0., nTime, nIter)
Option <- matrix(0., nTime, 2 * nIter)

# Brownian Drawings
W[,] <- rnorm(nIter * nTime, 0, sqrt(deltaT))
W_anti <- W * (-1)
Stock <- cbind(W, W_anti)

# Forward simulation of the underlying
Stock[1,] <- S0 * exp((r - q - 0.5 * Vol * Vol) * deltaT + Vol * Stock[1, ])
for(i in 2:nTime) {
  Stock[i, ] <- Stock[i-1, ] * exp((r - q - 0.5 * Vol * Vol) * deltaT + Vol * Stock[i, ])
}

# Backward LSM
Option[nTime, ] <- sapply(Stock[nTime, ], Payoff)
for(i in (nTime-1):1) {
  #Payoff computation at i
  ExercisePayoff <- sapply(Stock[i, ], Payoff)
  #Regression on ITM paths at i of E[Y = CV(i+1) | X = S(i)]
  Mask <- which(ExercisePayoff > 0.)
  Y <- Option[i + 1, Mask] * exp(-r * deltaT) 
  X <- Stock[i, Mask]
  PolynomialRegression <- lm(Y ~ poly(X, RegDeg))
  #Should we exercise ? Discount future cashflows - if any, and replace by today's CF if we should exercise.
  Option[i, ] = Option[i + 1, ] * exp(-r * deltaT)
  IsExercisedMask <- Mask[which(ExercisePayoff[Mask] > predict(PolynomialRegression, data.frame(X)))]
  Option[i, IsExercisedMask] = ExercisePayoff[IsExercisedMask]
}

#Option Value
OptionValue = sum(Option[1, ]) * exp(-r * deltaT) / (2 * nIter)
