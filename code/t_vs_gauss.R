# Define the integrand function
integrand <- function(x) {
  dnorm(x) * (dt(x, df = 51) - dnorm(x))^2
}

# Perform the numerical integration
result <- integrate(integrand, lower = -Inf, upper = Inf)

# Take the square root
sqrt(result$value)