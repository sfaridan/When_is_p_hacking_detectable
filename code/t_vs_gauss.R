# Define the integrand function
integrand <- function(x) {
  dnorm(x) * (dt(x, df = 185) - dnorm(x))^2
}
result <- sqrt(integrate(integrand, lower = -Inf, upper = Inf)$value)
result 

integrand <- function(x) {
  dnorm(x) * (dt(x, df = 287) - dnorm(x))^2
}
result <- sqrt(integrate(integrand, lower = -Inf, upper = Inf)$value)
result 

integrand <- function(x) {
  dnorm(x) * (dt(x, df = 261) - dnorm(x))^2
}
result <- sqrt(integrate(integrand, lower = -Inf, upper = Inf)$value)
result 



# Define the integrand function
integrand <- function(x) {
  dnorm(x-0.01) * (dt(x, df = 120) - dnorm(x))^2
}
sqrt(integrate(integrand, lower = -Inf, upper = Inf)$value)