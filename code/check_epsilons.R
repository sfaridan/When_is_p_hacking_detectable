
#old calculation
delta = (2*0.7)/1400
integrand <- function(x) {
  dnorm(x) * (dnorm(x-1-delta /2) -dnorm(x-1))^2
}
eps1_old_nearly_und<-  sqrt(integrate(integrand, lower = -Inf, upper = Inf)$value) 

delta = (2*6.5)/3000
integrand <- function(x) {
  dnorm(x) * (dnorm(x-1-delta /2) -dnorm(x-1))^2
}
eps1_old_power_comparison<-  sqrt(integrate(integrand, lower = -Inf, upper = Inf)$value) 

#Corrected eps1
## Standard normal density
phi <- function(t) dnorm(t)

## For given x,y compute:
## sqrt( ∫ φ(t) (φ(t - x) - φ(t - y))^2 dt )
kernel_diff <- function(x, y, rel.tol = 1e-8) {
  integrand <- function(t) {
    phi(t) * (phi(t - x) - phi(t - y))^2
  }
  val <- integrate(integrand, lower = -Inf, upper = Inf, rel.tol = rel.tol)$value
  sqrt(val)
}

## Approximate max_{|x - y| <= delta/2} of that quantity
## over a search box x,y ∈ [-L, L] with a grid of size nx × ny
max_diff_over_delta <- function(delta, L = 6, nx = 3000,ny=100,
                                rel.tol = 1e-8) {
  xs <- seq(-L, L, length.out = nx)
  max_val <- 0
  
  for (x in xs) {
    ys <- seq(x+delta/2, x-delta/2,length.out = ny)
    for (y in ys) {
     
        val <- kernel_diff(x, y, rel.tol = rel.tol)
        if (val > max_val) {
          max_val <- val
        }
      
    }
  }
  
  max_val
}

## Example: delta = 0.01

eps1_new_nearly_und <- max_diff_over_delta(delta = 2*6/1400, L = 6, nx = 1400) #L was hardcoded to 6
c(eps1_old_nearly_und,eps1_new_nearly_und)

eps1_new_power_comparison <- max_diff_over_delta(delta = 2*6.5/3500, L = 6.5, nx = 3000)
c(eps1_old_power_comparison,eps1_new_power_comparison)


#Try setting L=1 and grid to 1500 for the nearly undetectable sim
max_diff_over_delta(delta = 2*1/1400, L = 1, nx = 1400)



integrand2 <- function(x) {
  dnorm(x) * (dnorm(x-6.5))^2
}
eps2<-  sqrt(integrate(integrand2, lower = -10, upper = 10)$value)  #sqrt( dnorm(L/sigma_Y)/sigma_Y/sqrt(2*pi) )/2 #approx error from grid support
