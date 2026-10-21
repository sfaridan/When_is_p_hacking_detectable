lam <- 2                                        # publication propensity exp(lam*t)
e   <- seq(0, .05, .01)                         # p-value bin edges
q   <- qnorm(1 - e/2)                           # matching |t| cutoffs; q[1] = Inf
P   <- function(a, b) pnorm(b, lam) - pnorm(a, lam) + pnorm(-a, lam) - pnorm(-b, lam)
h   <- P(q[-1], q[-6]); h <- h/sum(h)           # p-curve given p < .05

barplot(100*h, names.arg = c(".00-.01", ".01-.02", ".02-.03", ".03-.04", ".04-.05"),
        ylab = "% of significant results", xlab = "p-value", ylim = c(0, 60),
        main = "Evidential value (H=0; lambda=2)")

t <- seq(-8, 12, .001)                          # the proposition, numerically
f <- exp(lam*t)*dnorm(t); f <- f/(sum(f)*.001)  # tilt the null t-curve
max(abs(f - dnorm(t, lam)))                     # ~1e-16: exactly N(2,1)