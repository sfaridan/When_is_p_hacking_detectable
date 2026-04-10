## Generates a figure with a suspicious-looking but unhacked t-curve
set.seed(123)

root <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"

setwd(paste0(root,"/figures"))

n <- 50000

# Mixture membership
g <- sample(1:2, size = n, replace = TRUE, prob = c(0.45, 0.55))

# Simulate signed t-like statistics from ordinary normals
t_signed <- numeric(n)
t_signed[g == 1] <- rnorm(sum(g == 1), mean = 0.2, sd = 1)
t_signed[g == 2] <- rnorm(sum(g == 2), mean = 3, sd = 1)

# Absolute value
t_abs <- abs(t_signed)

# Exact mixture density for |t|
x <- seq(0, 8, length.out = 2000)
dabsnorm <- function(x, mean = 0, sd = 1) {
  dnorm(x, mean = mean, sd = sd) + dnorm(-x, mean = mean, sd = sd)
}
dmix <- 0.65 * dabsnorm(x, 0.0, 0.70) +
  0.35 * dabsnorm(x, 2.2, 0.55)

op <- par(no.readonly = TRUE)
par(mar = c(4.2, 4.2, 3.2, 0.8), las = 1)

pdf("suspicious_but_not_hacked_tcurve.pdf", width = 7, height = 5)
hist(t_abs,
     breaks = seq(0, 20, by = 0.10),
     freq = FALSE,
     col = "gray85",
     border = "white",
     main = expression("Simulated t-curve"),
     xlab = expression("|t|"),
     ylab = "Density",
     xlim = c(0, 5))

#lines(x, dmix, lwd = 3, col = "black")

abline(v = 1.96, lty = 2, col = "gray40")
text(1.96, par("usr")[4] * 0.96, "1.96", pos = 4, cex = 0.85, col = "gray40")

dev.off()

