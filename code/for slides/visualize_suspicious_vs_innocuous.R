############################################################
## Intuition for convex projection in base R
## A non-technical picture:
##   - blue ellipse = set of "allowable" explanations
##   - red point    = observed object
##   - dark point   = projection (closest allowable point)
##   - dashed line  = distance minimized by projection
############################################################

root <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"

setwd(paste0(root,"/figures"))

# -----------------------------
# 1) Define a convex set: ellipse
# -----------------------------
theta <- seq(0, 2*pi, length.out = 500)

cx <- 0       # center x
cy <- 0       # center y
a  <- 2.8     # horizontal radius
b  <- 1.6     # vertical radius

x_ellipse <- cx + a * cos(theta)
y_ellipse <- cy + b * sin(theta)

# -----------------------------
# 2) Pick a point outside the set
# -----------------------------
x0 <- 3.8
y0 <- 2.4

x1 <- -1.5
y1 <- 1.6

# -----------------------------
# 3) Compute its projection onto the ellipse
#    (numerically: closest point on boundary)
# -----------------------------
dist2 <- (x_ellipse - x0)^2 + (y_ellipse - y0)^2
i_min <- which.min(dist2)
xp <- x_ellipse[i_min]
yp <- y_ellipse[i_min]

dist3 <- (x_ellipse - x1)^2 + (y_ellipse - y1)^2
i_min2 <- which.min(dist3)
xp2 <- x_ellipse[i_min2]
yp2 <- y_ellipse[i_min2]

# -----------------------------
# 4) Plot
# -----------------------------
op <- par(mar = c(3.5, 3.5, 3, 1))
pdf("suspicious_vs_innocuous.pdf", width = 8, height = 6)

plot(NA, NA,
     xlim = c(-3, 6), ylim = c(-2.5, 3.2),
     xlab = "", ylab = "", axes = FALSE,
     main = "The Projection Test")

# Fill the convex set
polygon(x_ellipse, y_ellipse,
        col = adjustcolor("skyblue2", alpha.f = 0.35),
        border = "steelblue4", lwd = 2)

# Dashed line from outside point to projection
segments(x0, y0, xp, yp, lty = 2, lwd = 2, col = "gray30")
#text((x0 + xp)/2-0.6, (y0 + yp)/2, "far", pos = 3, col = "gray30")

# Outside point
points(x0, y0, pch = 19, cex = 1.4, col = "firebrick3")

# Projection point
points(xp, yp, pch = 19, cex = 1.4, col = "navy")

# Labels
text(x0 + 0.15, y0 + 0.1, "Suspicous", col = "firebrick3", adj = 0)
#text(xp + 0.15, yp - 0.1, "Closest clean t-curve", col = "navy", adj = 0)
text(-0.2, 0, "All Possible Clean t-curves", col = "steelblue4")

# Optional: show center
#points(cx, cy, pch = 16, cex = 0.8, col = "steelblue4")

#second point
points(x1, y1, pch = 19, cex = 1.4, col = "firebrick3")
text(x1 + 0.15, y0 + 0.1, "Innocuous", col = "firebrick3", adj = 0)
points(xp2, yp2, pch = 19, cex = 1.4, col = "navy")
segments(x1, y1, xp2, yp2, lty = 2, lwd = 2, col = "gray30")
#text((x1 + xp2)/2+0.6, (y1 + yp2)/2, "close", pos = 3, col = "gray30")



box()
par(op)

dev.off()