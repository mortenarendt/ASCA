data(oecdpanel, package = "np")

# Bandwidth by CV for local constant -- use only two starts to reduce the
# computation time
out <- capture.output(
  bwOECD <- np::npregbw(formula = growth ~ factor(oecd) + ordered(year) +
                          initgdp + popgro + inv + humancap, data = oecdpanel,
                        regtype = "lc", nmulti = 2)
)
bwOECD

# Regression
fitOECD <- np::npreg(bwOECD)
summary(fitOECD)

# Plot marginal effects of each predictor on the response
par(mfrow = c(2, 3))
plot(fitOECD, plot.par.mfrow = FALSE)

par(mfrow = c(2, 3))
np::npplot(fitOECD, plot.par.mfrow = FALSE)


