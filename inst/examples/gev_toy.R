# Simulated monthly maxima with a seasonal location shift (240 months)
set.seed(11)
season <- factor(rep(c("AM", "JJ", "AS", "ON"), 60), levels = c("AM", "JJ", "AS", "ON"))
period <- factor(rep(c("pre", "post"), each = 120), levels = c("pre", "post"))
toy <- data.frame(season = season, period = period,
                  max_precip = extRemes::revd(240, loc = 17 + c(0, 4, 6, 2)[season],
                                              scale = 7, shape = 0.05))
