library(PtProcess)
library(Matrix)
library(lme4)
setwd("/home/jose/Documents/IPD440_Project/real/seismicHazardML/notebooks")
chile = read.csv("export_real.csv")
colnames(chile) <- c("Year", "Month", "Day", "Hour", "Minute", "Second", "Latitude", "Longitude", "Depth", "magnitude", "time")
chile$magnitude = chile$magnitude - 4.95

dmagn_mark2 <- function(x, data, params) {
    if (params[7] > 0) {
        lambda <- etas_gif(data, x[, "time"], params = params[1:5])
        y <- dgamma(x[, "magnitude"], shape = 1 + sqrt(lambda) * params[7],
                      rate = params[6], log = TRUE)
        } else y <- dexp(x[, "magnitude"], rate = params[6], log = TRUE)
        return(y)
}

rmagn_mark2 <- function(ti, data, params){
    if (params[7] > 0){
        lambda <- etas_gif(data, ti, params = params[1:5])
        y <- rgamma(1, shape = 1 + sqrt(lambda) * params[7], rate = params[6])
    } else y <- rexp(1, rate = params[6])
    return(list(Magnitude = y))
}

expmap <- function(y, p){
	y$params[1:5] <- exp(p)
	return(y)
}

allmap <- function(y, p){
	y$params <- exp(p)
	return(y)
}

# Null-model (Independent from the past)
TT <- c(0, as.numeric(tail(chile, n = 1)[11]))
params <- c(0.05, 3.1, 1.3, 0.02, 1.1, 1/mean(chile$magnitude), 0)
y <- mpp(data = chile, gif = etas_gif, mark = list(dmagn_mark2, rmagn_mark2), params = params, TT = TT, gmap = expression(params[1:5]), mmap = expression(params))
initial <- log(params[1:5])
w <- optim(initial, neglogLik, object = y, pmap = expmap, control = list(trace = 1, maxit = 100))
initial <- w$par
w <- nlm(neglogLik, initial, object = y, pmap = expmap, print.level = 2, iterlim = 500, typsize = initial)
x00 <- expmap(y, w$estimate)

# Full model (dependent of past)
#TT2 <- c(0, tail(chile, n = 1)[11])
initial <- log(c(0.05, 3.1, 13, 0.02, 1.1, 1/mean(chile$magnitude), 0.1))
z <- optim(initial, neglogLik, object = y, pmap = allmap, control = list(trace = 1, maxit = 200))
initial2 <- z$par
z <- nlm(neglogLik, initial, object = y, pmap = allmap, print.level = 2, iterlim = 500, typsize = initial)
x11 <- allmap(y, z$estimate)

# Store variables x00 and x11 which contain the GIF data (null-model and full-model respectibly)
saveRDS(x00, "null_model_mags_greater_than_5.rds")
saveRDS(x11, "full_model_mags_greater_than_5.rds")

# Plot GIF
par(mfrow=c(1,2)) # 1 row 2 columns
plot(x00, log = TRUE)
plot(x11, log = TRUE)
