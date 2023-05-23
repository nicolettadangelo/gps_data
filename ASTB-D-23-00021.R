rm(list = ls())

library(sf)
library(dplyr)
library(spatstat)
library(mgcv)

load("Palermo.RData")
# Palermo.RData contains:
# 1. df: the data (dataframe)
# 2. net: the road network of Palermo (linnet)
# 3. places: the observed point pattern of the turistic attractions (ppp)

# CREATE DATA FOR MODEL
dati.P <- lpp(X = cbind(df$x, df$y, as.factor(df$id)), net)
dati.P$data$marks <- as.factor(dati.P$data$marks)
df$x <- dati.P$data$x
df$y <- dati.P$data$y

P <- linequad(unmark(dati.P), nd = 200)
P # 278 data points, 10798 dummy points
W <- w.quad(P)
iz <- is.data(P)
Wdat <- W[iz]
Wdum <- W[!iz]

dati.d <- P$data
dummy <- P$dummy

## find the set of all possible marks
data.marks <- marks(dati.P)
markset <- levels(data.marks)
nmarks <- length(markset)

## replicate dummy points, one copy for each possible mark
## -> dummy x {1,..,K}
dumdum <- cartesian(dummy, markset)
Wdumdum <- rep.int(Wdum, nmarks)
Idumdum <- rep.int(npoints(dati.P) + seq_len(npoints(dummy)), nmarks)

## also make dummy marked points at same locations as data points
## but with different marks
dumdat <- cartesian(unmark(dati.d), markset)
Wdumdat <- rep.int(Wdat, nmarks)
Mdumdat <- marks(dumdat)
Idumdat <- rep.int(1:npoints(dati.P), nmarks)
Mrepdat <- rep.int(data.marks, nmarks)
ok <- (Mdumdat != Mrepdat)
dumdat <- dumdat[ok,]
Wdumdat <- Wdumdat[ok]
Idumdat <- Idumdat[ok]

## combine the two dummy patterns
dumb <- superimpose(dumdum, dumdat, W = dummy$window, check = FALSE)
Wdumb <- c(Wdumdum, Wdumdat)
Idumb <- c(Idumdum, Idumdat)

# might take a while..
Q <- superimpose(dati.P, dumb)

# merge data and dummy points
z <- c(rep(1, length(Wdat)), rep(0, length(Wdumb)))
y <- z / c(Wdat, Wdumb)

r <- 100
nu <- npoints(Q)
nx <- npoints(dati.P)
H <- matrix(NA, nu, nx)
# might take a while..
pd <- crossdist(Q, dati.P) 
for(i in 1:nu){
  for(j in 1:nx){
    if (pd[i, j] <= r)
      H[i, j] <- (1 - (pd[i, j] / r) ^ 2) ^ 2
    else
      H[i, j] <- 0
  }
}
vi <- rowSums(H)
w <- c(Wdat, Wdumb)

df_mod <- cbind(y, w, Q$data$x, Q$data$y, Q$data$marks, vi)
colnames(df_mod) <- c("y.pois", "w", "lat", "long", "id", "vi")

df_mod <- as.data.frame(df_mod)

DATI1id <- subset(df, !duplicated(df$id))
V2 <- DATI1id$id
for( i in 1:dim(DATI1id)[1]) {
  DATI1id$id[DATI1id$id == V2[i]] <- c(1:length(unique(df_mod$id)))[i]
}
df_mod_full <- inner_join(df_mod, DATI1id, by = "id")

udumb <- unmark(Q)
# might take a while..
uL <- lpp(udumb, L = net)
places_lpp <- as.lpp(places$data$x, places$data$y, L = net) 
df_mod_full$places_dist <- nncross(uL, places_lpp, what = "dist")
df_mod_full$places_dist[which(df_mod_full$places_dist == Inf)] <- NA
df_mod_full <- as.data.frame(df_mod_full)
df_mod_full$id <- as.factor(df_mod_full$id)

mod <- bam(y.pois ~ vi
            + s(lat, long, bs = "tp", k = 29)
            + places_dist
            + s(id, bs = 're'),
            weights = w,
            family = poisson, 
            data = df_mod_full,
            discrete = TRUE)

summary(mod)
gam.vcomp(mod)
exp(mod$coefficients[1]) * volume(net)  
