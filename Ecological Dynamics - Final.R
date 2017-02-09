rm(list=ls())
#install.packages("animation")
library(animation)
library(fields)
library(deSolve)
#install.packages("hash")
library(hash)

# ---- Periodic boundary conditions for neighbors ----
periodic <- function (x, dimension = 128) {
  if (x > dimension)
    x <- x - dimension
  else if (x < 1)
    x <- dimension + x
  return (x)  
}
# ---- Get 1-D position from 2-D coordinates ----
pos.from.coords <- function (x, y, dimension = 128) {
  pos <- y + (x-1) * dimension
  return(pos)
}
# ---- Get 2-D coordinates from 1-D position ----
coords.from.pos <- function (pos, dimension = 128) {
  x <- (pos - 1) %/% dimension + 1
  y <- (pos - 1) %% dimension + 1
  return(list(y=y, x=x))
}
# ---- Find neighbor coordinates ----
find.neighs4 <- function (x, y, dimension = 128) {
  xplus <- periodic(x+1, dimension)
  yplus <- periodic(y+1, dimension)
  xminus <- periodic(x-1, dimension)
  yminus <- periodic(y-1, dimension)
  
  # Return coordinates of neighbors
  return(list(x=c(x, x, xplus, xminus),
              y=c(yplus, yminus, y, y)))
}

# ---- Moore neighbor coordinates ----
find.neighs8 <- function (x, y, dimension = 128) {
  xplus <- periodic(x + 1, dimension)
  yplus <- periodic(y + 1, dimension)
  xminus <- periodic(x - 1, dimension)
  yminus <- periodic(y - 1, dimension)
  
  # Return coordinates of neighbors
  return(list(x=c(x, x, xplus, xminus, xminus, xminus, xplus, xplus),
              y=c(yplus, yminus, y, y, yminus, yplus, yplus, yminus)))
}

# ---- Moore neighbor positions ----
find.pos.neighs8 <- function (x, y, dimension = 128) {
  xplus <- periodic(x + 1, dimension)
  yplus <- periodic(y + 1, dimension)
  xminus <- periodic(x - 1, dimension)
  yminus <- periodic(y - 1, dimension)
  
  positions <- c(pos.from.coords(x = x, y = yplus, dimension = dimension),
                 pos.from.coords(x = x, y = yminus, dimension = dimension),
                 pos.from.coords(x = xplus, y = y, dimension = dimension),
                 pos.from.coords(x = xminus, y = y, dimension = dimension),
                 pos.from.coords(x = xminus, y = yminus, dimension = dimension),
                 pos.from.coords(x = xminus, y = yplus, dimension = dimension),
                 pos.from.coords(x = xplus, y = yplus, dimension = dimension),
                 pos.from.coords(x = xplus, y = yminus, dimension = dimension))
  
  # Return positions of neighbors
  return(positions)
}

# ---- Find neighbor positions ----
find.pos.neighs4 <- function (x, y, dimension = 128) {
  xplus <- periodic(x+1, dimension)
  yplus <- periodic(y+1, dimension)
  xminus <- periodic(x-1, dimension)
  yminus <- periodic(y-1, dimension)
  
  positions <- c(pos.from.coords(x = x, y = yplus, dimension = dimension),
                 pos.from.coords(x = x, y = yminus, dimension = dimension),
                 pos.from.coords(x = xplus, y = y, dimension = dimension),
                 pos.from.coords(x = xminus, y = y, dimension = dimension))
  
  # Return positions of neighbors
  return(positions)
}

# ---- Set neighbors for each cell ----
set.neighbors <- function (dimension = 128, size = dimension^2, n.neighs = 8) {
  neighs <- matrix(nrow = size, ncol = n.neighs + 1)
  for (p in 1:size) {
    # Find coordinates for position p
    p.coords <- coords.from.pos(p, dimension = dimension)
    # Find neighbors for position p
    all.neighs <- find.pos.neighs8(y = p.coords$y, 
                                   x = p.coords$x, 
                                   dimension = dimension)
    neighs[p, ] <- c(p, all.neighs)
  }
  return (neighs)
}

# ---- rules-IBM ----
rules.pred.prey <- function(lattice, a, parms) {
  # Empty = 1, # Prey = 2 - 101, # Pred = 102
  updated.lattice <- lattice
  with(parms, {
    rand.locs <- sample(1:size)
    for (i in rand.locs) {
      if (lattice[i] == 1) {
        # Empty
        colonizer <- strongest.colonizer(i, lattice, parms) # c(Type, Strength)
        if (colonizer[2] > runif(1)) {
          updated.lattice[i] <- colonizer[1]
        }
      } else if (2 <= lattice[i] & lattice[i] <= 101) {
        # Prey
        if (m > runif(1)) {
          updated.lattice[i] <- 1
        } else if (sum(lattice[neighs[i, 2:9]] == 102) > 0) {
          if (a > d.vec[lattice[i] - 1]) {
            updated.lattice[i] <- 102
          }
        }
      } else if (lattice[i] == 102) {
        # Predator
        if (m > runif(1)) {
          updated.lattice[i] <- 1
        }
      }
    }
    return(updated.lattice)
  })
}

strongest.colonizer <- function(i, lattice, parms) {
  n.freq <- hash() # Type -> Freq
  n.str <- hash() # Type -> Strength
  with(parms, {
    neigh.types <- as.character(lattice[neighs[i, 2:9]])
    # Set neighbor type frequencies
    for (n in neigh.types) {
      if (has.key(n, n.freq)) {
        .set(n.freq, keys=n, values=get(n, n.freq) + 1)
      } else {
        .set(n.freq, keys=n, values=1)
      }
    }
    # Determine strongest neighbor
    for (n in keys(n.freq)) {
      str <- r.vec[as.numeric(n)] * 0.125 * get(n, n.freq)
      .set(n.str, keys=n, values=str)
    }
    max.str <- max(values(n.str))
    n.str.inv <- invert(n.str)
    strongest <- sample(as.numeric(get(as.character(max.str), n.str.inv)), 1)
    return(c(strongest, max.str))
  })
}

# ---- theta-logistic curve ----
th.log <- function(x, theta = 1) {
  return(1 - x^theta)
}

# ---- relative abundances ----
rel.abundances <- function(lattice, t, size) {
  r.a.vec <- numeric(102)
  for (p in 1:102) {
    r.a.vec[p] <- sum(lattice[, , t] == p) / size * 100
  }
  return(r.a.vec)
}

# ---- Model parameters ----
dimension <- 64                         # Lattice dimension (1D)
size <- dimension^2                     # Lattice size (2D)
timeval <- 350                          # Total time of simulation
d.vec <- seq(0, .99, length.out = 100)  # Defensive abilities of prey
theta.vec <- c(.25, 1, 4)               # Theta to determine corresponding reproductive abilities
m <- 0.15                               # Mortality rate of prey/predator
a.vec <- seq(.1, 1, length.out = 10)    # Predation strength of predator

# ---- Initialize arrays ---
set.seed(1)
neighs <- set.neighbors(dimension = dimension)
lattice.1 <- array(sample(size = size, x = 1:102, prob = c(.40, rep(.005, 100), .10), replace = TRUE), 
                   dim = c(dimension, dimension, timeval))
lattice.2 <- array(sample(size = size, x = 1:102, prob = c(.40, rep(.005, 100), .10), replace = TRUE), 
                   dim = c(dimension, dimension, timeval))
lattice.3 <- array(sample(size = size, x = 1:102, prob = c(.40, rep(.005, 100), .10), replace = TRUE), 
                   dim = c(dimension, dimension, timeval))

rel.abun.1 <- array(dim = c(timeval, length(a.vec), 102))
rel.abun.2 <- array(dim = c(timeval, length(a.vec), 102))
rel.abun.3 <- array(dim = c(timeval, length(a.vec), 102))
for (a.i in 1:length(a.vec)) {
  rel.abun.1[1, a.i, ] <- rel.abundances(lattice.1, 1, size)
  rel.abun.2[1, a.i, ] <- rel.abundances(lattice.2, 1, size)
  rel.abun.3[1, a.i, ] <- rel.abundances(lattice.3, 1, size)
}

parms.1 <- list(d.vec = d.vec, r.vec = c(0, th.log(d.vec, theta.vec[1]), 0), m = m,
                dimension = dimension, size = size, neighs = neighs)
parms.2 <- list(d.vec = d.vec, r.vec = c(0, th.log(d.vec, theta.vec[2]), 0), m = m,
                dimension = dimension, size = size, neighs = neighs)
parms.3 <- list(d.vec = d.vec, r.vec = c(0, th.log(d.vec, theta.vec[3]), 0), m = m,
                dimension = dimension, size = size, neighs = neighs)

# ---- Launch simulation ---
start.time <- Sys.time()
print(start.time)
for (t in 2:timeval) {
  for (a.i in 1:length(a.vec)) {
    a = a.vec[a.i]
    lattice.1[ , , t] <- rules.pred.prey(lattice.1[ , , t - 1], a, parms = parms.1)
    lattice.2[ , , t] <- rules.pred.prey(lattice.2[ , , t - 1], a, parms = parms.2)
    lattice.3[ , , t] <- rules.pred.prey(lattice.3[ , , t - 1], a, parms = parms.3)
    rel.abun.1[t, a.i, ] <- rel.abundances(lattice.1, t, size)
    rel.abun.2[t, a.i, ] <- rel.abundances(lattice.2, t, size)
    rel.abun.3[t, a.i, ] <- rel.abundances(lattice.3, t, size)
  }
}
end.time <- Sys.time()
print(end.time - start.time)

# ---- Plot maps ----
make.plot <- function(lattice, dimension = 128, tmin = 1, tmax = timeval, 
                      max.val = 3, min.val = 1) {
  par(pty="s", mfrow = c(2, 1), tck = 0.01, mgp = c(1.5, 0.2, 0), 
      las = 0, mar = c(1, 1, 2, 1), oma = c(2, 1, 0, 0))
  for (t in tmin:tmax) {
    cat("Encoding time step: ", t, "\n")
    image(x = 1:dimension, y = 1:dimension, z = lattice[ , , t],  
          col = .cb[c("blue", "red", "black")], xpd = NA,
          xlab = "x-coordinates", ylab = "y-coordinates", 
          zlim = c(1, 3),
          main=paste0("SIR dynamics (time = ", t, ")"))
    
    if (t == tmin) {
      S <- sum(lattice[ , , t] == 1)
      I <- sum(lattice[ , , t] == 2)
      R <- sum(lattice[ , , t] == 3)
    } else {
      S <- apply(lattice[ , , tmin:t] == 1, MARGIN = 3, sum)
      I <- apply(lattice[ , , tmin:t] == 2, MARGIN = 3, sum)
      R <- apply(lattice[ , , tmin:t] == 3, MARGIN = 3, sum)
    }
    
    plot(tmin:t, S/(S + I + R), t ="l", xpd = NA,
         xlab = "Time", ylab = "Relative abundance", 
         lwd = 2, col = .cb["blue"], ylim = c(0, 1.1), xlim = c(tmin, tmax))
    lines(tmin:t, I/(S + I + R), lwd = 2, col = .cb["red"])
    lines(tmin:t, R/(S + I + R), lwd = 2, col = .cb["black"])
    legend(x = "topright", legend = c("S", "I", "R"), lty = 1, 
           col = .cb[c("blue", "red", "black")], bty = "n", lwd = 2)
  }
}

# ---- Setup animation ----
height <- 800
width <- 400
oopt <- ani.options(interval = 0, nmax = length(1:timeval), outdir = getwd(),
                    ani.dev = "png", ani.type = "png", autobrowse = FALSE,
                    ani.width = width, ani.height = height,
                    ffmpeg = paste0(getwd(), "/ffmpeg"))
saveVideo(make.plot(lattice), ffmpeg = "ffmpeg", 
          video.name = paste0("SIR_beta", betaval, ".mp4"),
          img.name = "SIR", interval = 0.1, width = width,
          height = height, clean = TRUE) #, other.opts = "-b 300k")
print(end.time - start.time)
