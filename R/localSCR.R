#' Function to define state-space for a spatially-explicit mark recapture
#' models
#'
#' Returns a list of a matrix or array object of grid
#' coordinates and an Extent object from the \code{raster} package as a state-space.
#'
#' @param X either a matrix or array representing the coordinates of traps in
#' UTMs. An array is used when traps are clustered over a survey area.
#' @param crs_ the UTM coordinate reference system (EPSG code) used for your
#' location provided as an integer (e.g., 32608 for WGS 84/UTM Zone 8N).
#' @param buff the distance (m or km) that the traps should be
#' buffered by as an integer. This is typically 3 times the sigma parameter.
#' @param res the grid cell resolution for the state-space.
#' @return a list of a matrix or array of grid coordinates \code{grid} and an
#' extent object \code{ext}. Note that a matrix object is returned if the coordinates
#' of the traps are a matrix (i.e., 2D), otherwise an array object is returned when trap
#' coordinates are in a 3D array.
#' @details This function supports spatial capture-recapture analysis by
#' creating two outputs that are used to define the state-space. If a habitat mask is
#' not used, then only the \code{Extent} object from the \code{raster} package is needed under a uniform state-space.
#' The matrix or array object can be used to develop a habitat mask in a uniform
#' state-space or as a discretized state-space.
#' @author Daniel Eacker
#' @seealso \code{\link[raster:extent]{raster::extent}}
#' @importFrom raster extent coordinates raster rasterFromXYZ ncell
#' @importFrom sf st_crs 
#' @examples
#' # simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' traps <- traps + runif(prod(dim(traps)),-20,20) # add some random noise to locations
#'
#' mysigma = 300 # simulate sigma of 300 m
#' mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#'
#' # create state-space
#' Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = 100)
#'
#' # make plot of grid and trap locations
#' par(mfrow=c(1,1))
#' plot(Grid$grid, pch=19)
#' points(traps, col="blue",pch=20)
#' @export
grid_classic <- function(X, crs_, buff, res){
  # need to determine if X is 2 or 3 dimensional (stop if not)
  if(length(dim(X))!=2 & length(dim(X))!=3){
    stop("Trapping grid must be only 2 or 3 dimensions")
  }
  if(res >= max(buff)){
    stop("Resolution should not be equal to or greater than buffer")
  }
  # for dim of length 2
  if(length(dim(X))==2){
    # get new extent
    ext <- raster::extent(c(apply(X, 2, min)-max(buff),apply(X, 2, max)+max(buff))[c(1,3,2,4)])
    rast <- raster::raster(ext = ext, crs = sf::st_crs(crs_)$crs, res = res)
    gridOut <- raster::coordinates(rast)
  }else
    # for dim of length 3 (assumes different sites are in dimension 3)
    if(length(dim(X))==3){
      ext_initial <- apply(X, 3, function(x) raster::extent(c(apply(x, 2, min)-max(buff),apply(x, 2, max)+max(buff))[c(1,3,2,4)]))
      centroids <- lapply(ext_initial, function(x) cbind(mean(c(x[2],x[1])),mean(c(x[4],x[3]))))
      x_lengths <- unlist(lapply(ext_initial, function(x) x[2]-x[1]))
      x_index <- which(x_lengths==max(x_lengths))[1]
      y_lengths <- unlist(lapply(ext_initial, function(x) x[4]-x[3]))
      y_index <- which(y_lengths==max(y_lengths))[1]
      ext_standard <-extent(c(-x_lengths[x_index]/2,x_lengths[x_index]/2,-y_lengths[y_index]/2,y_lengths[y_index]/2))
      coords <- raster::coordinates(raster(ext = ext_standard, crs = sf::st_crs(crs_)$crs, res = res))
      rast <- lapply(centroids, function(x) raster::rasterFromXYZ(cbind(coords[,1]+x[,1],coords[,2]+x[,2]),crs = sf::st_crs(crs_)$crs))
      gridOut <- lapply(rast, raster::coordinates)
      ext <- lapply(rast, function(x) raster::extent(x))
      gridOut <- array(unlist(gridOut),dim=c(raster::ncell(rast[[1]]),dim(X)[2:3])) # convert to array
    }
  return(list(grid = gridOut, ext = ext))
} # end 'grid_classic' function



#' Function to simulate basic spatially-explicit mark recapture data
#'
#' Returns a list of simulated data including the encounter
#' history, binary sex indicator, activity centers, and site identifier.
#'
#' @param X either a matrix or array object representing the coordinates of traps in
#' UTMs. An array is used when traps are clustered over a survey area.
#' @param ext an \code{Extent} object from the \code{raster} package. This is returned from \code{\link{grid_classic}}.
#' @param crs_ the UTM coordinate reference system (EPSG code) used for your
#' location provided as an integer (e.g., 32608 for WGS 84/UTM Zone 8N).
#' @param N simulated abundance as an integer.
#' @param sigma_ the scaling parameter of the bivariate normal kernel either
#' in meters or kilometers as an integer.
#' @param prop_sex the portion of females or males as a numeric value. This
#' will depend upon the indicator coding scheme used (e.g., females = 1 and
#' males = 0; then proportion of females in the simulation). Must be a numeric
#' value between 0 and 1. Note that 0 or 1 can be used if a non-sex-specific sigma is desired.
#' @param K the number of sampling occasions desired as an integer.
#' @param base_encounter the baseline encounter probability or rate as a numeric
#' value. Note that a probabilty is given for a \code{"binomial"} observation
#' distribution while a rate is given for a \code{"poisson"} distribution.
#' @param enc_dist either \code{"binomial"} or \code{"poisson"}. Default is
#' \code{"binomial"}.
#' @param hab_mask either \code{FALSE} (the default) or a matrix or arrary output from \code{\link{mask_polygon}}
#' or \code{\link{mask_raster}} functions.
#' @param setSeed the random number generater seed as an integer used in
#' simulations to obtain repeatable data simulations. Default is 500.
#' @return a list of a matrix or array of encounter histories \code{y}, a
#' vector or matrix of 0's and 1's for \code{sex}, a batrix of simulated activity
#' centers \code{s}, and when a 3-dimensional trap array is given, a vector
#' for the site identifier \code{site}.
#' @author Daniel Eacker
#' @details This function supports spatial capture-recapture (SCR) analysis by allowing for
#' easy simulation of data components used by nimble in Baysian SCR models. Note that the
#' output for the encounter histories \code{y} will be sorted by detected and not detected individuals.
#' @seealso \code{\link{grid_classic}}
#' @importFrom stats runif rbinom rpois
#' @importFrom scales rescale
#' @importFrom sf st_cast st_sfc st_multipoint st_distance
#' @examples
#' # simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' traps <- traps + runif(prod(dim(traps)),-20,20) # add some random noise to locations
#'
#' mysigma = 300 # simulate sigma of 300 m
#' mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#'
#' # Create grid and extent
#' Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = 100)
#'
#' # simulate SCR data
#' data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, sigma_ = c(300), prop_sex = 1, N = 200, K = 4, base_encounter = 0.25, enc_dist = "binomial", hab_mask = FALSE, setSeed = 50)
#'
#' # make simple plot
#' par(mfrow=c(1,1))
#' plot(Grid$grid, pch=20,ylab="Northing",xlab="Easting")
#' points(traps, col="blue",pch=20)
#' points(data3d$s,col="red",pch = 20)
#' points(data3d$s[which(apply(data3d$y,1,sum)!=0),],col="green",pch = 20)
#' @export
sim_classic <- function(X, ext, crs_, N, sigma_, prop_sex, K, base_encounter, enc_dist = "binomial", hab_mask = FALSE, setSeed = 500){
  if(is.null(setSeed)==FALSE){
    set.seed(setSeed)
  }
  if(length(dim(X))!=2 & length(dim(X))!=3){
    stop("Trapping grid must be only 2 or 3 dimensions")
  }
  if(length(sigma_)>2 | length(sigma_)<1){
    stop("Sigma can only be of length 1 or 2")
  }
  if(enc_dist != "poisson" & enc_dist != "binomial"){
    stop("Encounter model has to be either binomial or poisson")
  }
  # for dim of length 2
  if(length(dim(X))==2){
    xlim <- ext[1:2] # create x limits for state-space
    ylim <- ext[3:4] # create y limits for state-space
    if(isFALSE(hab_mask)){
      sx <- runif(N,xlim[1],xlim[2])
      sy <- runif(N,ylim[1],ylim[2])
      s <- cbind(sx,sy)
    } else
      if(isFALSE(hab_mask)==FALSE){
        sx <- runif(N,xlim[1],xlim[2])
        sx.rescale <- scales::rescale(sx, to = c(0,dim(hab_mask)[2]), from=xlim)
        sy <- runif(N,ylim[1],ylim[2])
        sy.rescale <- scales::rescale(sy, to = c(0,dim(hab_mask)[1]), from=ylim)
        pOK <- numeric(N)
        for(i in 1:N){
          pOK[i] <- hab_mask[(trunc(sy.rescale[i])+1),(trunc(sx.rescale[i])+1)]
          while(pOK[i]==0){
            sx[i] <- runif(1,xlim[1],xlim[2])
            sx.rescale[i] <- scales::rescale(sx[i], to = c(0,dim(hab_mask)[2]), from=xlim)
            sy[i] <- runif(1,ylim[1],ylim[2])
            sy.rescale[i] <- scales::rescale(sy[i], to = c(0,dim(hab_mask)[1]), from=ylim)
            pOK[i] <- hab_mask[(trunc(sy.rescale[i])+1),(trunc(sx.rescale[i])+1)]
          }
        }
        s <- cbind(sx,sy)
      }
    sex <- rbinom(N, 1, prop_sex) # set propsex to 1 if you only want to simulate for one sex
    sex_ <- sex + 1 # sex indicator for sigma
    if(length(sigma_) == 1){ # make sigma length 2 if only 1 given
      sigma_ <- c(sigma_,sigma_)
    }
    Dmat <- sf::st_distance(sf::st_cast(sf::st_sfc(sf::st_multipoint(s), crs =  crs_),"POINT"),sf::st_cast(sf::st_sfc(sf::st_multipoint(X), crs =  crs_),"POINT"))  # compute distance matrix for traps
    Dmat <- matrix(as.numeric(Dmat), nrow=nrow(Dmat), ncol=ncol(Dmat)) # convert manually to matrix
    prob <- base_encounter*exp(-Dmat^2/(2*sigma_[sex_]^2)) # get detection prob matrix
    Y3d <- array(NA, dim=c(N,nrow(X),K)) # encounter data
    if(enc_dist == "binomial"){
      for(i in 1:N){
        for(j in 1:nrow(X)){
          Y3d[i,j,1:K] <- rbinom(K,1,prob[i,j])
        }}
    }else
      if(enc_dist == "poisson"){
        for(i in 1:N){
          for(j in 1:nrow(X)){
            Y3d[i,j,1:K] <- rpois(K,prob[i,j])
          }}
      }
    Y3d <- Y3d[c(which(apply(Y3d,1,sum)!=0),which(apply(Y3d,1,sum)==0)),,] # organize encountered indivdiuals and then 0's (augmented)
    s <- s[c(which(apply(Y3d,1,sum)!=0),which(apply(Y3d,1,sum)==0)),] # organize simulated activity centers
    sex[which(apply(Y3d,1,sum)==0)] <- NA
  }else
    # for dim of length 3 (assumes different sites are in dimension 3)
    if(length(dim(X))==3){
      if(length(sigma_) == 1){ # make sigma length 2 if only 1 given
        sigma_ <- c(sigma_,sigma_)
      }
      sex <- matrix(NA, nrow=N, ncol=dim(X)[3])
      site <- matrix(NA, nrow=N, ncol=dim(X)[3])
      smat <- array(NA, dim=c(N,2,dim(X)[3]))
      xlim <-sapply(ext, function(x) x[1:2]) # create x limits for state-space
      ylim <- sapply(ext, function(x) x[3:4]) # create y limits for state-space
      Y4d <- array(NA, dim=c(N,dim(X)[1],K,dim(X)[3]))
      for(g in 1:dim(X)[3]){ # now loop over trap dimensions
        sex[,g] <- rbinom(N, 1, prop_sex) # set propsex to 1 if you only want to simulate for one sex
        sex_ <- sex + 1 # indicator for sigma
        if(isFALSE(hab_mask)){
          sx <- runif(N,xlim[1,g],xlim[2,g])
          sy <- runif(N,ylim[1,g],ylim[2,g])
          s <- cbind(sx,sy)
        } else
          if(isFALSE(hab_mask)==FALSE){
            sx <- runif(N,xlim[1,g],xlim[2,g])
            sx.rescale <- scales::rescale(sx, to = c(0,dim(hab_mask)[2]), from=xlim[,g])
            sy <- runif(N,ylim[1,g],ylim[2,g])
            sy.rescale <- scales::rescale(sy, to = c(0,dim(hab_mask)[1]), from=ylim[,g])
            pOK <- numeric(N)
            for(i in 1:N){
              pOK[i] <- hab_mask[(trunc(sy.rescale[i])+1),(trunc(sx.rescale[i])+1),g]
              while(pOK[i]==0){
                sx[i] <- runif(1,xlim[1,g],xlim[2,g])
                sx.rescale[i] <- scales::rescale(sx[i], to = c(0,dim(hab_mask)[2]), from=xlim[,g])
                sy[i] <- runif(1,ylim[1,g],ylim[2,g])
                sy.rescale[i] <- scales::rescale(sy[i], to = c(0,dim(hab_mask)[1]), from=ylim[,g])
                pOK[i] <- hab_mask[(trunc(sy.rescale[i])+1),(trunc(sx.rescale[i])+1),g]
              }
            }
            s <- cbind(sx,sy)
          }
        Dmat <- st_distance(st_cast(st_sfc(st_multipoint(s), crs =  crs_),"POINT"),st_cast(st_sfc(st_multipoint(X[,,g]), crs =  crs_),"POINT"))  # compute distance matrix for traps
        Dmat <- matrix(as.numeric(Dmat), nrow=nrow(Dmat), ncol=ncol(Dmat)) # convert manually to matrix
        prob <- base_encounter*exp(-Dmat^2/(2*sigma_[sex_[,g]]^2)) # get detection prob matrix
        Y <- array(NA, dim=c(N,dim(X)[1],K)) # encounter data
        if(enc_dist == "binomial"){
          for(i in 1:N){
            for(j in 1:dim(X)[1]){
              Y[i,j,1:K] <- rbinom(K,1,prob[i,j])
            }}
        }else
          if(enc_dist == "poisson"){
            for(i in 1:N){
              for(j in 1:dim(X)[1]){
                Y[i,j,1:K] <- rpois(K,prob[i,j])
              }}
          }
        Y <- Y[c(which(apply(Y,1,sum)!=0),which(apply(Y,1,sum)==0)),,] # organize encountered and then 0's (augmented)
        Y4d[,,,g] <- Y
        smat[1:N,1:2,g] <- s[c(which(apply(Y,1,sum)!=0),which(apply(Y,1,sum)==0)),]
        sex[which(apply(Y,1,sum)==0),g] <- NA
        site[,g] <- g
      }
    }
  if(length(dim(X))==2){
    dataList <- list(y=Y3d,sex=sex,s=s)
  }else
    if(length(dim(X))==3){
      t.data <- data.frame(ind=1:length(as.vector(sex)),sex=as.vector(sex),site=as.vector(site))
      t.data <- t.data[c(which(is.na(t.data$sex)==FALSE),which(is.na(t.data$sex)==TRUE)),]
      dataList <- list(y=Y4d,sex=t.data$sex,site=t.data$site,s=smat)
    }
  return(dataList)
} # end 'sim_encounter' function



#' Function to retrieve nimbleCode for spatially-explicit mark recapture models
#'
#' Creates a \code{nimbleCode} object from the \code{nimble} package.
#'
#' @param dim_y an integer of either 2 (the default), 3, or 4 that defines what dimensional format the encounter history data are in.
#' @param enc_dist either \code{"binomial"} or \code{"poisson"}. Default is
#' \code{"binomial"}.
#' @param sex_sigma a logical value indicating whether the scaling parameter ('sigma') is sex-specific
#' @param hab_mask a logical value indicating whether a habitat mask will be used. Default is \code{FALSE}.
#' @return a \code{nimbleCode} object from the \code{nimble} package.
#' @details This function provides templates that could be copied and easily modified to include further model
#'  complexity such as covariates explaining detection probability. The models include different encounter probability distributions, sex-specific scaling parameters, and habitat masking.
#' @author Daniel Eacker
#' @importFrom nimble nimbleCode
#' @examples
#' # get model for 2D encounter data, binomial encounter distribution, non-sex-specific scaling parameter, and no habitat mask
#' scr_model = get_classic(dim_y = 2,enc_dist = "binomial",sex_sigma = FALSE,hab_mask = FALSE)
#'
#' # inspect model
#' scr_model
#' @export
get_classic <- function(dim_y, enc_dist = "binomial",sex_sigma = FALSE,hab_mask = FALSE){
    M <- J <- s <- X <- p0 <- sigma <- n0 <- z <- A <- lam0 <- K <- sex <- nSites <-
    site <- pixelWidth <- psi <- prop.habitat <- NULL
  if(dim_y !=2 & dim_y!=3 & dim_y !=4){
    stop("dim_y must be either 2, 3, or 4")
  }
  if(enc_dist != "poisson" & enc_dist != "binomial"){
    stop("Encounter distribution has to be either binomial or poisson")
  }
  if(isFALSE(hab_mask)){ # determine if hab_mask is included
    if(dim_y == 2){
      if(enc_dist == "binomial" & sex_sigma  == FALSE){
        scrcode <- nimble::nimbleCode({
          p0 ~ dunif(0,1) # baseline encounter probability
          sigma ~ dunif(0, sigma_upper) # scaling parameter
          psi ~ dunif(0, 1) # inclusion prob
          for(i in 1:M){
            z[i]~dbern(psi)
            s[i,1] ~ dunif(x_lower, x_upper)
            s[i,2] ~ dunif(y_lower, y_upper)
            dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2)
            p[i,1:J] <- p0*exp(-dist[i,1:J]^2/(2*sigma^2))
          } # i individuals
          # use zeros trick for detected individuals to speed up the computation
          for(i in 1:n0){
            for(j in 1:J){
              y[i,j] ~ dbin(p[i,j],K)
            } # i individuals
          } # j traps
          for(i in (n0+1):M){
            zeros[i] ~ dbern((1 - prod(1 - p[i,1:J])^K)*z[i])
          } # i individuals
          N <- sum(z[1:M])
          D <- N/A
        })
      } else
        if(enc_dist == "poisson" & sex_sigma  == FALSE){
          scrcode <- nimble::nimbleCode({
            lam0 ~ dunif(0,lam0_upper) # baseline encounter probability
            sigma ~ dunif(0, sigma_upper) # scaling parameter
            psi ~ dunif(0, 1) # inclusion prob
            for(i in 1:M){
              z[i]~dbern(psi)
              s[i,1] ~ dunif(x_lower, x_upper)
              s[i,2] ~ dunif(y_lower, y_upper)
              dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2)
              lam[i,1:J] <- lam0*K*exp(-dist[i,1:J]^2/(2*sigma^2))
            } # i individuals
            # use zeros trick for individuals to speed up the computation
            for(i in 1:n0){
              for(j in 1:J){
                y[i,j] ~ dpois(lam[i,j])
              } # i individuals
            } # j traps
            for(i in (n0+1):M){
              zeros[i] ~ dpois(sum(lam[i,1:J])*z[i])
            }
            N <- sum(z[1:M])
            D <- N/A
          })
        }else
          if(enc_dist == "binomial" & sex_sigma  == TRUE){
            scrcode <- nimble::nimbleCode({
              p0 ~ dunif(0,1) # baseline encounter probability
              psi_sex ~ dunif(0,1) # probability sex = 1
              sigma[1] ~ dunif(0, sigma_upper) # scaling parameter, sex = 0
              sigma[2] ~ dunif(0, sigma_upper) # scaling parameter, sex = 1
              psi ~ dunif(0, 1) # inclusion prob
              for(i in 1:M){
                z[i]~dbern(psi)
                sex[i] ~ dbern(psi_sex)
                sx[i] <- sex[i] + 1
                s[i,1] ~ dunif(x_lower, x_upper)
                s[i,2] ~ dunif(y_lower, y_upper)
                dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2)
                p[i,1:J] <- p0*exp(-dist[i,1:J]^2/(2*sigma[sx[i]]^2))
              } # i individuals
              # use zeros trick for individuals to speed up the computation
              for(i in 1:n0){
                for(j in 1:J){
                  y[i,j] ~ dbin(p[i,j],K)
                } # i individuals
              } # j traps
              for(i in (n0+1):M){
                zeros[i] ~ dbern((1 - prod(1 - p[i,1:J])^K)*z[i])
              } # i individuals
              N <- sum(z[1:M])
              D <- N/A
            })
          } else
            if(enc_dist == "poisson" & sex_sigma  == TRUE){
              scrcode <- nimble::nimbleCode({
                lam0 ~ dunif(0,lam0_upper) # baseline encounter probability
                psi_sex ~ dunif(0,1) # probability sex = 1
                sigma[1] ~ dunif(0, sigma_upper) # scaling parameter, sex = 0
                sigma[2] ~ dunif(0, sigma_upper) # scaling parameter, sex = 1
                psi ~ dunif(0, 1) # inclusion prob
                for(i in 1:M){
                  z[i]~dbern(psi)
                  sex[i] ~ dbern(psi_sex)
                  sx[i] <- sex[i] + 1
                  s[i,1] ~ dunif(x_lower, x_upper)
                  s[i,2] ~ dunif(y_lower, y_upper)
                  dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2)
                  lam[i,1:J] <- lam0*K*exp(-dist[i,1:J]^2/(2*sigma[sx[i]]^2))
                } # i individuals
                # use zeros trick for individuals to speed up the computation
                for(i in 1:n0){
                  for(j in 1:J){
                    y[i,j] ~ dpois(lam[i,j])
                  } # i individuals
                } # j traps
                for(i in (n0+1):M){
                  zeros[i] ~ dpois(sum(lam[i,1:J])*z[i])
                }
                N <- sum(z[1:M])
                D <- N/A
              })
            }
      return(scrcode)
    } else    # End 2D models
      if(dim_y == 3){
        if(enc_dist == "binomial" & sex_sigma  == FALSE){
          scrcode <- nimble::nimbleCode({
            p0 ~ dunif(0,1) # baseline encounter probability
            sigma ~ dunif(0, sigma_upper) # scaling parameter
            psi ~ dunif(0, 1) # inclusion prob
            for(i in 1:M){
              z[i]~dbern(psi)
              s[i,1] ~ dunif(x_lower, x_upper)
              s[i,2] ~ dunif(y_lower, y_upper)
              dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2)
              for(k in 1:K){
                p[i,1:J,k] <- p0*exp(-dist[i,1:J]^2/(2*sigma^2))
              }
            } # i individuals
            # use zeros trick for individuals to speed up the computation
            for(i in 1:n0){
              for(j in 1:J){
                for(k in 1:K){
                  y[i,j,k] ~ dbin(p[i,j,k],1)
                } # k occasions
              } # j traps
            } # i individuals
            for(i in (n0+1):M){
              zeros[i] ~ dbern((1 - prod(1 - p[i,1:J,1:K]))*z[i])
            } # i individuals
            N <- sum(z[1:M])
            D <- N/A
          })
        } else
          if(enc_dist == "poisson" & sex_sigma  == FALSE){
            scrcode <- nimble::nimbleCode({
              lam0 ~ dunif(0,lam0_upper) # baseline encounter rate
              sigma ~ dunif(0, sigma_upper) # scaling parameter
              psi ~ dunif(0, 1) # inclusion prob
              for(i in 1:M){
                z[i]~dbern(psi)
                s[i,1] ~ dunif(x_lower, x_upper)
                s[i,2] ~ dunif(y_lower, y_upper)
                dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2)
                for(k in 1:K){
                  lam[i,1:J,k] <- lam0*exp(-dist[i,1:J]^2/(2*sigma^2))
                } # k occasions
              } # i individuals
              # use zeros trick for individuals to speed up the computation
              for(i in 1:n0){
                for(j in 1:J){
                  for(k in 1:K){
                    y[i,j,k] ~ dpois(lam[i,j,k])
                  } # occasions
                } # j traps
              } # i individuals
              for(i in (n0+1):M){
                zeros[i] ~ dpois(sum(lam[i,1:J,1:K])*z[i])
              }
              N <- sum(z[1:M])
              D <- N/A
            })
          }else
            if(enc_dist == "binomial" & sex_sigma  == TRUE){
              scrcode <- nimble::nimbleCode({
                p0 ~ dunif(0,1) # baseline encounter probability
                psi_sex ~ dunif(0,1) # probability sex = 1
                sigma[1] ~ dunif(0, sigma_upper) # scaling parameter, sex = 0
                sigma[2] ~ dunif(0, sigma_upper) # scaling parameter, sex = 1
                psi ~ dunif(0, 1) # inclusion prob
                for(i in 1:M){
                  z[i]~dbern(psi)
                  sex[i] ~ dbern(psi_sex)
                  sx[i] <- sex[i] + 1
                  s[i,1] ~ dunif(x_lower, x_upper)
                  s[i,2] ~ dunif(y_lower, y_upper)
                  dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2)
                  for(k in 1:K){
                    p[i,1:J,k] <- p0*exp(-dist[i,1:J]^2/(2*sigma[sx[i]]^2))
                  } # k occasions
                } # i individuals
                # use zeros trick for individuals to speed up the computation
                for(i in 1:n0){
                  for(j in 1:J){
                    for(k in 1:K){
                      y[i,j,k] ~ dbin(p[i,j,k],1)
                    } # k occasions
                  } # j traps
                } # i individuals
                for(i in (n0+1):M){
                  zeros[i] ~ dbern((1 - prod(1 - p[i,1:J,1:K]))*z[i])
                } # i individuals
                N <- sum(z[1:M])
                D <- N/A
              })
            } else
              if(enc_dist == "poisson" & sex_sigma  == TRUE){
                scrcode <- nimble::nimbleCode({
                  lam0 ~ dunif(0,lam0_upper) # baseline encounter rate
                  psi_sex ~ dunif(0,1) # probability sex = 1
                  sigma[1] ~ dunif(0, sigma_upper) # scaling parameter, sex = 0
                  sigma[2] ~ dunif(0, sigma_upper) # scaling parameter, sex = 1
                  psi ~ dunif(0, 1) # inclusion prob
                  for(i in 1:M){
                    z[i]~dbern(psi)
                    sex[i] ~ dbern(psi_sex)
                    sx[i] <- sex[i] + 1
                    s[i,1] ~ dunif(x_lower, x_upper)
                    s[i,2] ~ dunif(y_lower, y_upper)
                    dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2)
                    for(k in 1:K){
                      lam[i,1:J,k] <- lam0*exp(-dist[i,1:J]^2/(2*sigma[sx[i]]^2))
                    } # k occasions
                  } # i marked individuals
                  # use zeros trick for marked individuals to speed up the computation
                  for(i in 1:n0){
                    for(j in 1:J){
                      for(k in 1:K){
                        y[i,j,k] ~ dpois(lam[i,j,k])
                      } # occasions
                    } # j traps
                  } # i individuals
                  for(i in (n0+1):M){
                    zeros[i] ~ dpois(sum(lam[i,1:J,1:K])*z[i])
                  }
                  N <- sum(z[1:M])
                  D <- N/A
                })
              }
        return(scrcode)
      } else  # End 3D models
        if(dim_y == 4){
          if(enc_dist == "binomial" & sex_sigma  == FALSE){
            scrcode <- nimble::nimbleCode({
              sigma ~ dunif(0, sigma_upper) # scaling parameter
              psi ~ dunif(0, 1) # inclusion prob
              for(g in 1:nSites){
                p0[g] ~ dunif(0,1) # site-specific baseline encounter probability
              } # g sites
              for(i in 1:M){
                z[i]~dbern(psi)
                s[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
                s[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
                dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1,site[i]])^2 + (s[i,2]-X[1:J,2,site[i]])^2)
                for(k in 1:K){
                  p[i,1:J,k] <- p0[site[i]]*exp(-dist[i,1:J]^2/(2*sigma^2))
                }
              } # i individuals
              # use zeros trick for marked individuals to speed up the computation
              for(g in 1:nSites){
                for(i in 1:n0){
                  for(j in 1:J){
                    for(k in 1:K){
                      y[i,j,k,g] ~ dbin(p[i,j,k],1)
                    } # k occasions
                  } # j traps
                } # i individuals
              } # g sites
              for(i in (n0+1):M){
                zeros[i] ~ dbern((1 - prod(1 - p[i,1:J,1:K]))*z[i])
              } # i individuals
              N <- sum(z[1:M])
              D <- N/A
            })
          } else
            if(enc_dist == "poisson" & sex_sigma  == FALSE){
              scrcode <- nimble::nimbleCode({
                sigma ~ dunif(0, sigma_upper) # scaling parameter
                psi ~ dunif(0, 1) # inclusion prob
                for(g in 1:nSites){
                  lam0[g] ~ dunif(0,lam0_upper) # site-specific baseline encounter rate
                } # g sites
                for(i in 1:M){
                  z[i]~dbern(psi)
                  s[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
                  s[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
                  dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1,site[i]])^2 + (s[i,2]-X[1:J,2,site[i]])^2)
                  for(k in 1:K){
                    lam[i,1:J,k] <- lam0[site[i]]*exp(-dist[i,1:J]^2/(2*sigma^2))
                  } # k occasions
                } # i individuals
                # use zeros trick for marked individuals to speed up the computation
                for(g in 1:nSites){
                  for(i in 1:n0){
                    for(j in 1:J){
                      for(k in 1:K){
                        y[i,j,k,g] ~ dpois(lam[i,j,k])
                      } # occasions
                    } # j traps
                  } # i individuals
                } # g sites
                for(i in (n0+1):M){
                  zeros[i] ~ dpois(sum(lam[i,1:J,1:K])*z[i])
                }
                N <- sum(z[1:M])
                D <- N/A
              })
            }else
              if(enc_dist == "binomial" & sex_sigma  == TRUE){
                scrcode <- nimble::nimbleCode({
                  psi_sex ~ dunif(0,1) # probability sex = 1
                  sigma[1] ~ dunif(0, sigma_upper) # scaling parameter, sex = 0
                  sigma[2] ~ dunif(0, sigma_upper) # scaling parameter, sex = 1
                  psi ~ dunif(0, 1) # inclusion prob
                  for(g in 1:nSites){
                    p0[g] ~ dunif(0,1) # site-specific baseline encounter probability
                  } # g sites
                  for(i in 1:M){
                    z[i]~dbern(psi)
                    sex[i] ~ dbern(psi_sex)
                    sx[i] <- sex[i] + 1
                    s[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
                    s[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
                    dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1,site[i]])^2 + (s[i,2]-X[1:J,2,site[i]])^2)
                    for(k in 1:K){
                      p[i,1:J,k] <- p0[site[i]]*exp(-dist[i,1:J]^2/(2*sigma[sx[i]]^2))
                    } # k occasions
                  } # i marked individuals
                  # use zeros trick for marked individuals to speed up the computation
                  for(g in 1:nSites){
                    for(i in 1:n0){
                      for(j in 1:J){
                        for(k in 1:K){
                          y[i,j,k,g] ~ dbin(p[i,j,k],1)
                        } # k occasions
                      } # j traps
                    } # i individuals
                  } # g sites
                  for(i in (n0+1):M){
                    zeros[i] ~ dbern((1 - prod(1 - p[i,1:J,1:K]))*z[i])
                  } # i individuals
                  N <- sum(z[1:M])
                  D <- N/A
                })
              } else
                if(enc_dist == "poisson" & sex_sigma  == TRUE){
                  scrcode <- nimble::nimbleCode({
                    lam0 ~ dunif(0,lam0_upper) # baseline encounter rate
                    psi_sex ~ dunif(0,1) # probability sex = 1
                    sigma[1] ~ dunif(0, sigma_upper) # scaling parameter, sex = 0
                    sigma[2] ~ dunif(0, sigma_upper) # scaling parameter, sex = 1
                    psi ~ dunif(0, 1) # inclusion prob
                    for(i in 1:M){
                      z[i]~dbern(psi)
                      sex[i] ~ dbern(psi_sex)
                      sx[i] <- sex[i] + 1
                      s[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
                      s[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
                      dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1,site[i]])^2 + (s[i,2]-X[1:J,2,site[i]])^2)
                      for(k in 1:K){
                        lam[i,1:J,k] <- lam0*exp(-dist[i,1:J]^2/(2*sigma[sx[i]]^2))
                      } # k occasions
                    } # i individuals
                    # use zeros trick for marked individuals to speed up the computation
                    for(g in 1:nSites){
                      for(i in 1:n0){
                        for(j in 1:J){
                          for(k in 1:K){
                            y[i,j,k,g] ~ dpois(lam[i,j,k])
                          } # occasions
                        } # j traps
                      } # i individuals
                    } # g sites
                    for(i in (n0+1):M){
                      zeros[i] ~ dpois(sum(lam[i,1:J,1:K])*z[i])
                    }
                    N <- sum(z[1:M])
                    D <- N/A
                  })
                }
          return(scrcode)
        } # end 4D model
  } else
    if(hab_mask==TRUE){
      if(dim_y == 2){
        if(enc_dist == "binomial" & sex_sigma  == FALSE){
          scrcode <- nimble::nimbleCode({
            p0 ~ dunif(0,1) # baseline encounter probability
            sigma ~ dunif(0, sigma_upper) # scaling parameter
            sigma.pixel <- sigma / pixelWidth # scaled for habitat mask
            psi ~ dunif(0, 1) # inclusion prob
            for(i in 1:M){
              z[i]~dbern(psi)
              s[i,1] ~ dunif(x_lower, x_upper)
              s[i,2] ~ dunif(y_lower, y_upper)
              pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1)] # habitat check
              OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
              dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2)
              p[i,1:J] <- p0*exp(-dist[i,1:J]^2/(2*sigma.pixel^2))
            } # i individuals
            # use zeros trick for detected individuals to speed up the computation
            for(i in 1:n0){
              for(j in 1:J){
                y[i,j] ~ dbin(p[i,j],K)
              } # i individuals
            } # j traps
            for(i in (n0+1):M){
              zeros[i] ~ dbern((1 - prod(1 - p[i,1:J])^K)*z[i])
            } # i individuals
            N <- sum(z[1:M])
            D <- N/A
          })
        } else
          if(enc_dist == "poisson" & sex_sigma  == FALSE){
            scrcode <- nimble::nimbleCode({
              lam0 ~ dunif(0,lam0_upper) # baseline encounter probability
              sigma ~ dunif(0, sigma_upper) # scaling parameter
              sigma.pixel <- sigma / pixelWidth # scaled for habitat mask
              psi ~ dunif(0, 1) # inclusion prob
              for(i in 1:M){
                z[i]~dbern(psi)
                s[i,1] ~ dunif(x_lower, x_upper)
                s[i,2] ~ dunif(y_lower, y_upper)
                pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1)] # habitat check
                OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
                dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2)
                lam[i,1:J] <- lam0*K*exp(-dist[i,1:J]^2/(2*sigma.pixel^2))
              } # i individuals
              # use zeros trick for individuals to speed up the computation
              for(i in 1:n0){
                for(j in 1:J){
                  y[i,j] ~ dpois(lam[i,j])
                } # i individuals
              } # j traps
              for(i in (n0+1):M){
                zeros[i] ~ dpois(sum(lam[i,1:J])*z[i])
              }
              N <- sum(z[1:M])
              D <- N/A
            })
          }else
            if(enc_dist == "binomial" & sex_sigma  == TRUE){
              scrcode <- nimble::nimbleCode({
                p0 ~ dunif(0,1) # baseline encounter probability
                psi_sex ~ dunif(0,1) # probability sex = 1
                sigma[1] ~ dunif(0, sigma_upper) # scaling parameter, sex = 0
                sigma[2] ~ dunif(0, sigma_upper) # scaling parameter, sex = 1
                sigma.pixel[1] <- sigma[1] / pixelWidth # scaled for habitat mask
                sigma.pixel[2] <- sigma[2] / pixelWidth # scaled for habitat mask
                psi ~ dunif(0, 1) # inclusion prob
                for(i in 1:M){
                  z[i]~dbern(psi)
                  sex[i] ~ dbern(psi_sex)
                  sx[i] <- sex[i] + 1
                  s[i,1] ~ dunif(x_lower, x_upper)
                  s[i,2] ~ dunif(y_lower, y_upper)
                  pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1)] # habitat check
                  OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
                  dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2)
                  p[i,1:J] <- p0*exp(-dist[i,1:J]^2/(2*sigma.pixel[sx[i]]^2))
                } # i individuals
                # use zeros trick for individuals to speed up the computation
                for(i in 1:n0){
                  for(j in 1:J){
                    y[i,j] ~ dbin(p[i,j],K)
                  } # i individuals
                } # j traps
                for(i in (n0+1):M){
                  zeros[i] ~ dbern((1 - prod(1 - p[i,1:J])^K)*z[i])
                } # i individuals
                N <- sum(z[1:M])
                D <- N/A
              })
            } else
              if(enc_dist == "poisson" & sex_sigma  == TRUE){
                scrcode <- nimble::nimbleCode({
                  lam0 ~ dunif(0,lam0_upper) # baseline encounter probability
                  psi_sex ~ dunif(0,1) # probability sex = 1
                  sigma[1] ~ dunif(0, sigma_upper) # scaling parameter, sex = 0
                  sigma[2] ~ dunif(0, sigma_upper) # scaling parameter, sex = 1
                  sigma.pixel[1] <- sigma[1] / pixelWidth # scaled for habitat mask
                  sigma.pixel[2] <- sigma[2] / pixelWidth # scaled for habitat mask
                  psi ~ dunif(0, 1) # inclusion prob
                  for(i in 1:M){
                    z[i]~dbern(psi)
                    sex[i] ~ dbern(psi_sex)
                    sx[i] <- sex[i] + 1
                    s[i,1] ~ dunif(x_lower, x_upper)
                    s[i,2] ~ dunif(y_lower, y_upper)
                    pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1)] # habitat check
                    OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
                    dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2)
                    lam[i,1:J] <- lam0*K*exp(-dist[i,1:J]^2/(2*sigma.pixel[sx[i]]^2))
                  } # i individuals
                  # use zeros trick for individuals to speed up the computation
                  for(i in 1:n0){
                    for(j in 1:J){
                      y[i,j] ~ dpois(lam[i,j])
                    } # i individuals
                  } # j traps
                  for(i in (n0+1):M){
                    zeros[i] ~ dpois(sum(lam[i,1:J])*z[i])
                  }
                  N <- sum(z[1:M])
                  D <- N/A
                })
              }
        return(scrcode)
      } else    # End 2D models
        if(dim_y == 3){
          if(enc_dist == "binomial" & sex_sigma  == FALSE){
            scrcode <- nimble::nimbleCode({
              p0 ~ dunif(0,1) # baseline encounter probability
              sigma ~ dunif(0, sigma_upper) # scaling parameter
              sigma.pixel <- sigma / pixelWidth # scaled for habitat mask
              psi ~ dunif(0, 1) # inclusion prob
              for(i in 1:M){
                z[i]~dbern(psi)
                s[i,1] ~ dunif(x_lower, x_upper)
                s[i,2] ~ dunif(y_lower, y_upper)
                pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1)] # habitat check
                OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
                dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2)
                for(k in 1:K){
                  p[i,1:J,k] <- p0*exp(-dist[i,1:J]^2/(2*sigma.pixel^2))
                }
              } # i individuals
              # use zeros trick for individuals to speed up the computation
              for(i in 1:n0){
                for(j in 1:J){
                  for(k in 1:K){
                    y[i,j,k] ~ dbin(p[i,j,k],1)
                  } # k occasions
                } # j traps
              } # i individuals
              for(i in (n0+1):M){
                zeros[i] ~ dbern((1 - prod(1 - p[i,1:J,1:K]))*z[i])
              } # i individuals
              N <- sum(z[1:M])
              D <- N/A
            })
          } else
            if(enc_dist == "poisson" & sex_sigma  == FALSE){
              scrcode <- nimble::nimbleCode({
                lam0 ~ dunif(0,lam0_upper) # baseline encounter rate
                sigma ~ dunif(0, sigma_upper) # scaling parameter
                sigma.pixel <- sigma / pixelWidth # scaled for habitat mask
                psi ~ dunif(0, 1) # inclusion prob
                for(i in 1:M){
                  z[i]~dbern(psi)
                  s[i,1] ~ dunif(x_lower, x_upper)
                  s[i,2] ~ dunif(y_lower, y_upper)
                  pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1)] # habitat check
                  OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
                  dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2)
                  for(k in 1:K){
                    lam[i,1:J,k] <- lam0*exp(-dist[i,1:J]^2/(2*sigma.pixel^2))
                  } # k occasions
                } # i individuals
                # use zeros trick for individuals to speed up the computation
                for(i in 1:n0){
                  for(j in 1:J){
                    for(k in 1:K){
                      y[i,j,k] ~ dpois(lam[i,j,k])
                    } # occasions
                  } # j traps
                } # i individuals
                for(i in (n0+1):M){
                  zeros[i] ~ dpois(sum(lam[i,1:J,1:K])*z[i])
                }
                N <- sum(z[1:M])
                D <- N/A
              })
            }else
              if(enc_dist == "binomial" & sex_sigma  == TRUE){
                scrcode <- nimble::nimbleCode({
                  p0 ~ dunif(0,1) # baseline encounter probability
                  psi_sex ~ dunif(0,1) # probability sex = 1
                  sigma[1] ~ dunif(0, sigma_upper) # scaling parameter, sex = 0
                  sigma[2] ~ dunif(0, sigma_upper) # scaling parameter, sex = 1
                  sigma.pixel[1] <- sigma[1] / pixelWidth # scaled for habitat mask
                  sigma.pixel[2] <- sigma[2] / pixelWidth # scaled for habitat mask
                  psi ~ dunif(0, 1) # inclusion prob
                  for(i in 1:M){
                    z[i]~dbern(psi)
                    sex[i] ~ dbern(psi_sex)
                    sx[i] <- sex[i] + 1
                    s[i,1] ~ dunif(x_lower, x_upper)
                    s[i,2] ~ dunif(y_lower, y_upper)
                    pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1)] # habitat check
                    OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
                    dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2)
                    for(k in 1:K){
                      p[i,1:J,k] <- p0*exp(-dist[i,1:J]^2/(2*sigma.pixel[sx[i]]^2))
                    } # k occasions
                  } # i individuals
                  # use zeros trick for individuals to speed up the computation
                  for(i in 1:n0){
                    for(j in 1:J){
                      for(k in 1:K){
                        y[i,j,k] ~ dbin(p[i,j,k],1)
                      } # k occasions
                    } # j traps
                  } # i individuals
                  for(i in (n0+1):M){
                    zeros[i] ~ dbern((1 - prod(1 - p[i,1:J,1:K]))*z[i])
                  } # i individuals
                  N <- sum(z[1:M])
                  D <- N/A
                })
              } else
                if(enc_dist == "poisson" & sex_sigma  == TRUE){
                  scrcode <- nimble::nimbleCode({
                    lam0 ~ dunif(0,lam0_upper) # baseline encounter rate
                    psi_sex ~ dunif(0,1) # probability sex = 1
                    sigma[1] ~ dunif(0, sigma_upper) # scaling parameter, sex = 0
                    sigma[2] ~ dunif(0, sigma_upper) # scaling parameter, sex = 1
                    sigma.pixel[1] <- sigma[1] / pixelWidth # scaled for habitat mask
                    sigma.pixel[2] <- sigma[2] / pixelWidth # scaled for habitat mask
                    psi ~ dunif(0, 1) # inclusion prob
                    for(i in 1:M){
                      z[i]~dbern(psi)
                      sex[i] ~ dbern(psi_sex)
                      sx[i] <- sex[i] + 1
                      s[i,1] ~ dunif(x_lower, x_upper)
                      s[i,2] ~ dunif(y_lower, y_upper)
                      pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1)] # habitat check
                      OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
                      dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1])^2 + (s[i,2]-X[1:J,2])^2)
                      for(k in 1:K){
                        lam[i,1:J,k] <- lam0*exp(-dist[i,1:J]^2/(2*sigma.pixel[sx[i]]^2))
                      } # k occasions
                    } # i marked individuals
                    # use zeros trick for marked individuals to speed up the computation
                    for(i in 1:n0){
                      for(j in 1:J){
                        for(k in 1:K){
                          y[i,j,k] ~ dpois(lam[i,j,k])
                        } # occasions
                      } # j traps
                    } # i individuals
                    for(i in (n0+1):M){
                      zeros[i] ~ dpois(sum(lam[i,1:J,1:K])*z[i])
                    }
                    N <- sum(z[1:M])
                    D <- N/A
                  })
                }
          return(scrcode)
        } else  # End 3D models
          if(dim_y == 4){
            if(enc_dist == "binomial" & sex_sigma  == FALSE){
              scrcode <- nimble::nimbleCode({
                sigma ~ dunif(0, sigma_upper) # scaling parameter
                sigma.pixel <- sigma / pixelWidth # scaled for habitat mask
                psi ~ dunif(0, 1) # inclusion prob
                for(g in 1:nSites){
                  p0[g] ~ dunif(0,1) # site-specific baseline encounter probability
                } # g sites
                for(i in 1:M){
                  z[i]~dbern(psim[i])
                  psim[i] <- (1-(1-psi)^prop.habitat[site[i]]) # adjust psi for the proportion of available habitat at each site
                  s[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
                  s[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
                  pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1),site[i]] # habitat check
                  OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
                  dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1,site[i]])^2 + (s[i,2]-X[1:J,2,site[i]])^2)
                  for(k in 1:K){
                    p[i,1:J,k] <- p0[site[i]]*exp(-dist[i,1:J]^2/(2*sigma.pixel^2))
                  }
                } # i individuals
                # use zeros trick for marked individuals to speed up the computation
                for(g in 1:nSites){
                  for(i in 1:n0){
                    for(j in 1:J){
                      for(k in 1:K){
                        y[i,j,k,g] ~ dbin(p[i,j,k],1)
                      } # k occasions
                    } # j traps
                  } # i individuals
                } # g sites
                for(i in (n0+1):M){
                  zeros[i] ~ dbern((1 - prod(1 - p[i,1:J,1:K]))*z[i])
                } # i individuals
                N <- sum(z[1:M])
                D <- N/A
              })
            } else
              if(enc_dist == "poisson" & sex_sigma  == FALSE){
                scrcode <- nimble::nimbleCode({
                  sigma ~ dunif(0, sigma_upper) # scaling parameter
                  sigma.pixel <- sigma / pixelWidth # scaled for habitat mask
                  psi ~ dunif(0, 1) # inclusion prob
                  for(g in 1:nSites){
                    lam0[g] ~ dunif(0,lam0_upper) # site-specific baseline encounter rate
                  } # g sites
                  for(i in 1:M){
                    z[i]~dbern(psim[i])
                    psim[i] <- (1-(1-psi)^prop.habitat[site[i]]) # adjust psi for the proportion of available habitat at each site
                    s[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
                    s[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
                    pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1),site[i]] # habitat check
                    OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
                    dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1,site[i]])^2 + (s[i,2]-X[1:J,2,site[i]])^2)
                    for(k in 1:K){
                      lam[i,1:J,k] <- lam0[site[i]]*exp(-dist[i,1:J]^2/(2*sigma.pixel^2))
                    } # k occasions
                  } # i individuals
                  # use zeros trick for marked individuals to speed up the computation
                  for(g in 1:nSites){
                    for(i in 1:n0){
                      for(j in 1:J){
                        for(k in 1:K){
                          y[i,j,k,g] ~ dpois(lam[i,j,k])
                        } # occasions
                      } # j traps
                    } # i individuals
                  } # g sites
                  for(i in (n0+1):M){
                    zeros[i] ~ dpois(sum(lam[i,1:J,1:K])*z[i])
                  }
                  N <- sum(z[1:M])
                  D <- N/A
                })
              }else
                if(enc_dist == "binomial" & sex_sigma  == TRUE){
                  scrcode <- nimble::nimbleCode({
                    psi_sex ~ dunif(0,1) # probability sex = 1
                    sigma[1] ~ dunif(0, sigma_upper) # scaling parameter, sex = 0
                    sigma[2] ~ dunif(0, sigma_upper) # scaling parameter, sex = 1
                    sigma.pixel[1] <- sigma[1] / pixelWidth # scaled for habitat mask
                    sigma.pixel[2] <- sigma[2] / pixelWidth # scaled for habitat mask
                    psi ~ dunif(0, 1) # inclusion prob
                    for(g in 1:nSites){
                      p0[g] ~ dunif(0,1) # site-specific baseline encounter probability
                    } # g sites
                    for(i in 1:M){
                      z[i]~dbern(psim[i])
                      psim[i] <- (1-(1-psi)^prop.habitat[site[i]]) # adjust psi for the proportion of available habitat at each site
                      sex[i] ~ dbern(psi_sex)
                      sx[i] <- sex[i] + 1
                      s[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
                      s[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
                      pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1),site[i]] # habitat check
                      OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
                      dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1,site[i]])^2 + (s[i,2]-X[1:J,2,site[i]])^2)
                      for(k in 1:K){
                        p[i,1:J,k] <- p0[site[i]]*exp(-dist[i,1:J]^2/(2*sigma.pixel[sx[i]]^2))
                      } # k occasions
                    } # i marked individuals
                    # use zeros trick for marked individuals to speed up the computation
                    for(g in 1:nSites){
                      for(i in 1:n0){
                        for(j in 1:J){
                          for(k in 1:K){
                            y[i,j,k,g] ~ dbin(p[i,j,k],1)
                          } # k occasions
                        } # j traps
                      } # i individuals
                    } # g sites
                    for(i in (n0+1):M){
                      zeros[i] ~ dbern((1 - prod(1 - p[i,1:J,1:K]))*z[i])
                    } # i individuals
                    N <- sum(z[1:M])
                    D <- N/A
                  })
                } else
                  if(enc_dist == "poisson" & sex_sigma  == TRUE){
                    scrcode <- nimble::nimbleCode({
                      lam0 ~ dunif(0,lam0_upper) # baseline encounter rate
                      psi_sex ~ dunif(0,1) # probability sex = 1
                      sigma[1] ~ dunif(0, sigma_upper) # scaling parameter, sex = 0
                      sigma[2] ~ dunif(0, sigma_upper) # scaling parameter, sex = 1
                      sigma.pixel[1] <- sigma[1] / pixelWidth # scaled for habitat mask
                      sigma.pixel[2] <- sigma[2] / pixelWidth # scaled for habitat mask
                      psi ~ dunif(0, 1) # inclusion prob
                      for(i in 1:M){
                        z[i]~dbern(psim[i])
                        psim[i] <- (1-(1-psi)^prop.habitat[site[i]]) # adjust psi for the proportion of available habitat at each site
                        sex[i] ~ dbern(psi_sex)
                        sx[i] <- sex[i] + 1
                        s[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
                        s[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
                        pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1),site[i]] # habitat check
                        OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
                        dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1,site[i]])^2 + (s[i,2]-X[1:J,2,site[i]])^2)
                        for(k in 1:K){
                          lam[i,1:J,k] <- lam0*exp(-dist[i,1:J]^2/(2*sigma.pixel[sx[i]]^2))
                        } # k occasions
                      } # i individuals
                      # use zeros trick for marked individuals to speed up the computation
                      for(g in 1:nSites){
                        for(i in 1:n0){
                          for(j in 1:J){
                            for(k in 1:K){
                              y[i,j,k,g] ~ dpois(lam[i,j,k])
                            } # occasions
                          } # j traps
                        } # i individuals
                      } # g sites
                      for(i in (n0+1):M){
                        zeros[i] ~ dpois(sum(lam[i,1:J,1:K])*z[i])
                      }
                      N <- sum(z[1:M])
                      D <- N/A
                    })
                  }
            return(scrcode)
          } # end 4D model
    } # end models with hab_mask
} # End function 'get_classic"



#' Function to generate starting locations for activity centers
#'
#' Generate a matrix of intial starting locations, possibly accounting for habitat mask.
#'
#' @param y either a matrix or array of encounter history data, possiblity from \code{sim_classic()}.
#' @param M an integer of the total augmented population size (i.e., detected and augmented individuals)
#' @param X either a matrix or array representing the coordinates of traps in
#' UTMs. An array is used when traps are clustered over a survey area.
#' @param buff the distance (m or km) that the traps should be
#' buffered by as an integer. This is typically 3 times the sigma parameter.
#' @param hab_mask either \code{FALSE} (the default) or a matrix or arrary output from \code{\link{mask_polygon}}
#' or \code{\link{mask_raster}} functions.
#' @return a matrix of initial activity center coordinates with \code{M} rows and 2 columns.
#' @details This function generates initial activity center locations based on encounter histories, augmented population size, state-space buffer, and potentially a habitat mask. Note that mean trap detection locations are used for detected individuals while intial values are randomly drawn for augemented individuals. Also, a habitat check will be conducted for all locations when a habitat mask is included.
#' @author Daniel Eacker
#' @examples
#' # simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' traps <- traps + runif(prod(dim(traps)),-20,20) # add some random noise to locations
#'
#' mysigma = 300 # simulate sigma of 300 m
#' mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#'
#' # Create grid and extent
#' Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = 100)
#'
#' # simulate SCR data
#' data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, sigma_ = mysigma, prop_sex = 1, N = 200, K = 4, base_encounter = 0.25, enc_dist = "binomial", hab_mask = FALSE, setSeed = 50)
#'
#' # generate initial activity center coordinates for 2D trap array without habitat mask
#' s.st3d = initialize_classic(y=data3d$y, M=500, X=traps, buff = 3*mysigma, hab_mask = FALSE)
#'
#' # make simple plot
#' par(mfrow=c(1,1))
#' plot(Grid$grid, pch=20,ylab="Northing",xlab="Easting")
#' points(traps, col="blue",pch=20)
#' points(s.st3d, col="red",pch=20)
#' @export
initialize_classic <- function(y, M, X, buff, hab_mask = FALSE){
  if(length(dim(X))!=2 & length(dim(X))!=3){
    stop("Trapping grid must be only 2 or 3 dimensions")
  }
  # for dim of length 2
  if(length(dim(X))==2){
    n0 <- length(which(apply(y,1,sum)!=0))
    s.st <- matrix(NA, nrow=M, ncol=2)
    xlim <- c(min(X[,1] - max(buff)), max(X[,1] + max(buff))) # create x limits for state-space
    ylim <- c(min(X[,2] - max(buff)), max(X[,2] + max(buff))) # create y limits for state-space
    for(i in 1:n0){
      temp.y <- y[i,,]
      temp.X <- X[which(apply(temp.y,1,sum)!=0),]
      ntimes <- apply(temp.y,1,sum)[which(apply(temp.y,1,sum)!=0)]
      df <- data.frame(temp.X,ntimes)
      temp.X2 <- df[rep(seq_len(nrow(df)), df$ntimes),1:2]
      if(isFALSE(hab_mask)){ # if no habitat mask used
        if(length(ntimes)==1){
          s.st[i,1]<-temp.X[1]
          s.st[i,2]<-temp.X[2]
        }else
          if(length(ntimes)>1){
            s.st[i,1]<-mean(temp.X2[,1])
            s.st[i,2]<-mean(temp.X2[,2])
          }
      }else
        if(isFALSE(hab_mask)==FALSE){ # if habitat mask used
          if(length(ntimes)==1){
            s.st[i,1]<-temp.X[1]
            s.st[i,2]<-temp.X[2]
          }else
            if(length(ntimes)>1){
              s.st[i,1]<-mean(temp.X2[,1])
              s.st[i,2]<-mean(temp.X2[,2])
            }
          sx.rescale <- scales::rescale(s.st[i,1], to = c(0,dim(hab_mask)[2]), from=xlim)
          sy.rescale <- scales::rescale(s.st[i,2], to = c(0,dim(hab_mask)[1]), from=ylim)
          pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1)]
          while(pOK==0){ # if their mean activity center is not in suitable habitat, then randomly sample
            s.st[i,1]<-runif(1,xlim[1],xlim[2])
            s.st[i,2]<-runif(1,ylim[1],ylim[2])
            sx.rescale <- scales::rescale(s.st[i,1], to = c(0,dim(hab_mask)[2]), from=xlim)
            sy.rescale <- scales::rescale(s.st[i,2], to = c(0,dim(hab_mask)[1]), from=ylim)
            pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1)]
          }
        }# end habitat check
    } # end detected individuals
    for(i in (n0+1):M){
      if(isFALSE(hab_mask)){ # check for habitat mask for augmented individuals
        s.st[i,1]<-runif(1,xlim[1],xlim[2])
        s.st[i,2]<-runif(1,ylim[1],ylim[2])
      } else
        if(isFALSE(hab_mask)==FALSE){
          s.st[i,1]<-runif(1,xlim[1],xlim[2])
          s.st[i,2]<-runif(1,ylim[1],ylim[2])
          sx.rescale <- scales::rescale(s.st[i,1], to = c(0,dim(hab_mask)[2]), from=xlim)
          sy.rescale <- scales::rescale(s.st[i,2], to = c(0,dim(hab_mask)[1]), from=ylim)
          pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1)]
          while(pOK==0){
            s.st[i,1]<-runif(1,xlim[1],xlim[2])
            s.st[i,2]<-runif(1,ylim[1],ylim[2])
            sx.rescale <- scales::rescale(s.st[i,1], to = c(0,dim(hab_mask)[2]), from=xlim)
            sy.rescale <- scales::rescale(s.st[i,2], to = c(0,dim(hab_mask)[1]), from=ylim)
            pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1)]
          }
        } # end habitat check
    } # augmented individuals
  }  else
    if(length(dim(X))==3){
      n0 <- apply(y, 4, function(x) length(which(apply(x,1,sum)!=0)))
      n0.all <- sum(n0)
      s.st <- matrix(NA, nrow=M*dim(X)[3], ncol=2)
      s.st.array <- array(NA, dim=c(M,2,dim(X)[3]))
      for(g in 1:dim(X)[3]){
        xlim <- c(min(X[,1,g] - max(buff)), max(X[,1,g] + max(buff))) # create x limits for state-space
        ylim <- c(min(X[,2,g] - max(buff)), max(X[,2,g] + max(buff))) # create y limits for state-space
        for(i in 1:n0[g]){
          temp.y <- y[i,,,g]
          temp.X <- X[which(apply(temp.y,1,sum)!=0),,g]
          ntimes <- apply(temp.y,1,sum)[which(apply(temp.y,1,sum)!=0)]
          df <- data.frame(temp.X,ntimes)
          temp.X2 = df[rep(seq_len(nrow(df)), df$ntimes),1:2]
          if(isFALSE(hab_mask)){ # if no habitat mask used
            if(length(ntimes)==1){
              s.st.array[i,1,g]<-temp.X[1]
              s.st.array[i,2,g]<-temp.X[2]
            }else
              if(length(ntimes)>1){
                s.st.array[i,1,g]<-mean(temp.X2[,1])
                s.st.array[i,2,g]<-mean(temp.X2[,2])
              }
          }else
            if(isFALSE(hab_mask)==FALSE){ # if habitat mask used
              if(length(ntimes)==1){
                s.st.array[i,1,g]<-temp.X[1]
                s.st.array[i,2,g]<-temp.X[2]
              }else
                if(length(ntimes)>1){
                  s.st.array[i,1,g]<-mean(temp.X2[,1])
                  s.st.array[i,2,g]<-mean(temp.X2[,2])
                }
              sx.rescale <- scales::rescale(s.st.array[i,1,g], to = c(0,dim(hab_mask)[2]), from=xlim)
              sy.rescale <- scales::rescale(s.st.array[i,2,g], to = c(0,dim(hab_mask)[1]), from=ylim)
              pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1),g]
              while(pOK==0){ # if their mean acitivity center is not in suitable habitat, then randomly sample
                s.st.array[i,1,g]<-runif(1,xlim[1],xlim[2])
                s.st.array[i,2,g]<-runif(1,ylim[1],ylim[2])
                sx.rescale <- scales::rescale(s.st.array[i,1,g], to = c(0,dim(hab_mask)[2]), from=xlim)
                sy.rescale <- scales::rescale(s.st.array[i,2,g], to = c(0,dim(hab_mask)[1]), from=ylim)
                pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1),g]
              }
            } # end habitat check
        } # end detected individuals
      } # close g loop
      s.st[1:n0.all,1:2] <-  na.omit(apply(s.st.array, 2, rbind))
      s.st.array <- array(NA, dim=c(M,2,dim(X)[3]))
      for(g in 1:dim(X)[3]){
        xlim <- c(min(X[,1,g] - max(buff)), max(X[,1,g] + max(buff))) # create x limits for state-space
        ylim <- c(min(X[,2,g] - max(buff)), max(X[,2,g] + max(buff))) # create y limits for state-space
        for(i in (n0[g]+1):M){
          if(isFALSE(hab_mask)){ # check for habitat mask for augmented individuals
            s.st.array[i,1,g]<-runif(1,xlim[1],xlim[2])
            s.st.array[i,2,g]<-runif(1,ylim[1],ylim[2])
          } else
            if(isFALSE(hab_mask)==FALSE){
              s.st.array[i,1,g]<-runif(1,xlim[1],xlim[2])
              s.st.array[i,2,g]<-runif(1,ylim[1],ylim[2])
              sx.rescale <- scales::rescale(s.st.array[i,1,g], to = c(0,dim(hab_mask)[2]), from=xlim)
              sy.rescale <- scales::rescale(s.st.array[i,2,g], to = c(0,dim(hab_mask)[1]), from=ylim)
              pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1),g]
              while(pOK==0){
                s.st.array[i,1,g]=runif(1,xlim[1],xlim[2])
                s.st.array[i,2,g]=runif(1,ylim[1],ylim[2])
                sx.rescale <- scales::rescale(s.st.array[i,1,g], to = c(0,dim(hab_mask)[2]), from=xlim)
                sy.rescale <- scales::rescale(s.st.array[i,2,g], to = c(0,dim(hab_mask)[1]), from=ylim)
                pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1),g]
              }
            } # end habitat check
        } # end augmented individuals
      } # end trap array loop g
      s.st[(n0.all+1):(M*dim(X)[3]),1:2] <-  na.omit(apply(s.st.array, 2, rbind))
    } # end 3D initialize
  return(s.st)
} # end function 'initialize_classic'



#' Function to rescale trap coordinates, grid extent, and starting activity center coordinates
#'
#' Rescale inputs to prepare data for habitat mask to be used.
#'
#' @param X either a matrix or array representing the coordinates of traps in
#' UTMs. An array is used when traps are clustered over a survey area.
#' @param ext an \code{Extent} object from the \code{raster} package. This is returned from \code{\link{grid_classic}}.
#' @param s.st a matrix of starting activity center coordinates. This is returned from \code{\link{initialize_classic}}
#' buffered by as an integer. This is typically 3 times the sigma parameter.
#' @param site Either \code{NULL} (if a 2D trap array is used) or a vector of integers denoting which trap array an individual (either detected or augmented) belongs to. Note that \code{site} is provided from \code{\link{sim_classic}} when a 3D trap array is used. However, this \code{site} variable must be correctly augmented based on the total augmented population size (i.e., \code{M}).
#' @param hab_mask a matrix or arrary output from \code{\link{mask_polygon}} or \code{\link{mask_raster}} functions.
#' @return a list of rescaled trap coordinates, grid extents, and starting activtiy center coordinates.
#' @details This function is only meant to be used when habitat masking is incorporated into the model. The functions properly rescales inputs based on the dimensions of the habitat mask. Note that the \code{pixelWidth} needs to be included as an input in the model after inputs are rescaled to correctly estimate the scaling parameter (i.e., 'sigma').
#' @author Daniel Eacker
#' @importFrom sf st_polygon
#' @seealso \code{\link{mask_polygon}}, \code{\link{mask_raster}}
#' @examples
#'# simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' traps <- traps + runif(prod(dim(traps)),-20,20) # add some random noise to locations
#'
#' mysigma = 300 # simulate sigma of 300 m
#' mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#'
#' # create grid
#' Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = 100)
#'
#' # create polygon to use as a mask
#' library(sf)
#' poly = st_sfc(st_polygon(x=list(matrix(c(-1765,-1765,1730,-1650,1600,1650,0,1350,-800,1700,-1850,1000,-1765,-1765),ncol=2, byrow=TRUE))), crs =  mycrs)
#'
#' # make simple plot
#' par(mfrow=c(1,1))
#' plot(Grid$grid, pch=20)
#' points(traps, col="blue",pch=20)
#' plot(poly, add=TRUE)
#'
#' # create habitat mask
#' hab_mask = mask_polygon(poly = poly, grid = Grid$grid, crs_ = mycrs, prev_mask = NULL)
#'
#' # simulate SCR data
#' data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, sigma_ = mysigma, prop_sex = 1, N = 200, K = 4, base_encounter = 0.25, enc_dist = "binomial", hab_mask = hab_mask, setSeed = 50)
#'
#' # generate initial activity center coordinates for 2D trap array without habitat mask
#' s.st3d = initialize_classic(y=data3d$y, M=500, X=traps, buff = 3*mysigma, hab_mask = hab_mask)
#'
#' # build rescaled constants list for 2D trap array
#' constList = rescale_classic(X = traps, ext = Grid$ext, s.st = s.st3d, site = NULL, hab_mask = hab_mask)
#' str(constList)
#' @export
rescale_classic <- function(X, ext, s.st, site, hab_mask){
  if(length(dim(X))!=2 & length(dim(X))!=3){
    stop("Trapping grid must be only 2 or 3 dimensions")
  }
  if(is.null(hab_mask)){
    stop("Must include habitat mask for rescaling of inputs")
  }
  rescale_list <- list(X=X,ext=ext,s.st=s.st)
  # for dim of length 2
  if(length(dim(X))==2){
    rescale_list$X[,1] <- scales::rescale(X[,1], to = c(0,dim(hab_mask)[2]), from=ext[1:2])
    rescale_list$X[,2] <- scales::rescale(X[,2], to = c(0,dim(hab_mask)[1]), from=ext[3:4])
    rescale_list$s.st[,1] <- scales::rescale(s.st[,1], to = c(0,dim(hab_mask)[2]), from=ext[1:2])
    rescale_list$s.st[,2] <- scales::rescale(s.st[,2], to = c(0,dim(hab_mask)[1]), from=ext[3:4])
    rescale_list$ext[1:2] <- scales::rescale(ext[1:2], to = c(0,dim(hab_mask)[2]), from=ext[1:2])
    rescale_list$ext[3:4] <- scales::rescale(ext[3:4], to = c(0,dim(hab_mask)[1]), from=ext[3:4])

  }else
    if(length(dim(X))==3){
      for(g in 1:dim(X)[3]){ # loop over trap arrays
        rescale_list$X[,1,g] <- scales::rescale(X[,1,g], to = c(0,dim(hab_mask)[2]), from=ext[[g]][1:2])
        rescale_list$X[,2,g] <- scales::rescale(X[,2,g], to = c(0,dim(hab_mask)[1]), from=ext[[g]][3:4])
        rescale_list$s.st[which(site==g),1] <- scales::rescale(s.st[which(site==g),1], to = c(0,dim(hab_mask)[2]), from=ext[[g]][1:2])
        rescale_list$s.st[which(site==g),2] <- scales::rescale(s.st[which(site==g),2], to = c(0,dim(hab_mask)[1]), from=ext[[g]][3:4])
        rescale_list$ext[[g]][1:2] <- scales::rescale(ext[[g]][1:2], to = c(0,dim(hab_mask)[2]), from=ext[[g]][1:2])
        rescale_list$ext[[g]][3:4] <- scales::rescale(ext[[g]][3:4], to = c(0,dim(hab_mask)[1]), from=ext[[g]][3:4])
      }
    }
  return(rescale_list)
} # end function 'rescale_classic'



#' Function to create habitat mask matrix or array from polygon
#'
#' Creates a matrix or array to use as a habitat mask to account for unsuitable habitat
#'
#' @param poly a polygon created using the \code{sf} package of class \code{"sfc_POLYGON"}
#' @param grid a matrix or array object of the the state-space grid. This is returned from \code{\link{grid_classic}}.
#' @param crs_ the UTM coordinate reference system (EPSG code) used for your
#' location provided as an integer (e.g., 32608 for WGS 84/UTM Zone 8N).
#' @param prev_mask either \code{NULL} or a previously created habitat mask matrix or array from \code{\link{mask_polygon}} or \code{\link{mask_raster}}. This allows for habitat masks to be combined to account for different spatial features.
#' @return a matrix or array of 0's and 1's denoting unsuitable and suitable habitat respectively.
#' @details This function creates a habitat matrix or array depending upon whether a 2D (former) or 3D (latter) trap array is used. This matrix can be directly included as data in Bayesian SCR models run using \code{nimble}.
#' @author Daniel Eacker
#' @seealso \code{\link{mask_raster}}
#' @examples
#' # simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' traps <- traps + runif(prod(dim(traps)),-20,20) # add some random noise to locations
#'
#' mysigma = 300 # simulate sigma of 300 m
#' mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#'
#' # create state-space grid and extent
#' Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = 100)
#'
#' # create polygon to use as a mask
#' library(sf)
#' poly = st_sfc(st_polygon(x=list(matrix(c(-1765,-1765,1730,-1650,1600,1650,0,1350,-800,1700,-1850,1000,-1765,-1765),ncol=2, byrow=TRUE))), crs =  mycrs)
#'
#' # make simple plot
#' par(mfrow=c(1,2))
#' plot(Grid$grid, pch=20)
#' points(traps, col="blue",pch=20)
#' plot(poly, add=TRUE)
#'
#' # create habitat mask from polygon
#' hab_mask = mask_polygon(poly = poly, grid = Grid$grid, crs_ = mycrs, prev_mask = NULL)
#'
#' # make simple plot
#' library(raster)
#' plot(raster(apply(hab_mask,2,rev)))
#'
#' # make simple plot
#' poly2 = st_sfc(st_polygon(x=list(matrix(c(-1365,-1365,1730,-1650,1500,1550,0,1350,-800,1700,-1850,1000,-1365,-1365),ncol=2, byrow=TRUE))), crs =  mycrs)
#' plot(poly2, add=TRUE)
#'
#' # mask second polygon, building on previous habitat mask
#' hab_mask2 = mask_polygon(poly = poly2, grid = Grid$grid, crs_ = mycrs, prev_mask = hab_mask)
#'
#' # make simple plot
#' plot(Grid$grid, pch=20)
#' points(traps, col="blue",pch=20)
#' plot(poly, add=TRUE)
#' plot(poly2, add=TRUE)
#' plot(raster(apply(hab_mask2,2,rev)))
#'
#' # create an array of traps, as an approach where individuals will only be detected at one of the trap arrays (e.g., Furnas et al. 2018)
#' Xarray = array(NA, dim=c(nrow(traps),2,2))
#' Xarray[,,1]=traps
#' Xarray[,,2]=traps+4000 # shift trapping grid to new locations
#'
#' # Example of using habitat mask with 3D trap array (need polygon that masks both trapping extents)
#' GridX = grid_classic(X = Xarray, crs_ = mycrs, buff = 3*mysigma, res = 100)
#'
#' # make simple plot
#' par(mfrow=c(1,1))
#' plot(GridX$grid[,,1],xlim=c(-1600,6000),ylim=c(-1600,6000),col="darkgrey",pch=20,ylab="Northing",xlab="Easting")
#' points(Xarray[,,1],col="blue",pch=20)
#' points(GridX$grid[,,2],pch=20,col="darkgrey")
#' points(Xarray[,,2],col="blue",pch=20)
#'
#' poly = st_sfc(st_polygon(x=list(matrix(c(-1660,-1900,5730,-1050,5470,5150,0,6050,-1800,5700,-1660,-1900),ncol=2, byrow=TRUE))), crs =  mycrs)
#' plot(poly, add=TRUE)
#'
#' # get 3D habitat mask array for 3D grid
#' hab_mask = mask_polygon(poly = poly, grid = GridX$grid, crs_ = mycrs, prev_mask = NULL)
#'
#' # make simple plot
#' par(mfrow=c(1,2))
#' apply(hab_mask,3,function(x) plot(raster(apply(x,2,rev))))
#' @export
mask_polygon <- function(poly, grid, crs_, prev_mask){
  # need to determine if X is 2 or 3 dimensional (stop if not)
  if(length(dim(grid))!=2 & length(dim(grid))!=3){
    stop("Trapping grid must be only 2 or 3 dimensions")
  }
  # for dim of length 2
  if(length(dim(grid))==2){
    grid_pts <- sf::st_cast(sf::st_sfc(sf::st_multipoint(grid), crs =  crs_),"POINT")
    habitat_mask <- apply(matrix(ifelse(is.na(as.numeric(sf::st_intersects(grid_pts,poly))),0,
                                        as.numeric(sf::st_intersects(grid_pts,poly))),nrow=dim(raster::rasterFromXYZ(grid,crs=sf::st_crs(crs_)$crs))[1],
                                 ncol=dim(raster::rasterFromXYZ(grid,crs=sf::st_crs(crs_)$crs))[2], byrow=TRUE),2,rev)
    # to combine with a previous habitat mask
    if(is.null(prev_mask)==FALSE){
      # Check to see if dimensions of previous and current habitat
      if(all(dim(habitat_mask) == dim(prev_mask))==FALSE){
        stop("Dimension of previous habitat matrix does not equal dimension of current one")
      }
      habitat_mask <- habitat_mask * prev_mask
    }
  }else
    # for dim of length 3
    if(length(dim(grid))==3){
      grid_pts <- apply(grid, 3, function(x) sf::st_cast(sf::st_sfc(sf::st_multipoint(x), crs =  crs_),"POINT"))
      habitat_mask <- lapply(grid_pts, function(x) apply(matrix(ifelse(is.na(as.numeric(sf::st_intersects(x,poly))),0,
                                                                       as.numeric(sf::st_intersects(x,poly))),nrow=dim(raster::rasterFromXYZ(grid[,,1],crs=sf::st_crs(crs_)$crs))[1],
                                                                ncol=dim(raster::rasterFromXYZ(grid[,,1],crs=sf::st_crs(crs_)$crs))[2], byrow=TRUE),2,rev))
      habitat_mask <- array(unlist(habitat_mask),dim=c(dim(habitat_mask[[1]]),dim(grid)[3])) # convert to array
      # to combine with a previous habitat mask
      if(is.null(prev_mask)==FALSE){
        # Check to see if dimensions of previous and current habitat
        if(all(dim(habitat_mask) == dim(prev_mask))==FALSE){
          stop("Dimension of previous habitat matrix does not equal dimension of current one")
        }
        habitat_mask <- habitat_mask * prev_mask
      }
    }
  return(habitat_mask)
} # End function 'mask_polygon'





#' Function to create habitat mask matrix or array from raster
#'
#' Creates a matrix or array to use as a habitat mask to account for unsuitable habitat
#'
#' @param rast a raster layer created using the \code{raster} package of class \code{"RasterLayer"}
#' @param FUN a function that defines the criteria for suitable habitat.
#' @param grid a matrix or array object of the the state-space grid. This is returned from \code{\link{grid_classic}}.
#' @param crs_ the UTM coordinate reference system (EPSG code) used for your
#' location provided as an integer (e.g., 32608 for WGS 84/UTM Zone 8N).
#' @param prev_mask either \code{NULL} or a previously created habitat mask matrix or array from \code{\link{mask_polygon}} or \code{\link{mask_raster}}. This allows for habitat masks to be combined to account for different spatial features.
#' @return a matrix or array of 0's and 1's denoting unsuitable and suitable habitat respectively.
#' @details This function creates a habitat matrix or array depending upon whether a 2D (former) or 3D (latter) trap array is used. This matrix can be directly included as data in Bayesian SCR models run using \code{nimble}.
#' @author Daniel Eacker
#' @importFrom sp CRS
#' @importFrom sf st_intersects
#' @importFrom raster proj4string extract
#' @importFrom methods as
#' @seealso \code{\link{mask_polygon}}
#' @examples
#' # simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' traps <- traps + runif(prod(dim(traps)),-20,20) # add some random noise to locations
#'
#' mysigma = 300 # simulate sigma of 300 m
#' mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#'
#' # create state-space grid and extent
#' Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = 100)
#'
#' # run previous code used for mask_polygon() to create raster for example
#' library(sf)
#' poly = st_sfc(st_polygon(x=list(matrix(c(-1665,-1665,1730,-1650,1600,1650,0,1350,-800,1700,-1850,1000,-1665,-1665),ncol=2, byrow=TRUE))), crs =  mycrs)
#' hab_mask = mask_polygon(poly = poly, grid = Grid$grid, crs_ = mycrs, prev_mask = NULL)
#'
#' # create raster for demonstration purposes
#' library(raster)
#' rast <- raster(nrow=dim(hab_mask)[1], ncol=dim(hab_mask)[2],ext=Grid$ext,crs=CRS(st_crs(mycrs)$proj4string))
#' rast[] = apply(hab_mask,2,rev)
#'
#' # create habitat mask using raster
#' hab_mask_r = mask_raster(rast = rast, FUN = function(x){x==1}, grid = Grid$grid, crs_ = mycrs, prev_mask = NULL)
#'
#' # make simple plot
#' plot(raster(apply(hab_mask_r,2,rev))) # returns idential results as input rast (but this was just an example raster)
#'
#' # create an array of traps, as an approach where individuals will only be detected at one of the trap arrays (e.g., Furnas et al. 2018)
#' Xarray = array(NA, dim=c(nrow(traps),2,2))
#' Xarray[,,1]=traps
#' Xarray[,,2]=traps+4000 # shift trapping grid to new locations
#'
#' # create grid and extent for 3D trap array
#' GridX = grid_classic(X = Xarray, crs_ = mycrs, buff = 3*mysigma, res = 100)
#'
#' # make simple plot
#' par(mfrow=c(1,1))
#' plot(GridX$grid[,,1],xlim=c(-1600,6000),ylim=c(-1600,6000),col="darkgrey",pch=20,ylab="Northing",xlab="Easting")
#' points(Xarray[,,1],col="blue",pch=20)
#' points(GridX$grid[,,2],pch=20,col="darkgrey")
#' points(Xarray[,,2],col="blue",pch=20)
#'
#' # create polygon to use as a mask and covert to raster
#' poly = st_sfc(st_polygon(x=list(matrix(c(-1660,-1900,5730,-1050,5470,5650,0,6050,-1800,5700,-1660,-1900),ncol=2, byrow=TRUE))), crs =  mycrs)
#'
#' # add polygon to plot
#' plot(poly, add=TRUE)
#'
#' # make raster from polygon
#' rast = raster(xmn=-2000, xmx=6000, ymn=-2000, ymx=6500,res=100,crs=CRS(st_crs(mycrs)$proj4string))
#' rast[]=st_intersects(st_cast(st_sfc(st_multipoint(coordinates(rast)), crs =  mycrs),"POINT"),poly,sparse=FALSE)
#'
#' # make simple plot of raster
#' plot(rast)
#'
#' # get 3D habitat mask array for 3D grid
#' hab_mask = mask_raster(rast = rast, FUN = function(x){x==1}, grid = GridX$grid, crs_ = mycrs, prev_mask = NULL)
#  #make plot
#' par(mfrow=c(1,2))
#' apply(hab_mask,3,function(x) plot(raster(apply(x,2,rev))))
#' @export
mask_raster <- function(rast, FUN, grid, crs_, prev_mask){
  # need to determine if X is 2 or 3 dimensional (stop if not)
  if(length(dim(grid))!=2 & length(dim(grid))!=3){
    stop("Trapping grid must be only 2 or 3 dimensions")
  }
  if(identical(raster::proj4string(rast),sf::st_crs(crs_)$proj4string)==FALSE){
    stop("crs of raster layer must be the same as crs_")
  }
  # for dim of length 2
  if(length(dim(grid))==2){
    grid_pts <- sf::st_cast(sf::st_sfc(sf::st_multipoint(grid), crs =  crs_),"POINT")
    vals <- raster::extract(rast,methods::as(grid_pts,"Spatial"))
    rast_ind <- FUN(vals)
    habitat_mask <- apply(matrix(as.numeric(rast_ind),nrow=dim(raster::rasterFromXYZ(grid,crs=sf::st_crs(crs_)$crs))[1],
                                 ncol=dim(raster::rasterFromXYZ(grid,crs=sf::st_crs(crs_)$crs))[2], byrow=TRUE),2,rev)
    # to combine with a previous habitat mask
    if(is.null(prev_mask)==FALSE){
      # Check to see if dimensions of previous and current habitat
      if(all(dim(habitat_mask) == dim(prev_mask))==FALSE){
        stop("Dimension of previous habitat matrix does not equal dimension of current one")
      }
      habitat_mask <- habitat_mask * prev_mask
    }
  }else
    # for dim of length 3
    if(length(dim(grid))==3){
      grid_pts <- apply(grid, 3, function(x) sf::st_cast(sf::st_sfc(sf::st_multipoint(x), crs =  crs_),"POINT"))
      vals <- lapply(grid_pts,function(x) raster::extract(rast,methods::as(x,"Spatial")))
      rast_ind <- lapply(vals, function(x) FUN(x))
      habitat_mask <- lapply(rast_ind, function(x) apply(matrix(as.numeric(x),nrow=dim(raster::rasterFromXYZ(grid[,,1],crs=sf::st_crs(crs_)$crs))[1],
                                                                ncol=dim(raster::rasterFromXYZ(grid[,,1],crs=sf::st_crs(crs_)$crs))[2], byrow=TRUE),2,rev))
      habitat_mask <- array(unlist(habitat_mask),dim=c(dim(habitat_mask[[1]]),dim(grid)[3])) # convert to array
      # to combine with a previous habitat mask
      if(is.null(prev_mask)==FALSE){
        # Check to see if dimensions of previous and current habitat
        if(all(dim(habitat_mask) == dim(prev_mask))==FALSE){
          stop("Dimension of previous habitat matrix does not equal dimension of current one")
        }
        habitat_mask <- habitat_mask * prev_mask
      }
    }
  return(habitat_mask)
} # End function 'mask_raster'



#' Function to run spatially-explicit capture-recapture models in nimble using parallel processing
#'
#' A wrapper function to conduct Markov Chain Monte Carlo (MCMC) sampling using nimble
#'
#' @param model \code{nimbleCode} used to define model in \code{nimble} package, possibly generated from \code{\link{get_classic}}.
#' @param data a list of data inputs needed to run SCR models in \code{nimble}.
#' @param constants a list of constants needed to run SCR models in \code{nimble}.
#' @param inits starting values for stochastic parameters to begin MCMC sampling.
#' @param params a vector of character strings that define the parameters to trace in the MCMC simulation.
#' @param niter an integer value of the total number of MCMC iterations to run per chain.
#' @param nburnin an integer value of the number of MCMC iterations to discard as 'burnin'.
#' @param thin an integer value of the amount of thinning of the chain. For example, \code{thin=2} would retain every other MCMC sample.
#' @param nchains an integer value for the number of MCMC chains to run
#' @param parallel a logical value indicating whether MCMC chains shoud be run in parallel processing. Default is \code{FALSE}.
#' @param RNGseed an integer value specifying the random number generating seed using in parallel processing. This ensures that the MCMC samples will be the same during each run using the same data, etc.
#' @return a list of MCMC samples for each parameter traced with length equal to the number of chains run.
#' @details This function provides a wrapper to easily run Bayesian SCR models using \code{nimble}.
#' @importFrom tictoc tic toc
#' @importFrom graphics hist par abline lines
#' @importFrom parallel parLapply makeCluster stopCluster clusterSetRNGStream
#' @import nimble
#' @author Daniel Eacker
#' @examples
#' # simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' traps <- traps + runif(prod(dim(traps)),-20,20) # add some random noise to locations
#'
#' mysigma = 300 # simulate sigma of 300 m
#' mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#' pixelWidth = 100 # store pixelWidth
#'
# create grid and extent
#' Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = pixelWidth) # create slightly larger buffer area for example
#'
#' # create polygon to use as a mask
#' library(sf)
#' poly = st_sfc(st_polygon(x=list(matrix(c(-1765,-1765,1730,-1650,1600,1650,0,1350,-800,1700,-1850,1000,-1765,-1765),ncol=2, byrow=TRUE))), crs =  mycrs)
#'
#' # create habitat mask
#' hab_mask = mask_polygon(poly = poly, grid = Grid$grid, crs_ = mycrs, prev_mask = NULL)
#'
#' # simulate data for uniform state-space and habitat mask
#' data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, sigma_ = mysigma, prop_sex = 1, N = 200, K = 4, base_encounter = 0.15, enc_dist = "binomial",hab_mask = hab_mask, setSeed = 50)
#'
#' # get initial activity center starting values
#' s.st3d = initialize_classic(y=data3d$y, M=500, X=traps, buff = 3*mysigma, hab_mask = hab_mask)
#'
#' # rescale inputs
#' rescale_list = rescale_classic(X = traps, ext = Grid$ext, s.st = s.st3d, hab_mask = hab_mask)
#'
#' # store rescaled extent
#' ext = rescale_list$ext
#'
#' # prepare data
#' data = list(y=data3d$y)
#' data$y = data$y[which(apply(data$y, 1, sum)!=0),,] # remove augmented records (won't matter though)
#' data$y = apply(data$y, c(1,2), sum) # covert to 2d by summing over occasions
#'
#' # add rescaled traps
#' data$X = rescale_list$X
#'
#' # prepare constants
#' constants = list(M = 500,n0 = nrow(data$y),J=dim(data$y)[2], K=dim(data3d$y)[3],
#'                  x_lower = ext[1], x_upper = ext[2], y_lower = ext[3], y_upper = ext[4],
#'                  sigma_upper = 1000, A = sum(hab_mask)*(pixelWidth)^2,pixelWidth=pixelWidth)
#'
#' # add z and zeros vector data for latent inclusion indicator
#' data$z = c(rep(1,constants$n0),rep(NA,constants$M - constants$n0))
#' data$zeros =  c(rep(NA,constants$n0),rep(0,constants$M - constants$n0))
#'
#' # add hab_mask and OK for habitat check
#' data$hab_mask = hab_mask
#' data$OK = rep(1,constants$M)
#'
#' # get initial activity center starting values
#' s.st3d = rescale_list$s.st
#'
#' # define all initial values
#' inits = list(sigma = runif(1, 250, 350), s = s.st3d,psi=runif(1,0.2,0.3), p0 = runif(1, 0.05, 0.15),
#'              pOK = data$OK, z=c(rep(NA,constants$n0),rep(0,constants$M-constants$n0)))
#'
#' # parameters to monitor
#' params = c("sigma","psi","p0","N","D")
#'
#' # get model
#' scr_model = get_classic(dim_y = 2, enc_dist = "binomial",sex_sigma = FALSE,hab_mask=TRUE)
#'
#' # run model
#' library(tictoc)
#' tic() # track time elapsed
#' out = run_classic(model = scr_model, data=data, constants=constants, inits=inits, params = params,
#'                   niter = 1000, nburnin=500, thin=1, nchains=2, parallel=TRUE, RNGseed = 500)
#' toc()
#'
#' # summarize output
#' samples = do.call(rbind, out)
#' par(mfrow=c(1,1))
#' hist(samples[,which(dimnames(out[[1]])[[2]]=="N")], xlab = "Abundance", xlim = c(0,500), main="")
#' abline(v=200, col="red") # add line for simulated abundance
#'
#' # not run
#' #nimSummary(out)
#' @export
run_classic <- function(model, data, constants, inits, params,
                        niter = 1000, nburnin=100, thin=1, nchains=1, parallel=FALSE, RNGseed){
  # for parallel processing
  if(parallel == FALSE){
    SCRmodelR <- nimble::nimbleModel(code=model,data=data,constants=constants,inits=inits,check=FALSE,calculate=TRUE)
    SCRmodelR$initializeInfo()
    # compile model to C++#
    SCRmodelC <- nimble::compileNimble(SCRmodelR) # compile code
    mcmcspec<-nimble::configureMCMC(SCRmodelR, monitors=params) # MCMC configurations
    # block updating for marked individual activity centers
    mcmcspec$removeSamplers("s", print = FALSE)
    for(i in 1:constants$M){
      snew = paste("s[",i,","," 1:2","]",sep="")
      mcmcspec$addSampler(target = snew, type = 'RW_block', silent = TRUE)
    }
    scrMCMC <- nimble::buildMCMC(mcmcspec)
    CscrMCMC <- nimble::compileNimble(scrMCMC, project = SCRmodelR,resetFunctions = TRUE)
    results <- nimble::runMCMC(CscrMCMC, niter = niter, nburnin=nburnin, thin=thin, nchains=nchains)
    return(results)
  }else
    if(parallel == TRUE){
      run_parallel <- function(seed,data,model,constants,inits, params,
                               niter,nburnin,thin){
        SCRmodelR <- nimble::nimbleModel(code=model,data=data,constants=constants,inits=inits,check=FALSE,calculate=TRUE)
        SCRmodelR$initializeInfo()
        # compile model to C++#
        SCRmodelC <- nimble::compileNimble(SCRmodelR) # compile code
        mcmcspec<-nimble::configureMCMC(SCRmodelR, monitors=params) # MCMC configurations
        # block updating for marked individual activity centers
        mcmcspec$removeSamplers("s", print = FALSE)
        for(i in 1:constants$M){
          snew = paste("s[",i,","," 1:2","]",sep="")
          mcmcspec$addSampler(target = snew, type = 'RW_block')
        }
        scrMCMC <- nimble::buildMCMC(mcmcspec)
        CscrMCMC <- nimble::compileNimble(scrMCMC, project = SCRmodelR,resetFunctions = TRUE)
        results <- nimble::runMCMC(CscrMCMC, niter = niter, nburnin= nburnin,thin=thin,setSeed = seed)
        return(results)
      } # end parallel processing function
      if(nchains < 2){
        stop("Must have at least 2 chains to run parallel processing")
      }
      this_cluster <- parallel::makeCluster(nchains)
      parallel::clusterEvalQ(this_cluster, library("nimble")) 
      if(is.null(RNGseed)==FALSE){
        parallel::clusterSetRNGStream(cl = this_cluster, RNGseed)
      }
      chain_output <- parallel::parLapply(cl = this_cluster, fun = run_parallel,X=1:nchains,
                                          model = model, data=data,constants=constants, inits=inits, params = params,
                                          niter = niter, nburnin=nburnin, thin=thin)
      parallel::stopCluster(this_cluster)
      return(chain_output)
    }
} # End function 'run_classic'




#' Function to summarize MCMC chain output from nimble
#'
#' Summarizes lists of MCMC chain output from nimble
#'
#' @param d a list of MCMC samples from each chain returned from \code{\link{run_classic}}.
#' @param trace a logical value indicating whether or not traces of MCMC samples should be plotted. Default is \code{FALSE}.
#' @param plot_all a logical value indicating whether or not all parameters should be included in plots. This assumes that some parameters
#' are excluded in the summary table (i.e., \code{exclude.params != NULL}). Default is \code{FALSE}.
#' @param exclude.params either \code{NULL} or a scalar or vector containing parameter(s) to exclude from summary. Note that high dimensional parameters (e.g., \code{s[1, 1, 1]}) can be excluded by calling \code{exclude.params = "s"}. Default is \code{NULL}.
#' @param digits an integer value indicating how many digits the output should be rounded to.
#' @return a dataframe of summary statistics for MCMC samples.
#' @details This function summarizes Bayesian SCR models run using \code{nimble} including mean and quantiles of samples, as well as effective sample size and Rhat statistics. Note that \code{f0} is the proportion of samples that are greater than zero. Also, at least 2 chains must be run to use this function.
#' @author Daniel Eacker
#' @importFrom stringr str_extract str_locate
#' @importFrom coda effectiveSize mcmc.list as.mcmc gelman.diag
#' @importFrom stats na.omit sd quantile
#' @importFrom crayon red
#' @examples
#' # simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' traps <- traps + runif(prod(dim(traps)),-20,20) # add some random noise to locations
#'
#' mysigma = 300 # simulate sigma of 300 m
#' mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#' pixelWidth = 100 # store pixelWidth
#'
#' # create grid and extent
#' Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = pixelWidth)
#'
#' # simulate encounter data
#' data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, sigma = mysigma, prop_sex = 1, N = 200, K = 4, base_encounter = 0.15, enc_dist = "binomial", hab_mask = FALSE, setSeed = 200)
#'
#' # prepare data
#' data = list(y=data3d$y)
#' data$y = data$y[which(apply(data$y, 1, sum)!=0),,] # remove augmented records
#' data$y = apply(data$y, c(1,2), sum) # covert to 2d by summing over occasions
#' data$X = traps/1000 # rescale to kilometers
#' ext = as.vector(Grid$ext)/1000 # recale to kilometers
#'
#' # prepare constants
#' constants = list(M = 500,n0 = nrow(data$y),J=dim(data$y)[2], K=dim(data3d$y)[3],
#'                  x_lower = ext[1], x_upper = ext[2], y_lower = ext[3], y_upper = ext[4],
#'                  sigma_upper = 1, A = prod(ext[2]-ext[1],ext[4]-ext[3]))
#'
#' # add z and zeros vector data for latent inclusion indicator
#' data$z = c(rep(1,constants$n0),rep(NA,constants$M - constants$n0))
#' data$zeros =  c(rep(NA,constants$n0),rep(0,constants$M - constants$n0))
#'
#' # get initial activity center starting values
#' s.st3d = initialize_classic(y=data3d$y, M=500, X=traps, buff = 3*mysigma, hab_mask=FALSE)
#'
#' # define all initial values
#' inits = list(sigma = runif(1, 0.250, 0.350), s = s.st3d/1000,psi=runif(1,0.2,0.3), p0 = runif(1, 0.05, 0.15),
#'              z=c(rep(NA,constants$n0),rep(0,constants$M-constants$n0)))
#'
#' # parameters to monitor
#' params = c("sigma","psi","p0","N","D")
#'
#' # get model
#' scr_model = get_classic(dim_y = 2, enc_dist = "binomial",sex_sigma = FALSE)
#'
#' # run model
#' library(tictoc)
#' tic() # track time elapsed
#' out = run_classic(model = scr_model, data=data, constants=constants, inits=inits, params = params,
#'                   niter = 5000, nburnin=1000, thin=1, nchains=2, parallel=TRUE, RNGseed = 500)
#' toc()
#'
#' # summarize output
#' samples = do.call(rbind, out)
#' par(mfrow=c(1,1))
#' hist(samples[,which(dimnames(out[[1]])[[2]]=="N")], xlab = "Abundance", xlim = c(0,500), main="")
#' abline(v=200, col="red") # add line for simulated abundance
#'
#' # summarize output
#' nimSummary(out, trace=TRUE, plot_all=TRUE)
#' @export
nimSummary <- function(d, trace=FALSE, plot_all=FALSE, exclude.params = NULL, digits=3){
  if(is.null(exclude.params)==FALSE){
    tmp1 <- ifelse(is.na(as.numeric(stringr::str_extract(attributes(d[[1]])$dimnames[[2]],"[1-9]+"))),attributes(d[[1]])$dimnames[[2]],substr(attributes(d[[1]])$dimnames[[2]],1,as.numeric(stringr::str_locate(attributes(d[[1]])$dimnames[[2]], "\\[")[, 1])-1))
    d.remove <- lapply(d, function(x) which(tmp1 %in% exclude.params))
    d2 <- lapply(d, function(x) x[,-d.remove[[1]]])
  }else
    if(is.null(exclude.params)){
      if(dim(d[[1]])[2]>100){
       param_length_check <- readline(cat(crayon::red(paste("Do you really want to build a summary table for",ncol(d3),"parameters? (1 = Yes or 0 = No); note that this will take forever! "))))
       if(param_length_check=="0"){stop("Exiting function...",call. = FALSE)}
      }
      d2 <- d
      d.remove <- 0
    }
  if((length(attributes(d[[1]])$dimnames[[2]])-length(d.remove[[1]])==1)){
    d3 <- as.data.frame(do.call(c, d2))
    Means <- mean(d3[,1], na.rm=TRUE)
    SDs <- stats::sd(d3[,1], na.rm=TRUE)
    q2.5 <- stats::quantile(d3[,1], 0.025, na.rm=TRUE)
    q50 <- stats::quantile(d3[,1], 0.50, na.rm=TRUE)
    q97.5 <- stats::quantile(d3[,1], 0.975, na.rm=TRUE)
    over.zero <- round(mean(d3[,1]>0),2)
    n.eff <- coda::effectiveSize(mcmc.list(lapply(d2, coda::as.mcmc)))
    Rhat <- round(coda::gelman.diag(coda::mcmc.list(lapply(d2, coda::as.mcmc)), multivariate = FALSE)[[1]][,1],3)
  }else
    if((length(attributes(d[[1]])$dimnames[[2]])-length(d.remove[[1]])>1)){
      d3 <- do.call(rbind,d2)
      Means <- apply(d3, 2,function(x) mean(x,na.rm=TRUE))
      SDs <- apply(d3, 2,function(x) stats::sd(x,na.rm=TRUE))
      q2.5 <- apply(d3, 2,function(x) stats::quantile(x, 0.025,na.rm=TRUE))
      q50 <- apply(d3, 2,function(x) stats::quantile(x, 0.50,na.rm=TRUE))
      q97.5 <- apply(d3, 2,function(x) stats::quantile(x, 0.975,na.rm=TRUE))
      over.zero <- round(apply(d3, 2, function(x) mean(x>0,na.rm=TRUE)),2)
      n.eff <- coda::effectiveSize(coda::mcmc.list(lapply(d2, coda::as.mcmc)))
      Rhat <- round(coda::gelman.diag(coda::mcmc.list(lapply(d2, coda::as.mcmc)), multivariate = FALSE)[[1]][,1],3)
    }
  if(trace==TRUE  & (length(attributes(d[[1]])$dimnames[[2]])-length(d.remove[[1]])>1)){
     par(mfrow=c(3,2))
     g <- 1 # set g index
     if(plot_all){
     d2 <- d
     d3 <- do.call(rbind, d2) # allow all parameters to be plotted despite excluding them from the summary table
     }
     plot.seq <- seq(3,3000,3)  # probably unlikely to want to plot more than 3,000 parameters
     for(i in 1:3){ # plot first 3 variables to start out
      plot(1:dim(d2[[1]])[1],d2[[1]][,i],xlab="iteration",ylab=colnames(d3)[i],type="l",ylim=range(do.call(rbind, lapply(d2,function(x) apply(x, 2, range)))[,i]))
      for(j in 2:length(d2)){
        lines(1:dim(d2[[1]])[1],d2[[j]][,i],xlab="iteration",ylab=colnames(d3)[i],type="l",col="red")
      }
      hist(d3[,i],main="",xlab=colnames(d3)[i])
     }
      if(interactive() & ncol(d3) > 3){
        answer <- readline(cat(crayon::red("Plot next set of parameters? (1 = Yes or 0 = No) ")))
        while(answer == "1"){
          upper_index <- ifelse(plot.seq[g+1] > ncol(d3), ncol(d3), plot.seq[g+1])
          for(i in (plot.seq[g]+1):upper_index){
             plot(1:dim(d2[[1]])[1],d2[[1]][,i],xlab="iteration",ylab=colnames(d3)[i],type="l",ylim=range(do.call(rbind, lapply(d2,function(x) apply(x, 2, range)))[,i]))
              for(j in 2:length(d2)){
              lines(1:dim(d2[[1]])[1],d2[[j]][,i],xlab="iteration",ylab=colnames(d3)[i],type="l",col="red")
             }
             hist(d3[,i],main="",xlab=colnames(d3)[i])
            } # end i loop
          g <- g + 1 
          if(plot.seq[g+1] <= ncol(d3)){
          answer <- readline(cat(crayon::red("Plot next set of parameters? (1 = Yes or 0 = No) ")))
          }else
          if(plot.seq[g+1] > ncol(d3)){
          break
          }
        }
      } # interactive?
  }else
    if(trace==TRUE  & (length(attributes(d[[1]])$dimnames[[2]])-length(d.remove[[1]])==1)){
      par(mfrow=c(1,2))
      plot(1:length(d2[[1]]),d2[[1]],xlab="iteration",ylab=colnames(d3)[i],type="l",ylim=range(d3[,1]))
      for(j in 2:length(d2)){
        lines(1:length(d2[[j]]),d2[[j]],xlab="iteration",ylab=colnames(d3)[i],type="l",col="red")
      }
      hist(d3[,1],main="",xlab=attributes(d[[1]])$dimnames[[2]][-d.remove[[1]]])
    }
  tmp.frame = data.frame(post.mean=Means,post.sd=SDs,q2.5=q2.5,q50=q50,q97.5=q97.5,f0=over.zero,n.eff=n.eff,Rhat=Rhat)
  if(nrow(tmp.frame)==1){
    row.names(tmp.frame) = attributes(d[[1]])$dimnames[[2]][-d.remove[[1]]]
  }
  on.exit(par(mfrow=c(1,1)))
  return(round(tmp.frame, digits=digits))
}



