#' Function to define state-space for spatial capture-recapture models
#'
#' @description Returns a list of a matrix or array object of grid
#' coordinates and an Extent object from the \code{raster} package as a 
#' state-space.
#'
#' @param X Either a matrix or array representing the coordinates of traps in
#' UTMs. An array is used when traps are clustered over a survey area.
#' @param crs_ The UTM coordinate reference system (EPSG code) used for your
#' location provided as an integer (e.g., 32608 for WGS 84/UTM Zone 8N).
#' @param buff The distance (m or km) that the traps should be
#' buffered by as an integer. This is typically 3 times the sigma parameter.
#' @param res The grid cell resolution for the state-space.
#' @return A list of a matrix or array of grid coordinates \code{grid} and an
#' extent object \code{ext}. Note that a matrix object is returned if the 
#' coordinates of the traps are a matrix (i.e., 2D), otherwise an array object
#' is returned when trap coordinates are in a 3D array.
#' @details This function supports spatial capture-recapture analysis by
#' creating two outputs that are used to define the state-space. If a habitat 
#' mask is not used, then only the \code{Extent} object from the \code{raster}
#' package is needed under a uniform state-space.
#' The matrix or array object can be used to develop a habitat mask in a uniform
#' state-space or as a discretized state-space.
#' @author Daniel Eacker
#' @seealso \code{\link[raster:extent]{raster::extent}}
#' @importFrom raster extent coordinates raster rasterFromXYZ ncell proj4string
#' @importFrom sf st_crs 
#' @examples
#' # simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' # add some random noise to locations
#' traps <- traps + runif(prod(dim(traps)),-20,20) 
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
#' @name grid_classic
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
    ext <- raster::extent(c(apply(X, 2, min)-max(buff),apply(X, 2, max)+
           max(buff))[c(1,3,2,4)])
    rast <- raster::raster(ext = ext, crs = crs_, res = res)
    gridOut <- raster::coordinates(rast)
  }else
    # for dim of length 3 (assumes different sites are in dimension 3)
    if(length(dim(X))==3){
      ext_initial <- apply(X, 3, function(x) raster::extent(c(apply(x, 2, min)-
                          max(buff),apply(x, 2, max)+max(buff))[c(1,3,2,4)]))
      centroids <- lapply(ext_initial, function(x) cbind(mean(c(x[2],x[1])),
                          mean(c(x[4],x[3]))))
      x_lengths <- unlist(lapply(ext_initial, function(x) x[2]-x[1]))
      x_index <- which(x_lengths==max(x_lengths))[1]
      y_lengths <- unlist(lapply(ext_initial, function(x) x[4]-x[3]))
      y_index <- which(y_lengths==max(y_lengths))[1]
      ext_standard <-extent(c(-x_lengths[x_index]/2,x_lengths[x_index]/
                         2,-y_lengths[y_index]/2,y_lengths[y_index]/2))
      coords <- raster::coordinates(raster(ext = ext_standard, crs = crs_, 
                         res = res))
      rast <- lapply(centroids, function(x) raster::rasterFromXYZ(
                         cbind(coords[,1]+x[,1],coords[,2]+x[,2]),crs = crs_))
      gridOut <- lapply(rast, raster::coordinates)
      ext <- lapply(rast, function(x) raster::extent(x))
      gridOut <- array(unlist(gridOut),dim=c(raster::ncell(rast[[1]]),
                        dim(X)[2:3])) # convert to array
    }
  return(list(grid = gridOut, ext = ext))
} # End 'grid_classic' function



#' Function to simulate basic spatial capture-recapture data
#'
#' @description Returns a list of simulated data including the encounter
#' history, binary sex indicator, activity centers, and site identifier.
#'
#' @param X Either a matrix or array object representing the coordinates of traps in
#' UTMs. An array is used when traps are clustered over a survey area.
#' @param ext An \code{Extent} object from the \code{raster} package. This is 
#' returned from \code{\link{grid_classic}}.
#' @param crs_ The UTM coordinate reference system (EPSG code) used for your
#' location provided as an integer (e.g., 32608 for WGS 84/UTM Zone 8N).
#' @param N Simulated total abundance as an integer. 
#' @param sigma_ The scaling parameter of the bivariate normal kernel either
#' in meters or kilometers as an integer.
#' @param prop_sex The portion of females or males as a numeric value. This
#' will depend upon the indicator coding scheme used (e.g., females = 1 and
#' males = 0; then proportion of females in the simulation). Must be a numeric
#' value between 0 and 1. Note that 0 or 1 can be used if a non-sex-specific 
#' sigma is desired.
#' @param K The number of sampling occasions desired as an integer.
#' @param base_encounter The baseline encounter probability or rate as a numeric
#' value. Note that a probabilty is given for a \code{"binomial"} observation
#' distribution while a rate is given for a \code{"poisson"} distribution.
#' @param enc_dist Either \code{"binomial"} or \code{"poisson"}. Default is
#' \code{"binomial"}.
#' @param hab_mask Either \code{FALSE} (the default) or a matrix or arrary output from \code{\link{mask_polygon}}
#' or \code{\link{mask_raster}} functions.
#' @param setSeed The random number generater seed as an integer used in
#' simulations to obtain repeatable data simulations. Default is 500.
#' @return A list of a matrix or array of encounter histories \code{y}, a
#' vector or matrix of 0's and 1's for \code{sex}, a batrix of simulated activity
#' centers \code{s}, and when a 3-dimensional trap array is given, a vector
#' for the site identifier \code{site}.
#' @author Daniel Eacker
#' @details This function supports spatial capture-recapture (SCR) analysis by 
#' allowing for easy simulation of data components used by nimble in Baysian 
#' SCR models. Note that the output for the encounter histories \code{y} will be
#' sorted by detected and not detected individuals.
#' @seealso \code{\link{grid_classic}}
#' @importFrom stats runif rbinom rpois
#' @importFrom scales rescale
#' @importFrom sf st_cast st_sfc st_multipoint st_distance
#' @examples
#' # simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' 
#' # add some random noise to locations
#' traps <- traps + runif(prod(dim(traps)),-20,20) 
#'
#' mysigma = 300 # simulate sigma of 300 m
#' mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#'
#' # Create grid and extent
#' Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = 100)
#'
#' # simulate SCR data
#' data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, sigma_ = c(300), 
#' prop_sex = 1, N = 200, K = 4, base_encounter = 0.25, enc_dist = "binomial", 
#' hab_mask = FALSE, setSeed = 50)
#'
#' # make simple plot
#' par(mfrow=c(1,1))
#' plot(Grid$grid, pch=20,ylab="Northing",xlab="Easting")
#' points(traps, col="blue",pch=20)
#' points(data3d$s,col="red",pch = 20)
#' points(data3d$s[which(apply(data3d$y,1,sum)!=0),],col="green",pch = 20)
#' @name sim_classic
#' @export
sim_classic <- function(X, ext, crs_, N, sigma_, prop_sex, K, base_encounter, 
                enc_dist = "binomial", hab_mask = FALSE, setSeed = 500){
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
      sx <- stats::runif(N,xlim[1],xlim[2])
      sy <- stats::runif(N,ylim[1],ylim[2])
      s <- cbind(sx,sy)
    } else
      if(isFALSE(hab_mask)==FALSE){
        sx <- stats::runif(N,xlim[1],xlim[2])
        sx.rescale <- scales::rescale(sx, to = c(0,dim(hab_mask)[2]), from=xlim)
        sy <- stats::runif(N,ylim[1],ylim[2])
        sy.rescale <- scales::rescale(sy, to = c(0,dim(hab_mask)[1]), from=ylim)
        pOK <- numeric(N)
        for(i in 1:N){
          pOK[i] <- hab_mask[(trunc(sy.rescale[i])+1),(trunc(sx.rescale[i])+1)]
          while(pOK[i]==0){
            sx[i] <- stats::runif(1,xlim[1],xlim[2])
            sx.rescale[i] <- scales::rescale(sx[i], to = c(0,dim(hab_mask)[2]), 
                             from=xlim)
            sy[i] <- stats::runif(1,ylim[1],ylim[2])
            sy.rescale[i] <- scales::rescale(sy[i], to = c(0,dim(hab_mask)[1]), 
                              from=ylim)
            pOK[i] <- hab_mask[(trunc(sy.rescale[i])+1),(trunc(sx.rescale[i])+1)]
          }
        }
        s <- cbind(sx,sy)
      }
     # set propsex to 1 if you only want to simulate for one sex
    sex <- stats::rbinom(N, 1, prop_sex)
    sex_ <- sex + 1 # sex indicator for sigma
    if(length(sigma_) == 1){ # make sigma length 2 if only 1 given
      sigma_ <- c(sigma_,sigma_)
    }
     # compute distance matrix for traps
    Dmat <- sf::st_distance(sf::st_cast(sf::st_sfc(sf::st_multipoint(s), 
                                                   crs =  crs_),"POINT"),
                            sf::st_cast(sf::st_sfc(sf::st_multipoint(X), 
                                                   crs =  crs_),"POINT")) 
    # convert manually to matrix
    Dmat <- matrix(as.numeric(Dmat), nrow=nrow(Dmat), ncol=ncol(Dmat)) 
    # get detection prob matrix
    prob <- base_encounter*exp(-Dmat^2/(2*sigma_[sex_]^2)) 
    Y3d <- array(NA, dim=c(N,nrow(X),K)) # encounter data
    if(enc_dist == "binomial"){
      for(i in 1:N){
        for(j in 1:nrow(X)){
          Y3d[i,j,1:K] <- stats::rbinom(K,1,prob[i,j])
        }}
    }else
      if(enc_dist == "poisson"){
        for(i in 1:N){
          for(j in 1:nrow(X)){
            Y3d[i,j,1:K] <- stats::rpois(K,prob[i,j])
          }}
      }
    # organize encountered individuals and then 0's (augmented)
    Y3d <- Y3d[c(which(apply(Y3d,1,sum)!=0),which(apply(Y3d,1,sum)==0)),,] 
    # organize simulated activity centers
    s <- s[c(which(apply(Y3d,1,sum)!=0),which(apply(Y3d,1,sum)==0)),]
    sex[which(apply(Y3d,1,sum)==0)] <- NA
  }else
    # for dim of length 3 (assumes different sites are in dimension 3)
    if(length(dim(X))==3){
      if(length(sigma_) == 1){ # make sigma length 2 if only 1 given
        sigma_ <- c(sigma_,sigma_)
      }
      site <- sort(rep(1:dim(X)[3],N/dim(X)[3]))
      if(length(site) < N){ # check if site variable too short
        N <- length(site)
        warning("N was reduced from requested number of simulated individuals")
      }  
      sex <- stats::rbinom(N, 1, prop_sex) 
      sex_ <- sex + 1 # indicator for sigma
      xlim <-sapply(ext, function(x) x[1:2]) # create x limits for state-space
      ylim <- sapply(ext, function(x) x[3:4]) # create y limits for state-space
        if(isFALSE(hab_mask)){
          sx <- stats::runif(N,xlim[1,site],xlim[2,site])
          sy <- stats::runif(N,ylim[1,site],ylim[2,site])
          s <- cbind(sx,sy)
        } else
          if(isFALSE(hab_mask)==FALSE){
              pOK <- numeric(N)
              s <- matrix(NA, nrow=N, ncol=2)
            for(i in 1:N){
               s[i,1] <- runif(1,xlim[1,site[i]],xlim[2,site[i]])
               sx.rescale <- scales::rescale(s[i,1], to = c(0,dim(hab_mask)[2]),
                                             from=xlim[,site[i]])
               s[i,2] <- runif(1,ylim[1,site[i]],ylim[2,site[i]])
               sy.rescale <- scales::rescale(s[i,2], to = c(0,dim(hab_mask)[1]),
                                             from=ylim[,site[i]])
               pOK[i] <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1),
                                  site[i]]
              while(pOK[i]==0){
               s[i,1] <- runif(1,xlim[1,site[i]],xlim[2,site[i]])
               sx.rescale <- scales::rescale(s[i,1], to = c(0,dim(hab_mask)[2]),
                                             from=xlim[,site[i]])
               s[i,2] <- runif(1,ylim[1,site[i]],ylim[2,site[i]])
               sy.rescale <- scales::rescale(s[i,2], to = c(0,dim(hab_mask)[1]), 
                                             from=ylim[,site[i]])
               pOK[i] <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1),
                                  site[i]]
              }
            }
          }
            Y4d <- array(NA, dim=c(N,dim(X)[1],K)) # encounter data
            prob <- list()
            # estimate probability matrices by site to save run time
          for(g in 1:dim(X)[3]){ 
            Dmat <- sf::st_distance(sf::st_cast(sf::st_sfc(
              sf::st_multipoint(s[site==g,]), crs = crs_),"POINT"),
              sf::st_cast(sf::st_sfc(sf::st_multipoint(X[,,g]), crs =  crs_),
                      "POINT"))  
            # compute distance matrix for traps and  convert manually to matrix     
            Dmat <- matrix(as.numeric(Dmat), nrow=nrow(Dmat), ncol=ncol(Dmat))
             # get detection prob matrix
            prob[[g]] <- base_encounter*exp(-Dmat^2/(2*sigma_[sex_[site==g]]^2))
          }
            prob = do.call(rbind, prob)
      for(i in 1:N){
        if(enc_dist == "binomial"){
            for(j in 1:dim(X)[1]){
              Y4d[i,j,1:K] <- stats::rbinom(K,1,prob[i,j])
            }
        }else
          if(enc_dist == "poisson"){
              for(j in 1:dim(X)[1]){
                Y4d[i,j,1:K] <- stats::rpois(K,prob[i,j])
              }
          }
        }
     } # end 3D loop
  if(length(dim(X))==2){
    dataList <- list(y=Y3d,sex=sex,s=s)
  }else
    if(length(dim(X))==3){
      # organize data elements s, site, sex, and y
      site <- site[c(which(apply(Y4d,1,sum)!=0),which(apply(Y4d,1,sum)==0))]
      sex <- sex[c(which(apply(Y4d,1,sum)!=0),which(apply(Y4d,1,sum)==0))]
      sex[(length(which(apply(Y4d,1,sum)!=0))+1):N] <- NA
      s <- s[c(which(apply(Y4d,1,sum)!=0),which(apply(Y4d,1,sum)==0)),1:2]
      Y4d <- Y4d[c(which(apply(Y4d,1,sum)!=0),which(apply(Y4d,1,sum)==0)),,] 
      dataList <- list(y=Y4d,sex=sex,site=site,s=s)
    }
  return(dataList)
} # End 'sim_encounter' function



#' Function to retrieve nimbleCode for spatial capture-recapture models
#'
#' @description Creates model code using the \code{\link[nimble]{nimbleCode}} function.
#'
#' @param dim_y An integer of either 2 (the default) or that defines what 
#' dimensional format the encounter history data are in.
#' @param enc_dist Either \code{"binomial"} or \code{"poisson"}. Default is
#' \code{"binomial"}.
#' @param sex_sigma A logical value indicating whether the scaling parameter 
#' ('sigma') is sex-specific
#' @param hab_mask A logical value indicating whether a habitat mask will be 
#' used. Default is \code{FALSE}.
#' @param trapsClustered A logical value indicating if traps are clustered in 
#' arrays across the sampling area.
#' @return A \code{nimbleCode} object from the \code{nimble} package.
#' @details This function provides templates that could be copied and easily 
#' modified to include further model complexity such as covariates explaining 
#' detection probability. The models include different encounter probability 
#' distributions, sex-specific scaling parameters, and habitat masking.
#' @author Daniel Eacker
#' @examples
#' # get model for 2D encounter data, binomial encounter distribution, 
#' # non-sex-specific scaling parameter, and no habitat mask
#' scr_model = get_classic(dim_y = 2,enc_dist = "binomial",sex_sigma = FALSE,
#'                          hab_mask = FALSE,trapsClustered = FALSE)
#'
#' # inspect model
#' scr_model
#' @name get_classic
#' @export
get_classic <- function(dim_y, enc_dist = "binomial",sex_sigma = FALSE,
                        hab_mask = FALSE,trapsClustered = FALSE){
    M <- J <- s <- X <- p0 <- sigma <- n0 <- z <- A <- lam0 <- K <- sex <- 
      nSites <-  site <- pixelWidth <- psi <- prop.habitat <- NULL
  if(dim_y !=2 & dim_y!=3){
    stop("dim_y must be either 2 or 3")
  }
  if(enc_dist != "poisson" & enc_dist != "binomial"){
    stop("Encounter distribution has to be either binomial or poisson")
  }
 if(trapsClustered == FALSE){    
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
      }  # End 3D models
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
   }  # End 3D model
  } # end models with hab_mask
 }else # trapsClustered 
 if(trapsClustered){
 if(isFALSE(hab_mask)){ # determine if hab_mask is included
  if(dim_y == 2){
   if(enc_dist == "binomial" & sex_sigma  == FALSE){
      scrcode <- nimble::nimbleCode({
     for(g in 1:nSites){  
      p0[g] ~ dunif(0,1) # baseline encounter probability
    } # g sites
  sigma ~ dunif(0, sigma_upper) # scaling parameter
  psi ~ dunif(0, 1) # inclusion prob
 for(i in 1:M){
  z[i]~dbern(psi)
  s[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
  s[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
  dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1,site[i]])^2 + (s[i,2]-X[1:J,2,site[i]])^2)
  p[i,1:J] <- p0[site[i]]*exp(-dist[i,1:J]^2/(2*sigma^2))
  } # i individuals
  # use zeros trick for individuals to speed up the computation
  for(i in 1:n0){
    for(j in 1:J){
       y[i,j] ~ dbin(p[i,j],K)
    } # j traps
  } # i individuals
  for(i in (n0+1):M){
    zeros[i] ~ dbern((1 - prod(1 - p[i,1:J])^K)*z[i])
  } # i individuals
    N <- sum(z[1:M])
    D <- N/A
    })
 } else
 if(enc_dist == "poisson" & sex_sigma  == FALSE){
   scrcode <- nimble::nimbleCode({
  for(g in 1:nSites){  
    lam0[g] ~ dunif(0,lam0_upper) # baseline encounter rate
  }
    sigma ~ dunif(0, sigma_upper) # scaling parameter
    psi ~ dunif(0, 1) # inclusion prob
  for(i in 1:M){
   z[i]~dbern(psi)
   s[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
   s[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
   dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1,site[i]])^2 + (s[i,2]-X[1:J,2,site[i]])^2)
   lam[i,1:J] <- lam0[site[i]]*K*exp(-dist[i,1:J]^2/(2*sigma^2))
  } # i individuals
  # use zeros trick for individuals to speed up the computation
  for(i in 1:n0){
    for(j in 1:J){
        y[i,j] ~ dpois(lam[i,j])
    } # j traps
  } # i individuals
  for(i in (n0+1):M){
    zeros[i] ~ dpois(sum(lam[i,1:J])*z[i])
  } # individuals
  N <- sum(z[1:M])
  D <- N/A
  })
}else
if(enc_dist == "binomial" & sex_sigma  == TRUE){
 scrcode <- nimble::nimbleCode({
  for(g in 1:nSites){  
      p0[g] ~ dunif(0,1) # baseline encounter probability
  } # g sites
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
  p[i,1:J] <- p0[site[i]]*exp(-dist[i,1:J]^2/(2*sigma[sx[i]]^2))
  } # i individuals
  # use zeros trick for individuals to speed up the computation
 for(i in 1:n0){
  for(j in 1:J){
   y[i,j] ~ dbin(p[i,j],K)
  } # j traps
 } # i individuals
 for(i in (n0+1):M){
  zeros[i] ~ dbern((1 - prod(1 - p[i,1:J])^K)*z[i])
 } # i individuals
  N <- sum(z[1:M])
  D <- N/A
 })
} else
 if(enc_dist == "poisson" & sex_sigma  == TRUE){
  scrcode <- nimble::nimbleCode({
   for(g in 1:nSites){  
    lam0[g] ~ dunif(0,lam0_upper) # baseline encounter rate
   }
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
   lam[i,1:J,k] <- lam0[site[i]]*K*exp(-dist[i,1:J]^2/(2*sigma[sx[i]]^2))
  } # i marked individuals
 # use zeros trick for marked individuals to speed up the computation
  for(i in 1:n0){
    for(j in 1:J){
      y[i,j] ~ dpois(lam[i,j])
    } # j traps
  } # i individuals
  for(i in (n0+1):M){
     zeros[i] ~ dpois(sum(lam[i,1:J])*z[i])
  } # individuals
    N <- sum(z[1:M])
    D <- N/A
  })
}
  return(scrcode)
} else  # End 2D models
if(dim_y == 3){
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
  for(i in 1:n0){
    for(j in 1:J){
      for(k in 1:K){
        y[i,j,k] ~ dpois(lam[i,j,k])
      } # occasions
    } # j traps
  } # i individuals
 for(i in (n0+1):M){
  zeros[i] ~ dpois(sum(lam[i,1:J,1:K])*z[i])
 } # i individuals
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
   for(g in 1:nSites){
    lam0[g] ~ dunif(0,lam0_upper) # baseline encounter rate
   } # g sites
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
    lam[i,1:J,k] <- lam0[site[i]]*exp(-dist[i,1:J]^2/(2*sigma[sx[i]]^2))
  } # k occasions
} # i individuals
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
     } # i individuals
    N <- sum(z[1:M])
    D <- N/A
  })
}
return(scrcode)
} # end 3D model
} else
if(hab_mask==TRUE){
if(dim_y == 2){
if(enc_dist == "binomial" & sex_sigma  == FALSE){
scrcode <- nimble::nimbleCode({
for(g in 1:nSites){  
p0[g] ~ dunif(0,1) # baseline encounter probability
} # g sites
sigma ~ dunif(0, sigma_upper) # scaling parameter
sigma.pixel <- sigma / pixelWidth # scaled for habitat mask
psi ~ dunif(0, 1) # inclusion prob
for(i in 1:M){
z[i]~dbern(psim[i])
  # adjust psi for the proportion of available habitat at each site
psim[i] <- (1-(1-psi)^prop.habitat[site[i]]) 
s[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
s[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1),site[i]] # habitat check
OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1,site[i]])^2 + (s[i,2]-X[1:J,2,site[i]])^2)
p[i,1:J] <- p0[site[i]]*exp(-dist[i,1:J]^2/(2*sigma.pixel^2))
} # i individuals
# use zeros trick for individuals to speed up the computation
for(i in 1:n0){
for(j in 1:J){
  y[i,j] ~ dbin(p[i,j],K)
} # j traps
} # i individuals
for(i in (n0+1):M){
zeros[i] ~ dbern((1 - prod(1 - p[i,1:J])^K)*z[i])
} # i individuals
N <- sum(z[1:M])
D <- N/A
})
} else
if(enc_dist == "poisson" & sex_sigma  == FALSE){
scrcode <- nimble::nimbleCode({
for(g in 1:nSites){  
lam0[g] ~ dunif(0,lam0_upper) # baseline encounter rate
}
sigma ~ dunif(0, sigma_upper) # scaling parameter
sigma.pixel <- sigma / pixelWidth # scaled for habitat mask
psi ~ dunif(0, 1) # inclusion prob
for(i in 1:M){
z[i]~dbern(psim[i])
# adjust psi for the proportion of available habitat at each site
psim[i] <- (1-(1-psi)^prop.habitat[site[i]]) 
s[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
s[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1),site[i]] # habitat check
OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1,site[i]])^2 + (s[i,2]-X[1:J,2,site[i]])^2)
lam[i,1:J] <- lam0[site[i]]*K*exp(-dist[i,1:J]^2/(2*sigma.pixel^2))
} # i individuals
# use zeros trick for individuals to speed up the computation
for(i in 1:n0){
for(j in 1:J){
    y[i,j] ~ dpois(lam[i,j])
} # j traps
} # i individuals
for(i in (n0+1):M){
zeros[i] ~ dpois(sum(lam[i,1:J])*z[i])
} # individuals
N <- sum(z[1:M])
D <- N/A
})
}else
if(enc_dist == "binomial" & sex_sigma  == TRUE){
scrcode <- nimble::nimbleCode({
for(g in 1:nSites){  
 p0[g] ~ dunif(0,1) # baseline encounter probability
} # g sites
psi_sex ~ dunif(0,1) # probability sex = 1
sigma[1] ~ dunif(0, sigma_upper) # scaling parameter, sex = 0
sigma[2] ~ dunif(0, sigma_upper) # scaling parameter, sex = 1
sigma.pixel[1] <- sigma[1] / pixelWidth # scaled for habitat mask 
sigma.pixel[2] <- sigma[2] / pixelWidth # scaled for habitat mask              
psi ~ dunif(0, 1) # inclusion prob
for(i in 1:M){
  sex[i] ~ dbern(psi_sex)
  sx[i] <- sex[i] + 1
  z[i]~dbern(psim[i])
  # adjust psi for the proportion of available habitat at each site 
  psim[i] <- (1-(1-psi)^prop.habitat[site[i]])
  s[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
  s[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
  pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1),site[i]] # habitat check
  OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
  dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1,site[i]])^2 + (s[i,2]-X[1:J,2,site[i]])^2)
  p[i,1:J] <- p0[site[i]]*exp(-dist[i,1:J]^2/(2*sigma.pixel[sx[i]]^2))
} # i individuals
# use zeros trick for individuals to speed up the computation
for(i in 1:n0){
  for(j in 1:J){
      y[i,j] ~ dbin(p[i,j],K)
  } # j traps
} # i individuals
for(i in (n0+1):M){
    zeros[i] ~ dbern((1 - prod(1 - p[i,1:J])^K)*z[i])
} # i individuals
N <- sum(z[1:M])
D <- N/A
})
} else
if(enc_dist == "poisson" & sex_sigma  == TRUE){
scrcode <- nimble::nimbleCode({
for(g in 1:nSites){ 
  lam0[g] ~ dunif(0,lam0_upper) # baseline encounter rate
} # g sites 
  psi_sex ~ dunif(0,1) # probability sex = 1
  sigma[1] ~ dunif(0, sigma_upper) # scaling parameter, sex = 0
  sigma[2] ~ dunif(0, sigma_upper) # scaling parameter, sex = 1
  sigma.pixel[1] <- sigma[1] / pixelWidth # scaled for habitat mask 
  sigma.pixel[2] <- sigma[2] / pixelWidth # scaled for habitat mask   
  psi ~ dunif(0, 1) # inclusion prob
  for(i in 1:M){
  sex[i] ~ dbern(psi_sex)
  sx[i] <- sex[i] + 1
  z[i]~dbern(psim[i])
   # adjust psi for the proportion of available habitat at each site
  psim[i] <- (1-(1-psi)^prop.habitat[site[i]])
  s[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
  s[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
  pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1),site[i]] # habitat check
  OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
  dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1,site[i]])^2 + (s[i,2]-X[1:J,2,site[i]])^2)
  lam[i,1:J,k] <- lam0[site[i]]*K*exp(-dist[i,1:J]^2/(2*sigma.pixel[sx[i]]^2))
  } # i marked individuals
  # use zeros trick for marked individuals to speed up the computation
   for(i in 1:n0){
    for(j in 1:J){
        y[i,j] ~ dpois(lam[i,j])
    } # j traps
  } # i individuals
  for(i in (n0+1):M){
    zeros[i] ~ dpois(sum(lam[i,1:J,1:K])*z[i])
  } # individuals
  N <- sum(z[1:M])
  D <- N/A
})
}
return(scrcode)
} else  # End 2D models
if(dim_y == 3){
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
   # adjust psi for the proportion of available habitat at each site
  psim[i] <- (1-(1-psi)^prop.habitat[site[i]]) 
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
  sigma ~ dunif(0, sigma_upper) # scaling parameter
  sigma.pixel <- sigma / pixelWidth # scaled for habitat mask
  psi ~ dunif(0, 1) # inclusion prob
  for(g in 1:nSites){
    lam0[g] ~ dunif(0,lam0_upper) # site-specific baseline encounter rate
  } # g sites
  for(i in 1:M){
    z[i]~dbern(psim[i])
   # adjust psi for the proportion of available habitat at each site
  psim[i] <- (1-(1-psi)^prop.habitat[site[i]]) 
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
   # adjust psi for the proportion of available habitat at each site
  psim[i] <- (1-(1-psi)^prop.habitat[site[i]])
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
    for(g in 1:nSites){  
      lam0[g] ~ dunif(0,lam0_upper) # baseline encounter rate
    } # g sites  
psi_sex ~ dunif(0,1) # probability sex = 1
sigma[1] ~ dunif(0, sigma_upper) # scaling parameter, sex = 0
sigma[2] ~ dunif(0, sigma_upper) # scaling parameter, sex = 1
sigma.pixel[1] <- sigma[1] / pixelWidth # scaled for habitat mask
sigma.pixel[2] <- sigma[2] / pixelWidth # scaled for habitat mask
psi ~ dunif(0, 1) # inclusion prob
for(i in 1:M){
  z[i]~dbern(psim[i])
   # adjust psi for the proportion of available habitat at each site
  psim[i] <- (1-(1-psi)^prop.habitat[site[i]])
  sex[i] ~ dbern(psi_sex)
  sx[i] <- sex[i] + 1
  s[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
  s[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
  pOK[i] <- hab_mask[(trunc(s[i,2])+1),(trunc(s[i,1])+1),site[i]] # habitat check
  OK[i] ~ dbern(pOK[i]) # OK[i] <- 1, the ones trick
  dist[i,1:J] <- sqrt((s[i,1]-X[1:J,1,site[i]])^2 + (s[i,2]-X[1:J,2,site[i]])^2)
  for(k in 1:K){
    lam[i,1:J,k] <- lam0[site[i]]*exp(-dist[i,1:J]^2/(2*sigma.pixel[sx[i]]^2))
  } # k occasions
} # i individuals
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
} # end 3D model
    } # end models with hab_mask
  } # end trapsClustered
} # End function 'get_classic"



#' Function to generate starting locations for activity centers
#'
#' @description Generate a matrix of initial starting locations, possibly accounting 
#' for habitat mask.
#'
#' @param y Either a matrix or array of encounter history data, possiblity 
#' from \code{sim_classic()}.
#' @param M An integer of the total augmented population size (i.e., detected 
#' and augmented individuals). UTMs. An array is used when traps are clustered 
#' over a survey area.
#' @param X Either a matrix or array representing the coordinates of traps in
#' UTMs. An array is used when traps are clustered over a survey area.
#' @param buff The distance (m or km) that the traps should be
#' buffered by as an integer. This is typically 3 times the sigma parameter.
#' @param site Either \code{NULL} (if a 2D trap array is used) or a vector of 
#' integers denoting which trap array an individual (either detected or 
#' augmented) belongs to. Note that \code{site} is provided from 
#' \code{\link{sim_classic}} when a 3D trap array is used. However, this
#'  \code{site} variable must be correctly augmented based on the total 
#'  augmented population size (i.e., \code{M}).
#' @param hab_mask Either \code{FALSE} (the default) or a matrix or array output
#'  from \code{\link{mask_polygon}} or \code{\link{mask_raster}} functions.
#' @param all_random Logical. If \code{TRUE}, then encounter data \code{y} are 
#' ignored and all initial activity center starting locations are randomly chosen.
#' If \code{FALSE} (the default), then initial values will be the mean capture 
#' location for detected individuals and random locations for augmented individuals. 
#' @return A matrix of initial activity center coordinates with \code{M} rows 
#' and 2 columns.
#' @details This function generates initial activity center locations based 
#' on encounter histories, augmented population size, state-space buffer, 
#' and potentially a habitat mask. Note that mean trap detection locations 
#' are used for detected individuals while intial values are randomly drawn 
#' for augemented individuals. Also, a habitat check will be conducted for 
#' all locations when a habitat mask is included.
#' @author Daniel Eacker
#' @examples
#' # simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' 
#' # add some random noise to locations
#' traps <- traps + runif(prod(dim(traps)),-20,20) 
#'
#' mysigma = 300 # simulate sigma of 300 m
#' mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#'
#' # Create grid and extent
#' Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = 100)
#'
#' # simulate SCR data
#' data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, 
#' sigma_ = mysigma, prop_sex = 1, N = 200, K = 4, base_encounter = 0.25, 
#' enc_dist = "binomial", hab_mask = FALSE, setSeed = 50)
#'
#' # generate initial activity center coordinates for 2D trap array without 
#' # habitat mask
#' s.st3d = initialize_classic(y=data3d$y, M=500, X=traps, buff = 3*mysigma, 
#' hab_mask = FALSE, all_random = FALSE)
#'
#' # make simple plot
#' par(mfrow=c(1,1))
#' plot(Grid$grid, pch=20,ylab="Northing",xlab="Easting")
#' points(traps, col="blue",pch=20)
#' points(s.st3d, col="red",pch=20)
#' @name initialize_classic
#' @export
initialize_classic <- function(y, M, X, buff, site, hab_mask, all_random = FALSE){
  if(length(dim(X))!=2 & length(dim(X))!=3){
    stop("Trapping grid must be only 2 or 3 dimensions")
  }
  if(isFALSE(all_random)){
  # for dim of length 2
  if(length(dim(X))==2){
    n0 <- length(which(apply(y,1,sum)!=0))
    s.st <- matrix(NA, nrow=M, ncol=2)
    # create x limits for state-space
    xlim <- c(min(X[,1] - max(buff)), max(X[,1] + max(buff)))
     # create y limits for state-space
    ylim <- c(min(X[,2] - max(buff)), max(X[,2] + max(buff)))
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
        sx.rescale <- scales::rescale(s.st[i,1], to = c(0,dim(hab_mask)[2]), 
                                      from=xlim)
        sy.rescale <- scales::rescale(s.st[i,2], to = c(0,dim(hab_mask)[1]), 
                                      from=ylim)
        pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1)]
         # if their mean activity center is not in suitable habitat, 
         # then randomly sample
        while(pOK==0){
          s.st[i,1]<-runif(1,xlim[1],xlim[2])
          s.st[i,2]<-runif(1,ylim[1],ylim[2])
          sx.rescale <- scales::rescale(s.st[i,1], to = c(0,dim(hab_mask)[2]), 
                                        from=xlim)
          sy.rescale <- scales::rescale(s.st[i,2], to = c(0,dim(hab_mask)[1]), 
                                        from=ylim)
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
          sx.rescale <- scales::rescale(s.st[i,1], to = c(0,dim(hab_mask)[2]), 
                                        from=xlim)
          sy.rescale <- scales::rescale(s.st[i,2], to = c(0,dim(hab_mask)[1]), 
                                        from=ylim)
          pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1)]
          while(pOK==0){
            s.st[i,1]<-runif(1,xlim[1],xlim[2])
            s.st[i,2]<-runif(1,ylim[1],ylim[2])
            sx.rescale <- scales::rescale(s.st[i,1], to = c(0,dim(hab_mask)[2]), 
                                          from=xlim)
            sy.rescale <- scales::rescale(s.st[i,2], to = c(0,dim(hab_mask)[1]),
                                          from=ylim)
            pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1)]
          }
        } # end habitat check
    } # augmented individuals
  }  else
  if(length(dim(X))==3){
    if(length(dim(X))==3 & length(site)!=M){
    stop("Include 'site' variable with length M")
    }
    n0 <- length(which(apply(y,1,sum)!=0))
    s.st <- matrix(NA, nrow=M, ncol=2)
    for(i in 1:n0){
      # create x limits for state-space
    xlim <- c(min(X[,1,site[i]] - max(buff)), max(X[,1,site[i]] + max(buff))) 
    # create y limits for state-space
    ylim <- c(min(X[,2,site[i]] - max(buff)), max(X[,2,site[i]] + max(buff))) 
    temp.y <- y[i,,]
    temp.X <- X[which(apply(temp.y,1,sum)!=0),,site[i]]
    ntimes <- apply(temp.y,1,sum)[which(apply(temp.y,1,sum)!=0)]
    df <- data.frame(temp.X,ntimes)
    temp.X2 = df[rep(seq_len(nrow(df)), df$ntimes),1:2]
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
          sx.rescale <- scales::rescale(s.st[i,1], 
                        to = c(0,dim(hab_mask)[2]), from=xlim)
          sy.rescale <- scales::rescale(s.st[i,2], 
                        to = c(0,dim(hab_mask)[1]), from=ylim)
          pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1),
                          site[i]]
  # if their mean activity center is not in suitable habitat, 
  #then randomly sample
        while(pOK==0){
          s.st[i,1]<-runif(1,xlim[1],xlim[2])
          s.st[i,2]<-runif(1,ylim[1],ylim[2])
          sx.rescale <- scales::rescale(s.st[i,1], 
                                        to = c(0,dim(hab_mask)[2]), from=xlim)
          sy.rescale <- scales::rescale(s.st[i,2], 
                                        to = c(0,dim(hab_mask)[1]), from=ylim)
          pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1),
                          site[i]]
        }
      } # end habitat check
  } # end detected individuals
  for(i in (n0+1):M){
        xlim <- c(min(X[,1,site[i]] - max(buff)), max(X[,1,site[i]] + max(buff))) # create x limits for state-space
        ylim <- c(min(X[,2,site[i]] - max(buff)), max(X[,2,site[i]] + max(buff))) # create y limits for state-space
    if(isFALSE(hab_mask)){ # check for habitat mask for augmented individuals
        s.st[i,1]<-runif(1,xlim[1],xlim[2])
        s.st[i,2]<-runif(1,ylim[1],ylim[2])
    } else
      if(isFALSE(hab_mask)==FALSE){
        s.st[i,1] <-runif(1,xlim[1],xlim[2])
        s.st[i,2] <-runif(1,ylim[1],ylim[2])
        sx.rescale <- scales::rescale(s.st[i,1], to = c(0,dim(hab_mask)[2]), 
                                      from=xlim)
        sy.rescale <- scales::rescale(s.st[i,2], to = c(0,dim(hab_mask)[1]), 
                                      from=ylim)
        pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1),site[i]]
        while(pOK==0){
        s.st[i,1] <-runif(1,xlim[1],xlim[2])
        s.st[i,2] <-runif(1,ylim[1],ylim[2])
        sx.rescale <- scales::rescale(s.st[i,1], to = c(0,dim(hab_mask)[2]), 
                                      from=xlim)
        sy.rescale <- scales::rescale(s.st[i,2], to = c(0,dim(hab_mask)[1]), 
                                      from=ylim)
        pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1),site[i]]
        }
      } # end habitat check
   } # end augmented individuals
 } # end 3D initialize
} else # end all_random = FALSE    
 if(isFALSE(all_random) == FALSE){ 
   # for dim of length 2
  if(length(dim(X))==2){
    s.st <- matrix(NA, nrow=M, ncol=2)
    # create x limits for state-space
    xlim <- c(min(X[,1] - max(buff)), max(X[,1] + max(buff)))
     # create y limits for state-space
    ylim <- c(min(X[,2] - max(buff)), max(X[,2] + max(buff)))
    for(i in 1:M){
      if(isFALSE(hab_mask)){ # check for habitat mask for augmented individuals
        s.st[i,1]<-runif(1,xlim[1],xlim[2])
        s.st[i,2]<-runif(1,ylim[1],ylim[2])
      } else
        if(isFALSE(hab_mask)==FALSE){
          s.st[i,1]<-runif(1,xlim[1],xlim[2])
          s.st[i,2]<-runif(1,ylim[1],ylim[2])
          sx.rescale <- scales::rescale(s.st[i,1], to = c(0,dim(hab_mask)[2]), 
                                        from=xlim)
          sy.rescale <- scales::rescale(s.st[i,2], to = c(0,dim(hab_mask)[1]), 
                                        from=ylim)
          pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1)]
          while(pOK==0){
            s.st[i,1]<-runif(1,xlim[1],xlim[2])
            s.st[i,2]<-runif(1,ylim[1],ylim[2])
            sx.rescale <- scales::rescale(s.st[i,1], to = c(0,dim(hab_mask)[2]), 
                                          from=xlim)
            sy.rescale <- scales::rescale(s.st[i,2], to = c(0,dim(hab_mask)[1]),
                                          from=ylim)
            pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1)]
          }
        } # end habitat check
    } # augmented individuals
  }  else
  if(length(dim(X))==3){
    if(length(dim(X))==3 & length(site)!=M){
    stop("Include 'site' variable with length M")
    }
    s.st <- matrix(NA, nrow=M, ncol=2)
  for(i in 1:M){
        xlim <- c(min(X[,1,site[i]] - max(buff)), max(X[,1,site[i]] + max(buff))) # create x limits for state-space
        ylim <- c(min(X[,2,site[i]] - max(buff)), max(X[,2,site[i]] + max(buff))) # create y limits for state-space
    if(isFALSE(hab_mask)){ # check for habitat mask for augmented individuals
        s.st[i,1]<-runif(1,xlim[1],xlim[2])
        s.st[i,2]<-runif(1,ylim[1],ylim[2])
    } else
      if(isFALSE(hab_mask)==FALSE){
        s.st[i,1] <-runif(1,xlim[1],xlim[2])
        s.st[i,2] <-runif(1,ylim[1],ylim[2])
        sx.rescale <- scales::rescale(s.st[i,1], to = c(0,dim(hab_mask)[2]), 
                                      from=xlim)
        sy.rescale <- scales::rescale(s.st[i,2], to = c(0,dim(hab_mask)[1]), 
                                      from=ylim)
        pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1),site[i]]
        while(pOK==0){
        s.st[i,1] <-runif(1,xlim[1],xlim[2])
        s.st[i,2] <-runif(1,ylim[1],ylim[2])
        sx.rescale <- scales::rescale(s.st[i,1], to = c(0,dim(hab_mask)[2]), 
                                      from=xlim)
        sy.rescale <- scales::rescale(s.st[i,2], to = c(0,dim(hab_mask)[1]), 
                                      from=ylim)
        pOK <- hab_mask[(trunc(sy.rescale)+1),(trunc(sx.rescale)+1),site[i]]
        }
      } # end habitat check
   } # end augmented individuals
 } # end 3D initialize
} # end all_random = TRUE  
  return(s.st)
} # end function 'initialize_classic'



#' Function to rescale trap coordinates, grid extent, and starting activity 
#' center coordinates
#'
#' @description Rescale inputs to prepare data for habitat mask to be used.
#'
#' @param X Either a matrix or array representing the coordinates of traps in
#' UTMs. An array is used when traps are clustered over a survey area.
#' @param ext An \code{Extent} object from the \code{raster} package. This is 
#' returned from \code{\link{grid_classic}}.
#' @param s.st A matrix of starting activity center coordinates. This is 
#' returned from \code{\link{initialize_classic}}
#' buffered by as an integer. This is typically 3 times the sigma parameter.
#' @param site Either \code{NULL} (if a 2D trap array is used) or a vector of 
#' integers denoting which trap array an individual (either detected or 
#' augmented) belongs to. Note that \code{site} is provided from 
#' \code{\link{sim_classic}} when a 3D trap array is used. However, this 
#' \code{site} variable must be correctly augmented based on the total 
#' augmented population size (i.e., \code{M}).
#' @param hab_mask A matrix or arary output from \code{\link{mask_polygon}} or
#'  \code{\link{mask_raster}} functions.
#' @return A list of rescaled trap coordinates, grid extents, and starting 
#' activity center coordinates.
#' @details This function is only meant to be used when habitat masking is 
#' incorporated into the model. The functions properly rescales inputs based on 
#' the dimensions of the habitat mask. Note that the \code{pixelWidth} needs to 
#' be included as an input in the model after inputs are rescaled to correctly
#'  estimate the scaling parameter (i.e., 'sigma').
#' @author Daniel Eacker
#' @importFrom sf st_polygon
#' @seealso \code{\link{mask_polygon}}, \code{\link{mask_raster}}
#' @examples
#'# simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' 
#' # add some random noise to locations
#' traps <- traps + runif(prod(dim(traps)),-20,20) 
#'
#' mysigma = 300 # simulate sigma of 300 m
#' mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#'
#' # create grid
#' Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = 100)
#'
#' # create polygon to use as a mask
#' library(sf)
#' poly = st_sfc(st_polygon(x=list(matrix(c(-1765,-1765,1730,-1650,1600,1650,
#' 0,1350,-800,1700,-1850,1000,-1765,-1765),ncol=2, byrow=TRUE))), crs =  mycrs)
#'
#' # make simple plot
#' par(mfrow=c(1,1))
#' plot(Grid$grid, pch=20)
#' points(traps, col="blue",pch=20)
#' plot(poly, add=TRUE)
#'
#' # create habitat mask
#' hab_mask = mask_polygon(poly = poly, grid = Grid$grid, crs_ = mycrs, 
#' prev_mask = NULL)
#'
#' # simulate SCR data
#' data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, 
#' sigma_ = mysigma, prop_sex = 1, N = 200, K = 4, base_encounter = 0.25, 
#' enc_dist = "binomial", hab_mask = hab_mask, setSeed = 50)
#'
#' # generate initial activity center coordinates for 2D trap array without 
#' #habitat mask
#' s.st3d = initialize_classic(y=data3d$y, M=500, X=traps, buff = 3*mysigma, 
#' hab_mask = hab_mask)
#'
#' # build rescaled constants list for 2D trap array
#' constList = rescale_classic(X = traps, ext = Grid$ext, s.st = s.st3d, 
#' site = NULL, hab_mask = hab_mask)
#' str(constList)
#' @name rescale_classic
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
    rescale_list$X[,1] <- scales::rescale(X[,1], to = c(0,dim(hab_mask)[2]), 
                                          from=ext[1:2])
    rescale_list$X[,2] <- scales::rescale(X[,2], to = c(0,dim(hab_mask)[1]), 
                                          from=ext[3:4])
    rescale_list$s.st[,1] <- scales::rescale(s.st[,1], 
                              to = c(0,dim(hab_mask)[2]), from=ext[1:2])
    rescale_list$s.st[,2] <- scales::rescale(s.st[,2], 
                              to = c(0,dim(hab_mask)[1]), from=ext[3:4])
    rescale_list$ext[1:2] <- scales::rescale(ext[1:2], 
                              to = c(0,dim(hab_mask)[2]), from=ext[1:2])
    rescale_list$ext[3:4] <- scales::rescale(ext[3:4], 
                              to = c(0,dim(hab_mask)[1]), from=ext[3:4])

  }else
if(length(dim(X))==3){
  for(g in 1:dim(X)[3]){ # loop over trap arrays
  rescale_list$X[,1,g] <- scales::rescale(X[,1,g], 
                  to = c(0,dim(hab_mask)[2]), from=ext[[g]][1:2])
  rescale_list$X[,2,g] <- scales::rescale(X[,2,g], 
                  to = c(0,dim(hab_mask)[1]), from=ext[[g]][3:4])
  rescale_list$s.st[which(site==g),1] <- scales::rescale(s.st[which(site==g),1], 
                 to = c(0,dim(hab_mask)[2]), from=ext[[g]][1:2])
  rescale_list$s.st[which(site==g),2] <- scales::rescale(s.st[which(site==g),2], 
                  to = c(0,dim(hab_mask)[1]), from=ext[[g]][3:4])
    # check that starting coordinates are still in bounds
  rescale_list$s.st[which(site==g),1] <- 
      ifelse(rescale_list$s.st[which(site==g),1]>=dim(hab_mask)[2],
      rescale_list$s.st[which(site==g),1]-1,rescale_list$s.st[which(site==g),1])
  rescale_list$s.st[which(site==g),2] <- 
      ifelse(rescale_list$s.st[which(site==g),2]>=dim(hab_mask)[1],
      rescale_list$s.st[which(site==g),2]-1,rescale_list$s.st[which(site==g),2])
  rescale_list$ext[[g]][1:2] <- scales::rescale(ext[[g]][1:2], 
       to = c(0,dim(hab_mask)[2]), from=ext[[g]][1:2])
  rescale_list$ext[[g]][3:4] <- scales::rescale(ext[[g]][3:4], 
      to = c(0,dim(hab_mask)[1]), from=ext[[g]][3:4])
  }
    }
  return(rescale_list)
} # end function 'rescale_classic'



#' Function to create habitat mask matrix or array from polygon
#'
#' @description Creates a matrix or array to use as a habitat mask to account for unsuitable 
#' habitat.
#'
#' @param poly A polygon created using the \code{sf} package of class 
#' \code{"sfc_POLYGON"}
#' @param grid A matrix or array object of the the state-space grid. This 
#' is returned from \code{\link{grid_classic}}.
#' @param crs_ The UTM coordinate reference system (EPSG code) used for your
#' location provided as an integer (e.g., 32608 for WGS 84/UTM Zone 8N).
#' @param prev_mask Either \code{NULL} or a previously created habitat mask 
#' matrix or array from \code{\link{mask_polygon}} or \code{\link{mask_raster}}. 
#' This allows for habitat masks to be combined to account for different 
#' spatial features.
#' @return A matrix or array of 0's and 1's denoting unsuitable and suitable 
#' habitat respectively.
#' @details This function creates a habitat matrix or array depending upon 
#' whether a 2D (former) or 3D (latter) trap array is used. This matrix can be 
#' directly included as data in Bayesian SCR models run using \code{nimble}.
#' @author Daniel Eacker
#' @seealso \code{\link{mask_raster}}
#' @examples
#' # simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' 
#' # add some random noise to locations
#' traps <- traps + runif(prod(dim(traps)),-20,20) 
#'
#' mysigma = 300 # simulate sigma of 300 m
#' mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#'
#' # create state-space grid and extent
#' Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = 100)
#'
#' # create polygon to use as a mask
#' library(sf)
#' poly = st_sfc(st_polygon(x=list(matrix(c(-1765,-1765,1730,-1650,1600,1650,
#' 0,1350,-800,1700,-1850,1000,-1765,-1765),ncol=2, byrow=TRUE))), crs =  mycrs)
#'
#' # make simple plot
#' par(mfrow=c(1,2))
#' plot(Grid$grid, pch=20)
#' points(traps, col="blue",pch=20)
#' plot(poly, add=TRUE)
#'
#' # create habitat mask from polygon
#' hab_mask = mask_polygon(poly = poly, grid = Grid$grid, crs_ = mycrs, 
#' prev_mask = NULL)
#'
#' # make simple plot
#' library(raster)
#' plot(raster(apply(hab_mask,2,rev)))
#'
#' # make simple plot
#' poly2 = st_sfc(st_polygon(x=list(matrix(c(-1365,-1365,1730,-1650,1500,1550,
#' 0,1350,-800,1700,-1850,1000,-1365,-1365),ncol=2, byrow=TRUE))), crs =  mycrs)
#' plot(poly2, add=TRUE)
#'
#' # mask second polygon, building on previous habitat mask
#' hab_mask2 = mask_polygon(poly = poly2, grid = Grid$grid, crs_ = mycrs, 
#' prev_mask = hab_mask)
#'
#' # make simple plot
#' plot(Grid$grid, pch=20)
#' points(traps, col="blue",pch=20)
#' plot(poly, add=TRUE)
#' plot(poly2, add=TRUE)
#' plot(raster(apply(hab_mask2,2,rev)))
#'
#' # create an array of traps, as an approach where individuals will only be 
#' # detected at one of the trap arrays (e.g., Furnas et al. 2018)
#' Xarray = array(NA, dim=c(nrow(traps),2,2))
#' Xarray[,,1]=traps
#' Xarray[,,2]=traps+4000 # shift trapping grid to new locations
#'
#' # Example of using habitat mask with 3D trap array (need polygon that 
#' # masks both trapping extents)
#' GridX = grid_classic(X = Xarray, crs_ = mycrs, buff = 3*mysigma, res = 100)
#'
#' # make simple plot
#' par(mfrow=c(1,1))
#' plot(GridX$grid[,,1],xlim=c(-1600,6000),ylim=c(-1600,6000),col="darkgrey",
#' pch=20,ylab="Northing",xlab="Easting")
#' points(Xarray[,,1],col="blue",pch=20)
#' points(GridX$grid[,,2],pch=20,col="darkgrey")
#' points(Xarray[,,2],col="blue",pch=20)
#'
#' poly = st_sfc(st_polygon(x=list(matrix(c(-1660,-1900,5730,-1050,5470,5150,
#' 0,6050,-1800,5700,-1660,-1900),ncol=2, byrow=TRUE))), crs =  mycrs)
#' plot(poly, add=TRUE)
#'
#' # get 3D habitat mask array for 3D grid
#' hab_mask = mask_polygon(poly = poly, grid = GridX$grid, crs_ = mycrs, 
#' prev_mask = NULL)
#'
#' # make simple plot
#' par(mfrow=c(1,2))
#' apply(hab_mask,3,function(x) plot(raster(apply(x,2,rev))))
#' @name mask_polygon
#' @export
mask_polygon <- function(poly, grid, crs_, prev_mask){
# need to determine if X is 2 or 3 dimensional (stop if not)
if(length(dim(grid))!=2 & length(dim(grid))!=3){
stop("Trapping grid must be only 2 or 3 dimensions")
}
# for dim of length 2
if(length(dim(grid))==2){
grid_pts <- sf::st_cast(sf::st_sfc(sf::st_multipoint(grid), crs =  crs_),
                        "POINT")
habitat_mask <- apply(matrix(ifelse(is.na(as.numeric(sf::st_intersects(grid_pts,
    poly))),0,as.numeric(sf::st_intersects(grid_pts,poly))),
  nrow=dim(raster::rasterFromXYZ(grid,crs=crs_))[1],
    ncol=dim(raster::rasterFromXYZ(grid,crs=crs_))[2], byrow=TRUE),2,rev)
# to combine with a previous habitat mask
if(is.null(prev_mask)==FALSE){
  # Check to see if dimensions of previous and current habitat
  if(all(dim(habitat_mask) == dim(prev_mask))==FALSE){
    stop("Dimension of previous habitat matrix does not equal 
         dimension of current one")
  }
  habitat_mask <- habitat_mask * prev_mask
}
}else
# for dim of length 3
if(length(dim(grid))==3){
  grid_pts <- apply(grid, 3, function(x) sf::st_cast(sf::st_sfc(
    sf::st_multipoint(x), crs =  crs_),"POINT"))
  habitat_mask <- lapply(grid_pts, function(x) apply(matrix(ifelse(
    is.na(as.numeric(sf::st_intersects(x,poly))),0,
     as.numeric(sf::st_intersects(x,poly))),nrow=dim(raster::rasterFromXYZ(
       grid[,,1],crs=crs_))[1], ncol=dim(raster::rasterFromXYZ(grid[,,1],
                            crs=crs_))[2], byrow=TRUE),2,rev))
  habitat_mask <- array(unlist(habitat_mask),dim=c(dim(habitat_mask[[1]]),
                        dim(grid)[3])) # convert to array
  # to combine with a previous habitat mask
  if(is.null(prev_mask)==FALSE){
    # Check to see if dimensions of previous and current habitat
    if(all(dim(habitat_mask) == dim(prev_mask))==FALSE){
      stop("Dimension of previous habitat matrix does not equal 
           dimension of current one")
    }
    habitat_mask <- habitat_mask * prev_mask
  }
}
return(habitat_mask)
} # End function 'mask_polygon'





#' Function to create habitat mask matrix or array from raster
#'
#' @description Creates a matrix or array to use as a habitat mask to account for unsuitable 
#' habitat.
#'
#' @param rast A raster layer created using the \code{raster} package of class 
#' \code{"RasterLayer"}
#' @param FUN A function that defines the criteria for suitable habitat.
#' @param grid A matrix or array object of the the state-space grid. This is 
#' returned from \code{\link{grid_classic}}.
#' @param crs_ The UTM coordinate reference system (EPSG code) used for your
#' location provided as an integer (e.g., 32608 for WGS 84/UTM Zone 8N).
#' @param prev_mask Either \code{NULL} or a previously created habitat mask 
#' matrix or array from \code{\link{mask_polygon}} or \code{\link{mask_raster}}. 
#' This allows for habitat masks to be combined to account for different spatial
#'  features.
#' @return A matrix or array of 0's and 1's denoting unsuitable and suitable
#'  habitat respectively.
#' @details This function creates a habitat matrix or array depending upon 
#' whether a 2D (former) or 3D (latter) trap array is used. This matrix can be 
#' directly included as data in Bayesian SCR models run using \code{nimble}.
#' @author Daniel Eacker
#' @importFrom sp CRS proj4string
#' @importFrom sf st_intersects
#' @importFrom raster proj4string extract
#' @importFrom methods as
#' @seealso \code{\link{mask_polygon}}
#' @examples
#' # simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' 
#'  # add some random noise to locations
#' traps <- traps + runif(prod(dim(traps)),-20,20)
#'
#' mysigma = 300 # simulate sigma of 300 m
#' mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#'
#' # create state-space grid and extent
#' Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = 100)
#'
#' # run previous code used for mask_polygon() to create raster for example
#' library(sf)
#' poly = st_sfc(st_polygon(x=list(matrix(c(-1665,-1665,1730,-1650,1600,1650,
#' 0,1350,-800,1700,-1850,1000,-1665,-1665),ncol=2, byrow=TRUE))), crs =  mycrs)
#' hab_mask = mask_polygon(poly = poly, grid = Grid$grid, crs_ = mycrs, 
#' prev_mask = NULL)
#'
#' # create raster for demonstration purposes
#' library(raster)
#' rast <- raster(nrow=dim(hab_mask)[1], ncol=dim(hab_mask)[2],ext=Grid$ext,
#' crs=mycrs)
#' rast[] = apply(hab_mask,2,rev)
#'
#' # create habitat mask using raster
#' hab_mask_r = mask_raster(rast = rast, FUN = function(x){x==1}, 
#' grid = Grid$grid, crs_ = mycrs, prev_mask = NULL)
#'
#' # make simple plot
#' # returns identical results as input rast (but this was just an example raster)
#' plot(raster(apply(hab_mask_r,2,rev))) 
#'
#' # create an array of traps, as an approach where individuals will only be 
#' # detected at one of the trap arrays (e.g., Furnas et al. 2018)
#' Xarray = array(NA, dim=c(nrow(traps),2,2))
#' Xarray[,,1]=traps
#' Xarray[,,2]=traps+4000 # shift trapping grid to new locations
#'
#' # create grid and extent for 3D trap array
#' GridX = grid_classic(X = Xarray, crs_ = mycrs, buff = 3*mysigma, res = 100)
#'
#' # make simple plot
#' par(mfrow=c(1,1))
#' plot(GridX$grid[,,1],xlim=c(-1600,6000),ylim=c(-1600,6000),col="darkgrey",
#' pch=20,ylab="Northing",xlab="Easting")
#' points(Xarray[,,1],col="blue",pch=20)
#' points(GridX$grid[,,2],pch=20,col="darkgrey")
#' points(Xarray[,,2],col="blue",pch=20)
#'
#' # create polygon to use as a mask and covert to raster
#' poly = st_sfc(st_polygon(x=list(matrix(c(-1660,-1900,5730,-1050,5470,5650,
#' 0,6050,-1800,5700,-1660,-1900),ncol=2, byrow=TRUE))), crs =  mycrs)
#'
#' # add polygon to plot
#' plot(poly, add=TRUE)
#'
#' # make raster from polygon
#' rast = raster(xmn=-2000, xmx=6000, ymn=-2000, ymx=6500,res=100,crs=mycrs)
#' rast[]=st_intersects(st_cast(st_sfc(st_multipoint(coordinates(rast)), 
#' crs =  mycrs),"POINT"),poly,sparse=FALSE)
#'
#' # make simple plot of raster
#' plot(rast)
#'
#' # get 3D habitat mask array for 3D grid
#' hab_mask = mask_raster(rast = rast, FUN = function(x){x==1},grid = GridX$grid, 
#' crs_ = mycrs, prev_mask = NULL)
#  #make plot
#' par(mfrow=c(1,2))
#' apply(hab_mask,3,function(x) plot(raster(apply(x,2,rev))))
#' @name mask_raster
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
    grid_pts <- sf::st_cast(sf::st_sfc(sf::st_multipoint(grid), crs =  crs_),
                            "POINT")
    vals <- raster::extract(rast,methods::as(grid_pts,"Spatial"))
    rast_ind <- FUN(vals)
    habitat_mask <- apply(matrix(as.numeric(rast_ind),nrow=dim(
      raster::rasterFromXYZ(grid,crs=crs_))[1],
      ncol=dim(raster::rasterFromXYZ(grid,crs=crs_))[2], byrow=TRUE),2,rev)
    # to combine with a previous habitat mask
    if(is.null(prev_mask)==FALSE){
      # Check to see if dimensions of previous and current habitat
      if(all(dim(habitat_mask) == dim(prev_mask))==FALSE){
        stop("Dimension of previous habitat matrix does not equal dimension 
             of current one")
      }
      habitat_mask <- habitat_mask * prev_mask
    }
  }else
    # for dim of length 3
    if(length(dim(grid))==3){
      grid_pts <- apply(grid, 3, function(x) sf::st_cast(sf::st_sfc(
        sf::st_multipoint(x), crs =  crs_),"POINT"))
      vals <- lapply(grid_pts,function(x) raster::extract(rast,
                              methods::as(x,"Spatial")))
      rast_ind <- lapply(vals, function(x) FUN(x))
      habitat_mask <- lapply(rast_ind, function(x) 
        apply(matrix(as.numeric(x),nrow=dim(raster::rasterFromXYZ(grid[,,1],
        crs=crs_))[1],
          ncol=dim(raster::rasterFromXYZ(grid[,,1],crs=crs_))[2], byrow=TRUE),
        2,rev))
      habitat_mask <- array(unlist(habitat_mask),dim=c(dim(habitat_mask[[1]]),
                                  dim(grid)[3])) # convert to array
      # to combine with a previous habitat mask
      if(is.null(prev_mask)==FALSE){
        # Check to see if dimensions of previous and current habitat
        if(all(dim(habitat_mask) == dim(prev_mask))==FALSE){
          stop("Dimension of previous habitat matrix does not equal dimension 
               of current one")
        }
        habitat_mask <- habitat_mask * prev_mask
      }
    }
  return(habitat_mask)
} # End function 'mask_raster'



#' Function to run spatial capture-recapture models in nimble using 
#' parallel processing
#'
#' @description A wrapper function to conduct Markov Chain Monte Carlo (MCMC) sampling using
#'  nimble.
#'
#' @param model The \code{nimbleCode} used to define model in \code{nimble} package,
#'  possibly generated from \code{\link{get_classic}}.
#' @param data A list of data inputs needed to run SCR models in \code{nimble}.
#' @param constants A list of constants needed to run SCR models in \code{nimble}.
#' @param inits Starting values for stochastic parameters to begin MCMC sampling.
#' @param params A vector of character strings that define the parameters to 
#' trace in the MCMC simulation.
#' @param niter An integer value of the total number of MCMC iterations to run 
#' per chain.
#' @param nburnin An integer value of the number of MCMC iterations to discard 
#' as 'burnin'.
#' @param thin An integer value of the amount of thinning of the chain. For 
#' example, \code{thin=2} would retain every other MCMC sample.
#' @param nchains An integer value for the number of MCMC chains to run
#' @param parallel A logical value indicating whether MCMC chains shoud be run 
#' in parallel processing. Default is \code{FALSE}.
#' @param RNGseed An integer value specifying the random number generating seed 
#' using in parallel processing. 
#' This ensures that the MCMC samples will be the same during each run using the
#'  same data, etc. Default is \code{NULL}.
#' @param s_alias A character value used to identify the latent activity center
#'  coordinates used in the model. Default is \code{"s"}. Note that the length
#'  of \code{s_alias} must be either 1 (e.g., \code{"s"}) or 2 
#'  (e.g., \code{c("s","su")}).
#' @return A list of MCMC samples for each parameter traced with length equal to
#'  the number of chains run.
#' @details This function provides a wrapper to easily run Bayesian SCR models
#'  using \code{nimble}.
#' @importFrom tictoc tic toc
#' @importFrom graphics hist par abline lines
#' @importFrom parallel parLapply makeCluster stopCluster clusterSetRNGStream
#' @import nimble
#' @author Daniel Eacker
#' @examples
#'\dontrun{
#' # simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' 
#' # add some random noise to locations
#' traps <- traps + runif(prod(dim(traps)),-20,20) 
#'
#' mysigma = 300 # simulate sigma of 300 m
#' mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#' pixelWidth = 100 # store pixelWidth
#'
# create grid and extent
#' Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, 
#' res = pixelWidth) # create slightly larger buffer area for example
#'
#' # create polygon to use as a mask
#' library(sf)
#' poly = st_sfc(st_polygon(x=list(matrix(c(-1765,-1765,1730,-1650,1600,1650,
#' 0,1350,-800,1700,-1850,1000,-1765,-1765),ncol=2, byrow=TRUE))), crs =  mycrs)
#'
#' # create habitat mask
#' hab_mask = mask_polygon(poly = poly, grid = Grid$grid, crs_ = mycrs, 
#' prev_mask = NULL)
#'
#' # simulate data for uniform state-space and habitat mask
#' data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, sigma_ = 
#' mysigma, prop_sex = 1, N = 200, K = 4, base_encounter = 0.15, 
#' enc_dist = "binomial",hab_mask = hab_mask, setSeed = 50)
#'
#' # get initial activity center starting values
#' s.st3d = initialize_classic(y=data3d$y, M=500, X=traps, buff = 3*mysigma, 
#' hab_mask = hab_mask)
#'
#' # rescale inputs
#' rescale_list = rescale_classic(X = traps, ext = Grid$ext, s.st = s.st3d, 
#' hab_mask = hab_mask)
#'
#' # store rescaled extent
#' ext = rescale_list$ext
#'
#' # prepare data
#' data = list(y=data3d$y)
#' # remove augmented records
#' data$y = data$y[which(apply(data$y, 1, sum)!=0),,] 
#' data$y = apply(data$y, c(1,2), sum) # covert to 2d by summing over occasions
#'
#' # add rescaled traps
#' data$X = rescale_list$X
#'
#' # prepare constants
#' constants = list(M = 500,n0 = nrow(data$y),J=dim(data$y)[2], 
#'                  K=dim(data3d$y)[3],x_lower = ext[1], x_upper = ext[2], 
#'                  y_lower = ext[3], y_upper = ext[4],sigma_upper = 1000, 
#'                  A = sum(hab_mask)*(pixelWidth)^2,pixelWidth=pixelWidth)
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
#' inits = list(sigma = runif(1, 250, 350), s = s.st3d,psi=runif(1,0.2,0.3),
#'           p0 = runif(1, 0.05, 0.15),   pOK = data$OK, 
#'           z=c(rep(NA,constants$n0),rep(0,constants$M-constants$n0)))
#'
#' # parameters to monitor
#' params = c("sigma","psi","p0","N","D")
#'
#' # get model
#' scr_model = get_classic(dim_y = 2, enc_dist = "binomial",sex_sigma = FALSE,
#' hab_mask=TRUE)
#'
#' # run model
#' library(tictoc)
#' tic() # track time elapsed
#' out = run_classic(model = scr_model, data=data, constants=constants, 
#'       inits=inits, params = params, niter = 1000, nburnin=500, 
#'       thin=1, nchains=2, parallel=TRUE, RNGseed = 500)
#' toc()
#'
#' # summarize output
#' samples = do.call(rbind, out)
#' par(mfrow=c(1,1))
#' hist(samples[,which(dimnames(out[[1]])[[2]]=="N")], xlab = "Abundance", 
#' xlim = c(0,500), main="")
#' abline(v=200, col="red") # add line for simulated abundance
#'
#' # not run
#' #nimSummary(out)
#'}
#' @name run_classic
#' @export
run_classic <- function(model, data, constants, inits, params,
                        niter = 1000, nburnin=100, thin=1, nchains=1, 
                        parallel=FALSE, RNGseed=NULL,s_alias="s"){
  # for parallel processing
  if(parallel == FALSE){
    SCRmodelR <- nimble::nimbleModel(code=model,data=data,constants=constants,
                            inits=inits,check=FALSE,calculate=TRUE)
    SCRmodelR$initializeInfo()
    # compile model to C++#
    SCRmodelC <- nimble::compileNimble(SCRmodelR) # compile code
    # MCMC configurations
    mcmcspec<-nimble::configureMCMC(SCRmodelR, monitors=params) 
    # block updating for individual activity centers
    if(length(s_alias)!=1 & length(s_alias)!=2){
    stop("s_alias must be either a scalar or vector of length 2")
    }  
    if(length(s_alias)==1 && s_alias == "s"){
    mcmcspec$removeSamplers(s_alias, print = FALSE)
    for(i in 1:constants$M){
      snew = paste(s_alias,"[",i,","," 1:2","]",sep="")
      mcmcspec$addSampler(target = snew, type = 'RW_block', silent = TRUE)
    }
    }else
    if(length(s_alias)==1 && s_alias == "su"){ 
    # test out effect of AF_slice sampler to reduce autocorrlation and cross correlation
    for(i in 1:constants$m){
      snew = paste(s_alias,"[",i,","," 1:2","]",sep="")
      mcmcspec$addSampler(target = snew, type = 'RW_block', silent = TRUE)
    }
     mcmcspec$removeSamplers(c("psiu","lam0","sigma"))
     mcmcspec$addSampler(target = c("psiu","lam0","sigma"), type = "AF_slice", silent = TRUE)
    
    }else
    if(length(s_alias)==2){ # for  spatial mark-resight model
    mcmcspec$removeSamplers(s_alias[which(s_alias=="s")], print = FALSE)
    for(i in 1:constants$M){
      snew = paste(s_alias[which(s_alias=="s")],"[",i,","," 1:2","]",sep="")
      mcmcspec$addSampler(target = snew, type = 'RW_block', silent = TRUE)
    }
    mcmcspec$removeSamplers(s_alias[which(s_alias=="su")], print = FALSE)
    for(i in 1:constants$m){
      snew = paste(s_alias[which(s_alias=="su")],"[",i,","," 1:2","]",sep="")
      mcmcspec$addSampler(target = snew, type = 'RW_block', silent = TRUE)
    }
    }
    scrMCMC <- nimble::buildMCMC(mcmcspec)
    CscrMCMC <- nimble::compileNimble(scrMCMC, project = SCRmodelR,
                                      resetFunctions = TRUE)
    results <- nimble::runMCMC(CscrMCMC, niter = niter, nburnin=nburnin,
                               thin=thin, nchains=nchains)
    return(results)
  }else
    if(parallel == TRUE){
      run_parallel <- function(seed,data,model,constants,inits, params,
                               niter,nburnin,thin){
    SCRmodelR <- nimble::nimbleModel(code=model,data=data,
            constants=constants,inits=inits,check=FALSE,calculate=TRUE)
    SCRmodelR$initializeInfo()
    # compile model to C++#
    SCRmodelC <- nimble::compileNimble(SCRmodelR) # compile code
    # MCMC configurations
    mcmcspec<-nimble::configureMCMC(SCRmodelR, monitors=params) 
    # block updating for individual activity centers
    if(length(s_alias)!=1 & length(s_alias)!=2){
    stop("s_alias must be either a scalar or vector of length 2")
    }  
    if(length(s_alias)==1 && s_alias == "s"){
    mcmcspec$removeSamplers(s_alias, print = FALSE)
    for(i in 1:constants$M){
      snew = paste(s_alias,"[",i,","," 1:2","]",sep="")
      mcmcspec$addSampler(target = snew, type = 'RW_block', silent = TRUE)
    }
    }else
    if(length(s_alias)==1 && s_alias == "su"){ 
    # test out effect of AF_slice sampler to reduce autocorrlation and cross correlation
    for(i in 1:constants$m){
      snew = paste(s_alias,"[",i,","," 1:2","]",sep="")
      mcmcspec$addSampler(target = snew, type = 'RW_block', silent = TRUE)
    }
     mcmcspec$removeSamplers(c("psiu","lam0","sigma"))
     mcmcspec$addSampler(target = c("psiu","lam0","sigma"), type = "AF_slice", silent = TRUE)
    
    }else
    if(length(s_alias)==2){ # for  spatial mark-resight model
    mcmcspec$removeSamplers(s_alias[which(s_alias=="s")], print = FALSE)
    for(i in 1:constants$M){
      snew = paste(s_alias[which(s_alias=="s")],"[",i,","," 1:2","]",sep="")
      mcmcspec$addSampler(target = snew, type = 'RW_block', silent = TRUE)
    }
    mcmcspec$removeSamplers(s_alias[which(s_alias=="su")], print = FALSE)
    for(i in 1:constants$m){
      snew = paste(s_alias[which(s_alias=="su")],"[",i,","," 1:2","]",sep="")
      mcmcspec$addSampler(target = snew, type = 'RW_block', silent = TRUE)
    }
    }        
        scrMCMC <- nimble::buildMCMC(mcmcspec)
        CscrMCMC <- nimble::compileNimble(scrMCMC, project = SCRmodelR,
                                          resetFunctions = TRUE)
        results <- nimble::runMCMC(CscrMCMC, niter = niter, nburnin= nburnin,
                                   thin=thin,setSeed = seed)
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
      chain_output <- parallel::parLapply(cl = this_cluster, fun = run_parallel,
                      X=1:nchains,model = model, data=data,constants=constants, 
                      inits=inits,params = params, niter = niter, 
                      nburnin=nburnin, thin=thin)
      parallel::stopCluster(this_cluster)
      return(chain_output)
    }
} # End function 'run_classic'




#' Function to summarize MCMC chain output from nimble
#'
#' @description Summarizes lists of MCMC chain output from nimble
#'
#' @param d A list of MCMC samples from each chain returned from 
#' \code{\link{run_classic}}.
#' @param trace A logical value indicating whether or not traces of MCMC samples
#'  should be plotted. Default is \code{FALSE}.
#' @param plot_all A logical value indicating whether or not all parameters 
#' should be included in plots. This assumes that some parameters
#' are excluded in the summary table (i.e., \code{exclude_params != NULL}). 
#' Default is \code{FALSE}.
#' @param exclude_params Either \code{NULL} or a scalar or vector containing 
#' parameter(s) to exclude from summary. Note that high dimensional parameters 
#' (e.g., \code{s[1, 1, 1]}) can be excluded by calling \code{exclude_params =  
#' "s"}. Default is \code{NULL}.
#' @param digits An integer value indicating how many digits the output should 
#' be rounded to.
#' @return A dataframe of summary statistics for MCMC samples.
#' @details This function summarizes Bayesian SCR models run using \code{nimble} 
#' including mean and quantiles of samples, as well as effective sample size and
#'  Rhat statistics. Note that \code{f0} is the proportion of samples that are
#'   greater than zero. Also, at least 2 chains must be run to use this function.
#' @author Daniel Eacker
#' @importFrom stringr str_extract str_locate
#' @importFrom coda effectiveSize mcmc.list as.mcmc gelman.diag
#' @importFrom stats na.omit sd quantile
#' @importFrom crayon red
#' @importFrom viridisLite viridis
#' @examples
#'\dontrun{
#' # simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' 
#' # add some random noise to locations
#' traps <- traps + runif(prod(dim(traps)),-20,20) 
#'
#' mysigma = 300 # simulate sigma of 300 m
#' mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#' pixelWidth = 100 # store pixelWidth
#'
#' # create grid and extent
#' Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, 
#' res = pixelWidth)
#'
#' # simulate encounter data
#' data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, sigma = mysigma, 
#' prop_sex = 1, N = 200, K = 4, base_encounter = 0.15, enc_dist = "binomial", 
#' hab_mask = FALSE, setSeed = 200)
#'
#' # prepare data
#' data = list(y=data3d$y)
#' data$y = data$y[which(apply(data$y, 1, sum)!=0),,] # remove augmented records
#' data$y = apply(data$y, c(1,2), sum) # covert to 2d by summing over occasions
#' data$X = traps/1000 # rescale to kilometers
#' ext = as.vector(Grid$ext)/1000 # recale to kilometers
#'
#' # prepare constants
#' constants = list(M = 500,n0 = nrow(data$y),J=dim(data$y)[2], 
#'             K=dim(data3d$y)[3],x_lower = ext[1], x_upper = ext[2], 
#'            y_lower = ext[3], y_upper = ext[4],sigma_upper = 1, 
#'            A = prod(ext[2]-ext[1],ext[4]-ext[3]))
#'
#' # add z and zeros vector data for latent inclusion indicator
#' data$z = c(rep(1,constants$n0),rep(NA,constants$M - constants$n0))
#' data$zeros =  c(rep(NA,constants$n0),rep(0,constants$M - constants$n0))
#'
#' # get initial activity center starting values
#' s.st3d = initialize_classic(y=data3d$y, M=500, X=traps, buff = 3*mysigma, 
#' hab_mask=FALSE)
#'
#' # define all initial values
#' inits = list(sigma = runif(1, 0.250, 0.350), s = s.st3d/1000,
#'              psi=runif(1,0.2,0.3), p0 = runif(1, 0.05, 0.15),
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
#' out = run_classic(model = scr_model, data=data, constants=constants, 
#'                   inits=inits, params = params,niter = 5000, nburnin=1000,
#'                   thin=1, nchains=2, parallel=TRUE, RNGseed = 500)
#' toc()
#'
#' # summarize output
#' samples = do.call(rbind, out)
#' par(mfrow=c(1,1))
#' hist(samples[,which(dimnames(out[[1]])[[2]]=="N")], xlab = "Abundance", 
#' xlim = c(0,500), main="")
#' abline(v=200, col="red") # add line for simulated abundance
#'
#' # summarize output
#' nimSummary(out, trace=TRUE, plot_all=TRUE)
#'}
#' @name nimSummary
#' @export
nimSummary <- function(d, trace=FALSE, plot_all=FALSE, exclude_params = NULL, 
                       digits=3){
  if(is.null(exclude_params)==FALSE){
    tmp1 <- ifelse(is.na(as.numeric(stringr::str_extract(
      attributes(d[[1]])$dimnames[[2]],"[1-9]+"))),
      attributes(d[[1]])$dimnames[[2]],
      substr(attributes(d[[1]])$dimnames[[2]],1,
      as.numeric(stringr::str_locate(attributes(d[[1]])$dimnames[[2]],
                                     "\\[")[, 1])-1))
    d.remove <- lapply(d, function(x) which(tmp1 %in% exclude_params))
    d2 <- lapply(d, function(x) x[,-d.remove[[1]]])
  }else
if(is.null(exclude_params)){
  if(dim(d[[1]])[2]>100){
   param_length_check <- readline(cat(crayon::red(paste(
   "Do you really want to 
    build a summary table for",dim(d[[1]])[2],"parameters? (1 = Yes or 0 = No); 
    note that this will take forever! \n",
    "Recommend selecting 0 and setting exclude_params"))))
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
    Rhat <- round(coda::gelman.diag(coda::mcmc.list(lapply(d2, coda::as.mcmc)), 
    multivariate = FALSE)[[1]][,1],3)
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
      Rhat <- round(coda::gelman.diag(coda::mcmc.list(lapply(d2, coda::as.mcmc)),
      multivariate = FALSE)[[1]][,1],3)
    }
 if(trace){
  if((length(attributes(d[[1]])$dimnames[[2]])-length(d.remove[[1]])>1)){
     par(mfrow=c(3,2))
     g <- 1 # set g index
     if(plot_all){
     d2 <- d
     d3 <- do.call(rbind, d2)
     }
     plot.seq <- seq(3,3000,3)  
     cols <- rev(viridisLite::viridis(100)[seq(1,100,length.out=length(d2))])
     for(i in 1:3){ # plot first 3 variables to start out
      plot(1:dim(d2[[1]])[1],d2[[1]][,i],col=cols[1],xlab="iteration",
           ylab=colnames(d3)[i],type="l",ylim=range(do.call(rbind, 
           lapply(d2,function(x) apply(x, 2, range)))[,i]))
      for(j in 2:length(d2)){
        lines(1:dim(d2[[1]])[1],d2[[j]][,i],xlab="iteration",
            ylab=colnames(d3)[i],type="l",col=cols[j])
      }
      hist(d3[,i],main="",xlab=colnames(d3)[i])
     }
      if(interactive() & ncol(d3) > 3){
        answer <- readline(cat(crayon::red("Plot next set of parameters?\n 
 (1 = Yes or 0 = No) ")))
        while(answer == "1"){
          upper_index <- ifelse(plot.seq[g+1] > ncol(d3), ncol(d3), 
                                plot.seq[g+1])
          for(i in (plot.seq[g]+1):upper_index){
             plot(1:dim(d2[[1]])[1],d2[[1]][,i],col=cols[1],xlab="iteration",
            ylab=colnames(d3)[i],type="l",ylim=range(do.call(rbind, 
            lapply(d2,function(x) apply(x, 2, range)))[,i]))
              for(j in 2:length(d2)){
              lines(1:dim(d2[[1]])[1],d2[[j]][,i],xlab="iteration",
            ylab=colnames(d3)[i],type="l",col=cols[j])
             }
             hist(d3[,i],main="",xlab=colnames(d3)[i])
            } # end i loop
          g <- g + 1 
          if(plot.seq[g+1] <= ncol(d3)){
        answer <- readline(cat(crayon::red("Plot next set of parameters?\n 
 (1 = Yes or 0 = No) ")))
          }else
          if(plot.seq[g+1] > ncol(d3)){
          break
          }
        }
      } # is interactive
  }else
    if((length(attributes(d[[1]])$dimnames[[2]])-length(d.remove[[1]])==1)){
      cols <- rev(viridisLite::viridis(100)[seq(1,100,length.out=length(d2))])
      par(mfrow=c(1,2))
      plot(1:length(d2[[1]]),d2[[1]],xlab="iteration",ylab=colnames(d3)[i],
           type="l",ylim=range(d3[,1]),col=cols[1])
      for(j in 2:length(d2)){
        lines(1:length(d2[[j]]),d2[[j]],xlab="iteration",
              ylab=colnames(d3)[i],type="l",col=cols[j])
      }
      hist(d3[,1],main="",xlab=attributes(d[[1]])$dimnames[[2]][-d.remove[[1]]])
    }
     on.exit(par(mfrow=c(1,1))) # reset graphical frame
  } # end if trace
  tmp.frame = data.frame(post.mean=Means,post.sd=SDs,q2.5=q2.5,q50=q50,
                         q97.5=q97.5,f0=over.zero,n.eff=n.eff,Rhat=Rhat)
  if(nrow(tmp.frame)==1){
    row.names(tmp.frame) = attributes(d[[1]])$dimnames[[2]][-d.remove[[1]]]
  }
  return(round(tmp.frame, digits=digits))
} # End function 'nimSummary'



#' Function to generate realized density surface from MCMC output
#'
#' @description Streamlined construction of realized density surface from posterior
#' samples of latent indicator variable (\code{z}) and activity center
#' coordinates (\code{s})
#'
#' @param samples Either a matrix (for a single MCMC chain) or list of posterior
#' samples from multiple chains from MCMC sampling; possibly returned from 
#' \code{\link{run_classic}}.
#' @param grid A matrix or array object of the the state-space grid. This is 
#' returned from \code{\link{grid_classic}}.
#' @param crs_ The UTM coordinate reference system (EPSG code) used for your
#' location provided as an integer (e.g., 32608 for WGS 84/UTM Zone 8N).
#' @param site Site identifier variable included for detected and augmented 
#' individuals used as a constant in model runs. 
#' @param hab_mask Either \code{FALSE} (the default) or a matrix or arary 
#' output from \code{\link{mask_polygon}}
#' or \code{\link{mask_raster}} functions.
#' @param s_alias A character value used to identify the latent activity center 
#' coordinates used in the model. Default is \code{"s"}.
#' @param z_alias A character value used to identify the latent inclusion 
#' indicator used in the model. Default is \code{"z"}.
#' @importFrom sp SpatialPoints
#' @importFrom raster cellFromXY
#' @author Daniel Eacker
#' @return a \code{raster} object or a list of \code{raster} objects if 
#' state-space grid is an array. 
#' @details This function automates the construction of realized density 
#' surfaces from MCMC samples.
#' @examples
#' \dontrun{
#' # simulate a single trap array with random positional noise
#' x <- seq(-800, 800, length.out = 5)
#' y <- seq(-800, 800, length.out = 5)
#' traps <- as.matrix(expand.grid(x = x, y = y))
#' 
#' # add some random noise to locations
#' traps <- traps + runif(prod(dim(traps)),-20,20) 
#' 
#' mysigma = 300 # simulate sigma of 300 m
#' mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#' pixelWidth = 100 # store pixelWidth
#' 
#' # create grid and extent
#' Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, 
#' res = pixelWidth)
#' 
#' # simulate encounter data
#' data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, sigma = mysigma,
#'  prop_sex = 1, N = 200, K = 4, base_encounter = 0.15, enc_dist = "binomial",
#'  hab_mask = FALSE, setSeed = 200)
#' 
#' # prepare data
#' data = list(y=data3d$y)
#' data$y = data$y[which(apply(data$y, 1, sum)!=0),,] # remove augmented records
#' data$y = apply(data$y, c(1,2), sum) # covert to 2d by summing over occasions
#' data$X = traps # rescale to kilometers
#' ext = as.vector(Grid$ext) # recale to kilometers
#' 
#' # prepare constants
#' constants = list(M = 500,n0 = nrow(data$y),J=dim(data$y)[2], 
#'             K=dim(data3d$y)[3],x_lower = ext[1], x_upper = ext[2], 
#'             y_lower = ext[3], y_upper = ext[4],sigma_upper = 1000, 
#'             A = prod(ext[2]-ext[1],ext[4]-ext[3]))
#' 
#' # add z and zeros vector data for latent inclusion indicator
#' data$z = c(rep(1,constants$n0),rep(NA,constants$M - constants$n0))
#' data$zeros =  c(rep(NA,constants$n0),rep(0,constants$M - constants$n0))
#' 
#' # get initial activity center starting values
#' s.st3d = initialize_classic(y=data3d$y, M=500, X=traps, buff = 3*mysigma, 
#' hab_mask=FALSE)
#' 
#' # define all initial values
#' inits = list(sigma = runif(1, 250, 350), s = s.st3d,psi=runif(1,0.2,0.3), 
#'         p0 = runif(1, 0.05, 0.15),
#'          z=c(rep(NA,constants$n0),rep(0,constants$M-constants$n0)))
#' 
#' # parameters to monitor
#' params = c("sigma","psi","p0","N","D","s","z")
#' 
#' # get model
#' scr_model = get_classic(dim_y = 2, enc_dist = "binomial",sex_sigma = FALSE,
#' trapsClustered = FALSE)
#' 
#' # run model
#' library(tictoc())
#' tic() # track time elapsed
#' out = run_classic(model = scr_model, data=data, constants=constants, 
#'          inits=inits, params = params, niter = 5000, nburnin=1000, 
#'          thin=1, nchains=2, parallel=TRUE, RNGseed = 500)
#' toc()
#' 
#' # summarize output (exclude lengthy parameters "s" and "z")
#' nimSummary(out, exclude_params = c("s","z"), trace = TRUE)
#' 
#' library(tictoc)       
#' tic() # track time
#' r = realized_density(samples = out, grid = Grid$grid, crs_ = mycrs, 
#' site = NULL, hab_mask = FALSE)       
#' toc()      
#' 
#' # load virdiis color pallete library      
#' library(viridis)
#' library(raster)
#'
#' # make simple raster plot
#' label = expression("Realized density (activity centers/100 m"^2*")")
#' plot(r, col=viridis(100),main=label)
#'}
#' @name realized_density
#' @export
realized_density <- function(samples, grid, crs_, site, hab_mask, 
                             s_alias = "s", z_alias = "z"){
   if(length(samples) > 1){
    samples <- do.call(rbind, samples) # rbind mcmc samples
   }
   n.iter <- dim(samples)[1] # store total number of iterations
   # need to determine if state-space grid is 2 or 3 dimensional (stop if not)
   if(length(dim(grid))!=2 & length(dim(grid))!=3){
     stop("State-space grid must be only 2 or 3 dimensions")
   } 
   if(length(dim(grid))==2){
# process z indicator variable and sx and sy activity center coordinates as 
     # vectors for speed
       z <- samples[,grep(paste0(z_alias,"\\["), colnames(samples))]
       z_vec <- as.vector(z)
       sx <- samples[,grep(paste0(s_alias,"\\["), 
             colnames(samples))][,1:dim(z)[2]] # x-coordinate
       sx_vec <- as.vector(sx)[z_vec==1]
       sy <- samples[,grep(paste0(s_alias,"\\["), 
            colnames(samples))][,(dim(z)[2]+1):(dim(z)[2]*2)] # y-coordinates
       sy_vec <- as.vector(sy)[z_vec==1]
     if(isFALSE(hab_mask)){
         ac_pts <- sp::SpatialPoints(cbind(sx_vec, sy_vec),
                  proj4string = sp::CRS(sf::st_crs(crs_)$proj4string))
         r <- raster::rasterFromXYZ(grid,crs = crs_)
         tab <- table(raster::cellFromXY(r, ac_pts))
         r[as.numeric(names(tab))] <- tab/n.iter
     } else 
     if(isFALSE(hab_mask)==FALSE){
         r <- raster::rasterFromXYZ(grid,crs = crs_)
         rescale.sx <- scales::rescale(sx_vec, to = raster::extent(r)[1:2],
                                       from=c(0,dim(hab_mask)[2]))
         rescale.sy <- scales::rescale(sy_vec, to = raster::extent(r)[3:4], 
                                       from=c(0,dim(hab_mask)[1]))
         ac_pts <- sp::SpatialPoints(cbind(rescale.sx, rescale.sy),
                        proj4string = sp::CRS(sf::st_crs(crs_)$proj4string))
         tab <- table(raster::cellFromXY(r, ac_pts))
         r[as.numeric(names(tab))] <- tab/n.iter
     } # end habitat mask
   } else # end 2D grid 
 if(length(dim(grid))==3){
   if(isFALSE(hab_mask)){
       r <- apply(grid,3,function(x) raster::rasterFromXYZ(x,crs = crs_))
 # process z indicator variable and sx and sy activity center coordinates 
       # as vectors for speed
         z <- samples[,grep(paste0(z_alias,"\\["), 
               colnames(samples))] # indicator variable
         sx <- samples[,grep(paste0(s_alias,"\\["), 
               colnames(samples))][,1:dim(z)[2]] # x-coordinate
         sy <- samples[,grep(paste0(s_alias,"\\["), 
              colnames(samples))][,(dim(z)[2]+1):(dim(z)[2]*2)] # y-coordinates
   for(g in 1:dim(grid)[3]){
     z_site <- z[,site==g]
     z_vec <- as.vector(z_site)
     sx_site <- sx[,site==g]
     sx_vec <- as.vector(sx_site)[z_vec==1]
     sy_site <- sy[,site==g]             
     sy_vec <- as.vector(sy_site)[z_vec==1]
     ac_pts <- sp::SpatialPoints(cbind(sx_vec, sy_vec),
                    proj4string = sp::CRS(sf::st_crs(crs_)$proj4string))
     tab <- table(raster::cellFromXY(r[[g]], ac_pts))
     r[[g]][as.numeric(names(tab))] <- tab/n.iter 
   } # end g loop   
       } else
 if(isFALSE(hab_mask)==FALSE){  
     r <- apply(grid,3,function(x) raster::rasterFromXYZ(x,crs = crs_))
     for(g in 1:dim(grid)[3]){
       # process z indicator variable and sx and sy activity center 
       # coordinates as vectors for speed
       z <- samples[,grep(paste0(z_alias,"\\["), 
                  colnames(samples))] # indicator variable
       z_site <- z[,site==g]
       z_vec <- as.vector(z_site)
       sx <- samples[,grep(paste0(s_alias,"\\["),
                  colnames(samples))][,1:dim(z)[2]] # x-coordinate
       sx_site <- sx[,site==g]
       sx_vec <- as.vector(sx_site)[z_vec==1]
       sy <- samples[,grep(paste0(s_alias,"\\["), 
               colnames(samples))][,(dim(z)[2]+1):(dim(z)[2]*2)] # y-coordinates
       sy_site <- sy[,site==g]             
       sy_vec <- as.vector(sy_site)[z_vec==1]
       rescale.sx <- scales::rescale(sx_vec, 
                  to = raster::extent(r[[g]])[1:2], from=c(0,dim(hab_mask)[2]))
       rescale.sy <- scales::rescale(sy_vec, to = raster::extent(r[[g]])[3:4],
                  from=c(0,dim(hab_mask)[1]))
       ac_pts <- sp::SpatialPoints(cbind(rescale.sx, rescale.sy),
                  proj4string = sp::CRS(sf::st_crs(crs_)$proj4string))
       tab <- table(raster::cellFromXY(r[[g]], ac_pts))
       r[[g]][as.numeric(names(tab))] <- tab/n.iter 
     } # end g loop       
 } # end habitat mask    
}
     return(r)
} # End function 'realized_density'




#' Function to efficiently edit rows of model code generated from nimble
#'
#' @description Allows for efficient editing of model code produced by \code{nimbleCode()} 
#' function
#'
#' @param model The \code{nimbleCode()} used to define model in \code{nimble} package,
#' possibly generated from \code{\link{get_classic}} function.
#' @param append_code Either \code{NULL} or model code produced from 
#' \code{nimbleCode()}  or \code{\link{get_classic}} function. Note that if
#' \code{remove_line = NULL}, then code will be appended just after existing model
#' code; otherwise specify the lines to insert new code into by setting 
#' \code{append_line}.
#' @param append_line Either \code{NULL} or an integer value as a scalar or vector
#' defining which positions to insert new lines of code in the model. Note that if multiple 
#' lines of new code are to be inserted on the same line, then just use \code{rep(44,3)} for 
#' example if the new code had 3 lines to insert (not counting \code{"{"} and \code{"}"}). The lines
#' to append to should be based on the original model given to \code{model}.
#' @param remove_line Either \code{NULL} or an integer value as a scalar or vector
#' defining which lines of code to remove from the original model. Set to \code{NULL} when only 
#' appending to and not replacing code in previous model file (i.e., \code{model}).
#' @param write Logical. If \code{TRUE}, then a text file is written to the 
#' working directory called "new_model.txt". Otherwise, model is written to temp
#' file and then deleted. Default is \code{FALSE}.
#' @return A model description that can be run  using \code{\link{run_classic}}. 
#' @author Daniel Eacker
#' @importFrom R.utils insert
#' @examples
#' # get model
#' scr_model = get_classic(dim_y = 2, enc_dist = "binomial",sex_sigma = TRUE,
#' hab_mask=TRUE,trapsClustered = TRUE)
#' 
#' # create new nimbleCode to use for replacement in 'scr_model'
#' p0_prior = nimble::nimbleCode({
#'    p0[g] ~ dbeta(1,1)
#' })
#' 
#' # replace line 3 of old model code with 'p0_prior' 
#' new_model = customize_model(model = scr_model, append_code = p0_prior, 
#'                                    append_line = 3, remove_line = 3)
#'                                    
#' # inspect new model code
#' new_model
#' @name customize_model
#' @export 
customize_model <- function(model,append_code = NULL,append_line = NULL, 
                        remove_line = NULL,write = FALSE){
  # read in main model
     model_list <- as.list(model)
     file_list <- list()
     length_lines <- numeric(length(model_list))
   for(i in 1:length(model_list)){
     tmodel <- strsplit(as.character(model[i]),"\n")[[1]]
     file_list[[i]] <- tempfile(fileext = ".txt")
     writeLines(tmodel, file_list[[i]])
     length_lines[i] <- length(readLines(file_list[[i]],encoding="UTF-8"))
   }
   main_model <- character(sum(length_lines))
   for(i in 1:length(model_list)){
    main_model[which(main_model=="")[1]:ifelse(length_lines[i]==1,(which(main_model=="")[1]),
                  (which(main_model=="")[1])+length_lines[i]-1)] =
                  readLines(file_list[[i]],encoding="UTF-8")
   }
   main_model[length(main_model)+1] <- "}"
   main_model[c(1,length(main_model))] <-  c("nimble::nimbleCode({","})") 
   unlink(unlist(file_list)) # unlink temp files
  # check for removal of lines and do this first
  if(isFALSE(is.null(remove_line))){
      remove_model <- main_model[-remove_line]
  }
  # just append new code after old code
  if(isFALSE(is.null(append_code)) & is.null(append_line)){ 
       model_list <- as.list(append_code)
       file_list <- list()
       length_lines <- numeric(length(model_list))
   for(i in 1:length(model_list)){
     tmodel <- strsplit(as.character(append_code[i]),"\n")[[1]]
     file_list[[i]] <- tempfile(fileext = ".txt")
     writeLines(tmodel, file_list[[i]])
     length_lines[i] <- length(readLines(file_list[[i]],encoding="UTF-8"))
   }
   append_model <- character(sum(length_lines))
   for(i in 1:length(model_list)){
    append_model[which(append_model=="")[1]:ifelse(length_lines[i]==1,(which(append_model=="")[1]),
                  (which(append_model=="")[1])+length_lines[i]-1)] =
                  readLines(file_list[[i]],encoding="UTF-8")
   }
   append_model[length(append_model)+1] <- "}"
   unlink(unlist(file_list))
   if(is.null(remove_line)){
    updated_model <- c("nimble::nimbleCode({",
          main_model[-c(1,length(main_model))],
          append_model[-c(1,length(append_model))],
                  "})") 
   }else
   if(isFALSE(is.null(remove_line))){ 
   updated_model <- c("nimble::nimbleCode({",
          remove_model[-c(1,length(remove_model))],
          append_model[-c(1,length(append_model))],
                  "})")
   } # check for line removal
  }else # add new code to specific index position in character vector
  if(isFALSE(is.null(append_code)) & isFALSE(is.null(append_line))){
       model_list <- as.list(append_code)
       file_list <- list()
       length_lines <- numeric(length(model_list))
   for(i in 1:length(model_list)){
     tmodel <- strsplit(as.character(append_code[i]),"\n")[[1]]
     file_list[[i]] <- tempfile(fileext = ".txt")
     writeLines(tmodel, file_list[[i]])
     length_lines[i] <- length(readLines(file_list[[i]],encoding="UTF-8"))
   }
   append_model <- character(sum(length_lines))
   for(i in 1:length(model_list)){
    append_model[which(append_model=="")[1]:ifelse(length_lines[i]==1,(which(append_model=="")[1]),
                  (which(append_model=="")[1])+length_lines[i]-1)] =
                  readLines(file_list[[i]],encoding="UTF-8")
   }
   append_model[length(append_model)+1] <- "}"
   append_model[c(1,length(append_model))] <-  c("nimble::nimbleCode({","})") 
   unlink(unlist(file_list))
   # check if length(append_line) == length(append_model[-c(1,length(append_model))])
   if(length(append_line) != length(append_model[-c(1,length(append_model))])){
     stop(paste("length of append_line must equal lines of append_code minus outside braces, i.e.,","{","}"))
   }
    if(is.null(remove_line)){
        updated_model <- character(length(main_model)+length(append_model)-2)
        updated_model[append_line] <- append_model[-c(1,length(append_model))]
        updated_model[(1:(length(updated_model)-1))[-c(1,append_line)]] <- main_model[-c(1,length(main_model))]
        updated_model[c(1,length(updated_model))] <- c("nimble::nimbleCode({","})") 
   }else
   if(isFALSE(is.null(remove_line))){ 
        updated_model <- character(length(main_model))
        updated_model[(1:length(main_model))[-remove_line]] <- remove_model
        if(min(append_line) <= length(main_model)){
        updated_model <- R.utils::insert(updated_model, ats = append_line, 
                                        append_model[-c(1,length(append_model))])
        updated_model <- updated_model[which(updated_model!="")]
        }else
        if(min(append_line) > length(main_model)){
        updated_model <- c(updated_model[-length(updated_model)],
                      rep("",min(append_line) - length(main_model)),"})")
        updated_model <- R.utils::insert(updated_model, ats = append_line, 
                                        append_model[-c(1,length(append_model))])
        updated_model <- updated_model[which(updated_model!="")]  
        }
   } # check for line removal
  } else # end append_code and append_line = TRUE
    # if nothing is done to code, just return original code
  if(is.null(remove_line) & is.null(append_code) & is.null(append_line)){
    updated_model <- main_model
     warning("Returning same model code as input into function")
  }else
   # if only removal of lines
  if(isFALSE(is.null(remove_line)) & is.null(append_code) & is.null(append_line)){
    updated_model <- remove_model
  }
  # if write to file
  if(write){
    local.path <- paste0(getwd(),"/new_model.txt")
    writeLines(updated_model,local.path,useBytes = FALSE) 
  }
  txtPath1 <- tempfile(fileext = ".txt")
  writeLines(updated_model,txtPath1,useBytes = FALSE) 
  return(source(txtPath1)$value)
 # on.exit(unlink(txtPath1))
} # End function 'customize_model'



#' Function to retrieve nimbleCode for spatial count models
#'
#' @description Creates model code using the \code{\link[nimble]{nimbleCode}} function.
#'
#' @param occ_specific Logical. If \code{FALSE}, the encounter rate will
#' not include an occasion-specific loop in the detection function; otherwise, 
#' the model will include a for loop for occasions (K) in the detection function.
#' Default is \code{FALSE}.
#' @param hab_mask A logical value indicating whether a habitat mask will be 
#' used. Default is \code{FALSE}.
#' @param trapsClustered A logical value indicating if traps are clustered in 
#' arrays across the sampling area.
#' @return Model code created from \code{\link[nimble]{nimbleCode}}.
#' @details This function provides templates for unmarked models that can be 
#' easily modified to include further model complexity such as covariates 
#' explaining detection probability. The models include habitat masking.
#' @author Daniel Eacker
#' @examples
#' # get spatial count model with non-occasion-specific detection
#' # function, single scaling parameter, no habitat mask, and no clustering
#' unmarked_model = get_unmarked(occ_specific=FALSE,hab_mask = FALSE, 
#'                               trapsClustered = FALSE)
#'
#' # inspect model
#' unmarked_model
#' @name get_unmarked
#' @export
get_unmarked <- function(occ_specific = FALSE,
                        hab_mask = FALSE,trapsClustered = FALSE){
      m <- J <- su <- X <- sigma <- n0 <- zu <- A <- lam0 <- K <- 
      nSites <-  site <- pixelWidth <- psiu <- prop.habitat <- site_indexL <-
      site_indexU <- NULL
 if(trapsClustered == FALSE){    
  if(isFALSE(hab_mask)){ # determine if hab_mask is included
    if(isFALSE(occ_specific)){
          scrcode <- nimble::nimbleCode({
            lam0 ~ dunif(0,lam0_upper) # baseline encounter probability
            sigma ~ dunif(0, sigma_upper) # scaling parameter
            psiu ~ dunif(0, 1) # inclusion prob
            for(i in 1:m){ 
              zu[i]~dbern(psiu)
              su[i,1] ~ dunif(x_lower, x_upper)
              su[i,2] ~ dunif(y_lower, y_upper)
              distu[i,1:J] <- sqrt((su[i,1]-X[1:J,1])^2 + (su[i,2]-X[1:J,2])^2)
              lamu[i,1:J] <- lam0*exp(-distu[i,1:J]^2/(2*sigma^2))*zu[i]
            } # i individuals
            # unmarked spatial count model likelihood
            for(j in 1:J){
              bigLambda[j] <- sum(lamu[1:m,j])
            for(k in 1:K){
              n[j,k] ~ dpois(bigLambda[j])
            } # k occasions
           } # j traps
            N <- sum(zu[1:m])
            D <- N/A
          })
      return(scrcode)
    } else    # End occ_specific = FALSE
      if(isFALSE(occ_specific)==FALSE){
            scrcode <- nimble::nimbleCode({
              lam0 ~ dunif(0,lam0_upper) # baseline encounter rate
              sigma ~ dunif(0, sigma_upper) # scaling parameter
              psiu ~ dunif(0, 1) # inclusion prob
              for(i in 1:m){
                zu[i]~dbern(psiu)
                su[i,1] ~ dunif(x_lower, x_upper)
                su[i,2] ~ dunif(y_lower, y_upper)
                distu[i,1:J] <- sqrt((su[i,1]-X[1:J,1])^2 + (su[i,2]-X[1:J,2])^2)
                for(k in 1:K){
                  lamu[i,1:J,k] <- lam0*exp(-distu[i,1:J]^2/(2*sigma^2))*zu[i]
                } # k occasions
              } # i individuals
              # unmarked spatial count model likelihood
              for(j in 1:J){
                bigLambda[j] <- sum(lamu[1:m,j,1:K]) 
              for(k in 1:K){
                n[j,k] ~ dpois(bigLambda[j])
               } # k occasions
              } # j traps
              N <- sum(zu[1:m])
              D <- N/A
            })
        return(scrcode)
      } 
  } else # End no habitat mask
if(hab_mask==TRUE){
 if(isFALSE(occ_specific)){
     scrcode <- nimble::nimbleCode({
      lam0 ~ dunif(0,lam0_upper) # baseline encounter probability
      sigma ~ dunif(0, sigma_upper) # scaling parameter
      sigma.pixel <- sigma / pixelWidth # scaled for habitat mask
      psiu ~ dunif(0, 1) # inclusion prob
    for(i in 1:m){
      zu[i]~dbern(psiu)
      su[i,1] ~ dunif(x_lower, x_upper)
      su[i,2] ~ dunif(y_lower, y_upper)
      pOKu[i] <- hab_mask[(trunc(su[i,2])+1),(trunc(su[i,1])+1)] # habitat check
      OKu[i] ~ dbern(pOKu[i]) # OKu[i] <- 1, the ones trick
      distu[i,1:J] <- sqrt((su[i,1]-X[1:J,1])^2 + (su[i,2]-X[1:J,2])^2)
      lamu[i,1:J] <- lam0*exp(-distu[i,1:J]^2/(2*sigma.pixel^2))*zu[i]
    } # i individuals
          # unmarked spatial count model likelihood
          for(j in 1:J){
            bigLambda[j] <- sum(lamu[1:m,j]) 
          for(k in 1:K){
            n[j,k] ~ dpois(bigLambda[j])
           } # k occasions
          } # j traps
          N <- sum(zu[1:m])
          D <- N/A
        })
  return(scrcode)
  } else   
     if(isFALSE(occ_specific)==FALSE){
        scrcode <- nimble::nimbleCode({
          lam0 ~ dunif(0,lam0_upper) # baseline encounter rate
          sigma ~ dunif(0, sigma_upper) # scaling parameter
          sigma.pixel <- sigma / pixelWidth # scaled for habitat mask
          psiu ~ dunif(0, 1) # inclusion prob
        for(i in 1:m){
          zu[i]~dbern(psiu)
          su[i,1] ~ dunif(x_lower, x_upper)
          su[i,2] ~ dunif(y_lower, y_upper)
          pOKu[i] <- hab_mask[(trunc(su[i,2])+1),(trunc(su[i,1])+1)] # habitat check
          OKu[i] ~ dbern(pOKu[i]) # OKu[i] <- 1, the ones trick
          distu[i,1:J] <- sqrt((su[i,1]-X[1:J,1])^2 + (su[i,2]-X[1:J,2])^2)
        for(k in 1:K){
          lamu[i,1:J,k] <- lam0*exp(-distu[i,1:J]^2/(2*sigma.pixel^2))*zu[i]
        } # k occasions
      } # i individuals
        # unmarked spatial count model likelihood
        for(j in 1:J){
          bigLambda[j] <- sum(lamu[1:m,j,1:K]) 
        for(k in 1:K){
          n[j,k] ~ dpois(bigLambda[j])
         } # k occasions
        } # j traps
        N <- sum(zu[1:m])
        D <- N/A
      })
   return(scrcode)
   }  # end occasion-specific
  } # end models with hab_mask
 }else # end trapsClustered == FALSE
 if(trapsClustered){
 if(isFALSE(hab_mask)){ # determine if hab_mask is included
 if(isFALSE(occ_specific)){
   scrcode <- nimble::nimbleCode({
  for(g in 1:nSites){  
    lam0[g] ~ dunif(0,lam0_upper) # baseline encounter rate
  }
    sigma ~ dunif(0, sigma_upper) # scaling parameter
    psiu ~ dunif(0, 1) # inclusion prob
  for(i in 1:m){
   zu[i]~dbern(psiu)
   su[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
   su[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
   distu[i,1:J] <- sqrt((su[i,1]-X[1:J,1,site[i]])^2 + (su[i,2]-X[1:J,2,site[i]])^2)
   lamu[i,1:J] <- lam0[site[i]]*exp(-distu[i,1:J]^2/(2*sigma^2))*zu[i]
  } # i individuals
  # unmarked spatial count model likelihood
  for(g in 1:nSites){
  for(j in 1:J){
    bigLambda[j,g] <- sum(lamu[site_indexL[g]:site_indexU[g],j]) 
  for(k in 1:K){
    n[j,k,g] ~ dpois(bigLambda[j,g])
   } # k occasions
  } # j traps
  } # g sites
  N <- sum(zu[1:m])
  D <- N/A
  })
  return(scrcode)
} else # end occasion-specific == FALSE
if(isFALSE(occ_specific)==FALSE){ # occasion-specific models
scrcode <- nimble::nimbleCode({
sigma ~ dunif(0, sigma_upper) # scaling parameter
psiu ~ dunif(0, 1) # inclusion prob
for(g in 1:nSites){
  lam0[g] ~ dunif(0,lam0_upper) # site-specific baseline encounter rate
} # g sites
for(i in 1:m){
  zu[i]~dbern(psiu)
  su[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
  su[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
  distu[i,1:J] <- sqrt((su[i,1]-X[1:J,1,site[i]])^2 + (su[i,2]-X[1:J,2,site[i]])^2)
 for(k in 1:K){
 lamu[i,1:J,k] <- lam0*exp(-distu[i,1:J]^2/(2*sigma^2))*zu[i]
  } # k occasions
} # i marked individuals
  # unmarked spatial count model likelihood
  for(g in 1:nSites){
  for(j in 1:J){
    bigLambda[j,g] <- sum(lamu[site_indexL[g]:site_indexU[g],j,1:K]) 
  for(k in 1:K){
    n[j,k,g] ~ dpois(bigLambda[j,g])
   } # k occasions
  } # j traps
  } # g sites
 N <- sum(zu[1:m])
 D <- N/A
})
return(scrcode)
} # end occasion-specific models
} else # end hab_mask == FALSE
if(hab_mask==TRUE){
if(isFALSE(occ_specific)){
scrcode <- nimble::nimbleCode({
for(g in 1:nSites){  
lam0[g] ~ dunif(0,lam0_upper) # baseline encounter rate
}
sigma ~ dunif(0, sigma_upper) # scaling parameter
sigma.pixel <- sigma / pixelWidth # scaled for habitat mask
psiu ~ dunif(0, 1) # inclusion prob
for(i in 1:m){
zu[i]~dbern(psiuu[i])
# adjust psiu for the proportion of available habitat at each site
psiuu[i] <- (1-(1-psiu)^prop.habitat[site[i]]) 
su[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
su[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
pOKu[i] <- hab_mask[(trunc(su[i,2])+1),(trunc(su[i,1])+1),site[i]] # habitat check
OKu[i] ~ dbern(pOKu[i]) # OKu[i] <- 1, the ones trick
distu[i,1:J] <- sqrt((su[i,1]-X[1:J,1,site[i]])^2 + (su[i,2]-X[1:J,2,site[i]])^2)
lamu[i,1:J] <- lam0[site[i]]*exp(-distu[i,1:J]^2/(2*sigma.pixel^2))*zu[i]
} # i individuals
# unmarked spatial count model likelihood
for(g in 1:nSites){
for(j in 1:J){
  bigLambda[j,g] <- sum(lamu[site_indexL[g]:site_indexU[g],j]) 
for(k in 1:K){
  n[j,k,g] ~ dpois(bigLambda[j,g])
 } # k occasions
} # j traps
} # g sites
N <- sum(zu[1:m])
D <- N/A
})
return(scrcode)
} else  # End isFALSE(occ_specific)
if(isFALSE(occ_specific)==FALSE){
scrcode <- nimble::nimbleCode({
  sigma ~ dunif(0, sigma_upper) # scaling parameter
  sigma.pixel <- sigma / pixelWidth # scaled for habitat mask
  psiu ~ dunif(0, 1) # inclusion prob
  for(g in 1:nSites){
    lam0[g] ~ dunif(0,lam0_upper) # site-specific baseline encounter rate
  } # g sites
  for(i in 1:m){
    zu[i]~dbern(psiuu[i])
   # adjust psiu for the proportion of available habitat at each site
  psiuu[i] <- (1-(1-psiu)^prop.habitat[site[i]]) 
  su[i,1] ~ dunif(x_lower[site[i]], x_upper[site[i]])
  su[i,2] ~ dunif(y_lower[site[i]], y_upper[site[i]])
  pOKu[i] <- hab_mask[(trunc(su[i,2])+1),(trunc(su[i,1])+1),site[i]] # habitat check
  OKu[i] ~ dbern(pOKu[i]) # OKu[i] <- 1, the ones trick
  distu[i,1:J] <- sqrt((su[i,1]-X[1:J,1,site[i]])^2 + (su[i,2]-X[1:J,2,site[i]])^2)
    for(k in 1:K){
      lamu[i,1:J,k] <- lam0[site[i]]*exp(-distu[i,1:J]^2/(2*sigma.pixel^2))*zu[i]
    } # k occasions
  } # i individuals
# unmarked spatial count model likelihood
for(g in 1:nSites){
for(j in 1:J){
  bigLambda[j,g] <- sum(lamu[site_indexL[g]:site_indexU[g],j,1:K]) 
for(k in 1:K){
  n[j,k,g] ~ dpois(bigLambda[j,g])
 } # k occasions
} # j traps
} # g sites
N <- sum(zu[1:m])
D <- N/A
})
return(scrcode)
} # end isFALSE(occ_specific)==FALSE
    } # end models with hab_mask
  } # end trapsClustered
} # End function 'get_unmarked"
