% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/localSCR.R
\name{localize_classic}
\alias{localize_classic}
\title{Prepare classic SCR data components for local modeling approach}
\usage{
localize_classic(
  y,
  grid_ind,
  X,
  crs_,
  sigma_,
  s.st,
  site = NULL,
  hab_mask = FALSE
)
}
\arguments{
\item{y}{Either a matrix or array of encounter history data, possibly 
from \code{sim_classic()}.}

\item{grid_ind}{A matrix representing an individual state-space grid. This 
is returned from \code{\link{grid_classic}}.}

\item{X}{Either a matrix or array representing the coordinates of traps in
UTMs. An array is used when traps are clustered over a survey area.}

\item{crs_}{The UTM coordinate reference system (EPSG code) used for your
location provided as an integer (e.g., 32608 for WGS 84/UTM Zone 8N).}

\item{sigma_}{The scaling parameter of the bivariate normal kernel either
in meters or kilometers as an integer. Note that if \code{sigma_} is
sex-specific, the maximum value will be used.}

\item{s.st}{A matrix of starting activity center coordinates. This is 
returned from \code{\link{initialize_classic}}.}

\item{site}{Either \code{NULL} (if a 2D trap array is used) or a vector of 
integers denoting which trap array an individual (either detected or 
augmented) belongs to. Note that \code{site} is provided from 
\code{\link{sim_classic}} when a 3D trap array is used. However, this
 \code{site} variable must be correctly augmented based on the total 
 augmented population size (i.e., \code{M}).}

\item{hab_mask}{Either \code{FALSE} (the default) or a matrix or array output
from \code{\link{mask_polygon}} or \code{\link{mask_raster}} functions.}
}
\value{
A list of data components needed to for classic SCR models in a local
approach. Specifically, the function returns:
\itemize{
    \item{\code{y}}{ Individual-specific encounter history that only considers traps 
    within a distance of 9 times \code{sigma_}.}
    \item{\code{X}}{ Individual-specific trap array that only provides coordinates 
    for traps within a distance of 9 times \code{sigma_}.}
    \item{\code{grid}}{ A matrix or array of state-space grid coordinates that 
    encompasses all individual state-space grids.}
    \item{\code{prop_habitat}}{ The proportion of suitable habitat for each individual. 
    Note that this is only returned when a habitat mask is used.}
    \item{\code{ext_mat}}{ A matrix of individual state-space grid extents.}
    \item{\code{ext}}{ An \code{Extent} object from the \code{raster} package for the 
    scaled-up state-space grid.}
    \item{\code{Jind}}{ A vector with the number of traps each individual is exposed to.}
    \item{\code{s.st}}{ A matrix of starting activity center coordinates.}
}
}
\description{
A function to scale up from individual-level state-space to 
study-area level state-space using classic SCR data components.
}
\details{
This function converts classic SCR data to a format that is used in 
local evaluations of the individual state-space.
}
\examples{
\dontrun{
# simulate a single trap array with random positional noise
x <- seq(-1600, 1600, length.out = 6)
y <- seq(-1600, 1600, length.out = 6)
traps <- as.matrix(expand.grid(x = x, y = y))
# add some random noise to locations
traps <- traps + runif(prod(dim(traps)),-20,20) 
mysigma = 300 # simulate sigma of 300 m
mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
pixelWidth = 100 # grid resolution

# Create initial grid and extent (use a slightly bigger buffer to match
# scaled-up state-space below)
Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3.7*mysigma, res = pixelWidth)

# Simulated abundance
Nsim = 250

# simulate SCR data
data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, 
sigma_ = mysigma, prop_sex = 1, N = Nsim, K = 4, base_encounter = 0.15, 
enc_dist = "binomial", hab_mask = FALSE, setSeed = 50)

# generate initial activity center coordinates for 2D trap array without 
# habitat mask
s.st = initialize_classic(y=data3d$y, M=500, X=traps, ext = Grid$ext, 
hab_mask = FALSE, all_random = FALSE)

# now use grid_classic to create an individual-level state-space (with origin 0, 0)
Grid_ind = grid_classic(X = matrix(c(0,0),nrow=1), crs_ = mycrs, 
          buff = 3*mysigma, res = 100)

# now localize the data components created above
local_list = localize_classic(y = data3d$y, grid_ind = Grid_ind$grid, X=traps, 
                           crs_ = mycrs, sigma_ = mysigma, s.st = s.st,
                           site = NULL, hab_mask = FALSE)
# inspect local list
str(local_list)
}
}
\author{
Daniel Eacker
}
