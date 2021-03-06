% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/localSCR.R
\name{sim_classic}
\alias{sim_classic}
\title{Simulate basic spatial capture-recapture data}
\usage{
sim_classic(
  X,
  ext,
  crs_,
  N,
  sigma_,
  prop_sex,
  K,
  base_encounter,
  enc_dist = "binomial",
  hab_mask = FALSE,
  setSeed = 500
)
}
\arguments{
\item{X}{Either a matrix or array object representing the coordinates of traps in
UTMs. An array is used when traps are clustered over a survey area.}

\item{ext}{An \code{Extent} object from the \code{raster} package. This is 
returned from \code{\link{grid_classic}}.}

\item{crs_}{The UTM coordinate reference system (EPSG code) used for your
location provided as an integer (e.g., 32608 for WGS 84/UTM Zone 8N).}

\item{N}{Simulated total abundance as an integer.}

\item{sigma_}{The scaling parameter of the bivariate normal kernel either
in meters or kilometers as an integer.}

\item{prop_sex}{The portion of females or males as a numeric value. This
will depend upon the indicator coding scheme used (e.g., females = 1 and
males = 0; then proportion of females in the simulation). Must be a numeric
value between 0 and 1. Note that 0 or 1 can be used if a non-sex-specific 
sigma is desired.}

\item{K}{The number of sampling occasions desired as an integer.}

\item{base_encounter}{The baseline encounter probability or rate as a numeric
value. Note that a probabilty is given for a \code{"binomial"} observation
distribution while a rate is given for a \code{"poisson"} distribution.}

\item{enc_dist}{Either \code{"binomial"} or \code{"poisson"}. Default is
\code{"binomial"}.}

\item{hab_mask}{Either \code{FALSE} (the default) or a matrix or arrary output from \code{\link{mask_polygon}}
or \code{\link{mask_raster}} functions.}

\item{setSeed}{The random number generater seed as an integer used in
simulations to obtain repeatable data simulations. Default is 500.}
}
\value{
\itemize{
 \item{\code{y}}{ A list of a matrix or array of encounter histories.}
 \item{\code{sex}}{ A vector or matrix of 0's and 1's for sex identification.}
 \item{\code{s}}{ A matrix of simulated activity centers.}
 \item{\code{site}}{ A vector for the site identifier.}
}
}
\description{
Returns a list of simulated data including the encounter
history, binary sex indicator, activity centers, and site identifier.
}
\details{
This function supports spatial capture-recapture (SCR) analysis by 
allowing for easy simulation of data components used by nimble in Baysian 
SCR models. Note that the output for the encounter histories \code{y} will be
sorted by detected and not detected individuals.
}
\note{
The \code{site} identfier is only returned when a 3-dimensional 
trap array is provided.
}
\examples{
# simulate a single trap array with random positional noise
x <- seq(-800, 800, length.out = 5)
y <- seq(-800, 800, length.out = 5)
traps <- as.matrix(expand.grid(x = x, y = y))

# add some random noise to locations
traps <- traps + runif(prod(dim(traps)),-20,20) 

mysigma = 300 # simulate sigma of 300 m
mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N

# Create grid and extent
Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = 100)

# simulate SCR data
data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, sigma_ = c(300), 
prop_sex = 1, N = 200, K = 4, base_encounter = 0.25, enc_dist = "binomial", 
hab_mask = FALSE, setSeed = 50)

# make simple plot
par(mfrow=c(1,1))
plot(Grid$grid, pch=20,ylab="Northing",xlab="Easting")
points(traps, col="blue",pch=20)
points(data3d$s,col="red",pch = 20)
points(data3d$s[which(apply(data3d$y,1,sum)!=0),],col="green",pch = 20)
}
\seealso{
\code{\link{grid_classic}}
}
\author{
Daniel Eacker
}
