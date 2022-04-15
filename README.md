
<!-- README.md is generated from README.Rmd. Please edit that file -->

# localSCR

<!-- badges: start -->

<!-- badges: end -->

The goal of ‘localSCR’ is to provide user-friendly functions to
implement spatially-explicit mark recapture (SECR) models using the
‘nimble’ package (de Valpine et a. 2022) in R. The package currently
has functions to 1) assist with defining the state-space grid and extent
for a given 2-dimensional or 3-dimensional trap array (i.e., when traps
are clustered in space), 2) simulate data under different encounter
distributions and other parameters, 3) create habitat masks from either
raster data or spatial polygons, 4) provide template SECR models that
are easily customizable, 5) fit and summarize SECR models using ‘nimble’
(de Valpine et a. 2022) with options for parallel processing, and 6)
create realized density surfaces from MCMC output. Future functionality
will include discrete state-space models and implementing localized
approaches as in Milleret et al. (2019) and Woodruff et al. (2020). The
package also uses block updating of x and y activity center coordinates,
vectorization of traps in derivations, and separating the data
augmentation process (see Chandler 2018) to decrease run time.

Another useful package is ‘nimbleSCR’ (Bischof et al. 2021) that
implements custom sampling distributions to increase sampling speed and
efficiency and recent methods of localized activity area approaches (see
Milleret et al. 2019) along with other custom functionality. These
custom distributions from ‘nimbleSCR’ can be used within our template
models, but the idea with the ‘localSCR’ package is to provide a more
user-friendly implemention of SECR models in ‘nimble’ similar to the
functionality provided by the ‘jagsUI’ package (Keller 2019) for JAGS
(Plummer 2017). Our intention is to ease the transition for users that
were previously using JAGS or those that have yet to switch over to
Bayesian methods from frequentist because of the complexity of
implementing these approaches.

## Installation

You can install the development version of ‘localSCR’ like so:

``` r
library(remotes)
install_github("sitkensis22/localSCR")
#> Skipping install of 'localSCR' from a github remote, the SHA1 (3499af64) has not changed since last install.
#>   Use `force = TRUE` to force installation
```

Be sure to see important information about using ‘nimble’ on your
computer (including installing rtools): <https://r-nimble.org/download>

## Example

This example includes constructing the state-space, simulating SECR (or
SCR) data, and the entire workflows under a uniform state-space
assumption with a habitat mask for 2- and 3-dimensional trap arrays:

``` r

# load package
library(localSCR)
```

### (1) Simulate a single trap array with random positional noise and create state-space

``` r

# simulate a single trap array with random positional noise
x <- seq(-800, 800, length.out = 5)
y <- seq(-800, 800, length.out = 5)
traps <- as.matrix(expand.grid(x = x, y = y))
traps <- traps + runif(prod(dim(traps)),-20,20) # add some random noise to locations

mysigma = 300 # simulate sigma of 300 m
mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N

# create state-space
Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = 100)

# make plot of grid and trap locations
par(mfrow=c(1,1))
plot(Grid$grid, pch=20)
points(traps, col="blue",pch=19)
```

<p align="center">

<img src="man/figures/Fig1.png" />

</p>

### (2) Simulate SCR data and make a plot of it.

``` r

# simulate SCR data
data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, sigma_ = mysigma, prop_sex = 1,
N = 200, K = 4, base_encounter = 0.25, enc_dist = "binomial", hab_mask = FALSE, setSeed = 50)

# inspect simulated data
str(data3d)
#> List of 3
#>  $ y  : int [1:200, 1:25, 1:4] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ sex: int [1:200] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ s  : num [1:200, 1:2] 718.2 -210.1 -1023.9 917.9 48.5 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "sx" "sy"

# make simple plot
par(mfrow=c(1,1))
plot(Grid$grid, pch=20,ylab="Northing",xlab="Easting")
points(traps, col="blue",pch=20)
points(data3d$s,col="red",pch = 20) # all simulated activity centers
points(data3d$s[which(apply(data3d$y,1,sum)!=0),],col="green",pch = 20) # detected individuals
```

<p align="center">

<img src="man/figures/Fig2.png" />

</p>

### (3) Workflow for simple SCR model with sex-specific sigma, binomial encounter distribution, and habitat mask

``` r

# simulate a single trap array with random positional noise
x <- seq(-800, 800, length.out = 5)
y <- seq(-800, 800, length.out = 5)
traps <- as.matrix(expand.grid(x = x, y = y))
traps <- traps + runif(prod(dim(traps)),-20,20) # add some random noise to locations

mysigma = c(220, 300) # simulate sex-specific
mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
pixelWidth = 100 # store pixelWidth

# create state-space grid and extent
Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*max(mysigma), res = pixelWidth) 

# create polygon for mask
library(sf)
#> Linking to GEOS 3.9.1, GDAL 3.2.1, PROJ 7.2.1; sf_use_s2() is TRUE
poly = st_sfc(st_polygon(x=list(matrix(c(-1765,-1765,1730,-1650,1600,1650,0,1350,-800,1700,
-1850,1000,-1765,-1765),ncol=2, byrow=TRUE))), crs =  mycrs)

# create habitat mask
hab_mask = mask_polygon(poly = poly, grid = Grid$grid, crs_ = mycrs, prev_mask = NULL)

# simulate data for uniform state-space and habitat mask
data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, sigma_ = mysigma, prop_sex = 0.7,
N = 200, K = 4, base_encounter = 0.15, enc_dist = "binomial",hab_mask = hab_mask, setSeed = 50)

# total augmented population size 
M = 400

# get initial activity center starting values
s.st3d = initialize_classic(y=data3d$y, M=M, X=traps, buff = 3*max(mysigma), hab_mask = hab_mask)

# make simple plot
par(mfrow=c(1,1))
plot(Grid$grid, pch=20, xlim=c(-2000,2000),ylim=c(-2000,2000))
points(traps, col="blue",pch=19)
plot(poly, add=TRUE)
points(s.st3d,col="red",pch=20) # all initalized activity centers
```

<p align="center">

<img src="man/figures/Fig3.png" />

</p>

``` r
# rescale inputs
rescale_list = rescale_classic(X = traps, ext = Grid$ext, s.st = s.st3d, hab_mask = hab_mask)

# store rescaled extent
ext = rescale_list$ext

# prepare data
data = list(y=data3d$y)
data$y = data$y[which(apply(data$y, 1, sum)!=0),,] # remove augmented records 
data$y = apply(data$y, c(1,2), sum) # covert to 2d by summing over occasions

# add rescaled traps
data$X = rescale_list$X

# prepare constants (note get density in activity center/100 m2 rather than activity centers/m2)
constants = list(M = M,n0 = nrow(data$y),J=dim(data$y)[2], K=dim(data3d$y)[3],
x_lower = ext[1], x_upper = ext[2], y_lower = ext[3], y_upper = ext[4],
sigma_upper = 1000, A = (sum(hab_mask)*(pixelWidth/100)^2),pixelWidth=pixelWidth)

# augment sex
data$sex = c(data3d$sex,rep(NA,constants$M-length(data3d$sex)))

# add z and zeros vector data for latent inclusion indicator
data$z = c(rep(1,constants$n0),rep(NA,constants$M - constants$n0))
data$zeros =  c(rep(NA,constants$n0),rep(0,constants$M - constants$n0))

# add hab_mask and OK for habitat check
data$hab_mask = hab_mask
data$OK = rep(1,constants$M)

# get initial activity center starting values
s.st3d = rescale_list$s.st

# define all initial values
inits = list(sigma = runif(2, 250, 350), s = s.st3d,psi=runif(1,0.2,0.3),
p0 = runif(1, 0.05, 0.15),pOK=data$OK,z=c(rep(NA,constants$n0),rep(0,constants$M-constants$n0)))

# parameters to monitor
params = c("sigma","psi","p0","N","D","psi_sex","s","z")

# get model
scr_model = get_classic(dim_y = 2, enc_dist = "binomial",sex_sigma = TRUE,hab_mask=TRUE,trapsClustered=FALSE)

# run model
library(tictoc)
tic() # track time elapsed
out = run_classic(model = scr_model, data=data, constants=constants,
inits=inits, params = params,niter = 10000, nburnin=1000, thin=1, nchains=2, parallel=TRUE, RNGseed = 500)
toc()
#> 159.38 sec elapsed

# summarize output
samples = do.call(rbind, out)
par(mfrow=c(1,1))
hist(samples[,which(dimnames(out[[1]])[[2]]=="N")], xlab = "Abundance", xlim = c(0,500), main="")
abline(v=200, col="red") # add line for simulated abundance
```

<p align="center">

<img src="man/figures/Fig4.png" />

</p>

``` r

# summarize MCMC samples (exclude parameters and don't plot)
nimSummary(out, exclude_params = c("s","z"), trace=FALSE)
#>          post.mean post.sd    q2.5     q50   q97.5 f0   n.eff  Rhat
#> D            0.176   0.021   0.140   0.174   0.222  1 528.137 1.054
#> N          187.648  22.117 150.000 186.000 237.000  1 528.137 1.054
#> p0           0.163   0.029   0.113   0.161   0.226  1 473.142 1.003
#> psi          0.469   0.060   0.367   0.464   0.601  1 433.830 1.057
#> psi_sex      0.660   0.068   0.513   0.665   0.782  1 371.001 1.029
#> sigma[1]   242.149  30.309 191.889 239.089 311.141  1 410.051 1.036
#> sigma[2]   289.258  26.584 246.359 286.377 349.280  1 307.972 1.031

# make realized density plot 
r = realized_density(samples = out, grid = Grid$grid, crs_ = mycrs, site = NULL, hab_mask = hab_mask)       
      
# load virdiis color palette and raster libraries      
library(viridis)
#> Loading required package: viridisLite
library(raster)
#> Loading required package: sp

# make simple raster plot
plot(r, col=viridis(100),main=expression("Realized density (activity centers/100 m"^2*")"),
     ylab="Northing",xlab="Easting")
```

<p align="center">

<img src="man/figures/Fig5.png" />

</p>

### (4) Workflow for simple SCR model with sex-specific sigma, binomial encounter distribution, and habitat mask using a 3D trap array or clustered traps.

``` r

# simulate a single trap array with random positional noise
x <- seq(-800, 800, length.out = 5)
y <- seq(-800, 800, length.out = 5)
traps <- as.matrix(expand.grid(x = x, y = y))
traps <- traps + runif(prod(dim(traps)),-20,20) # add some random noise to locations

mysigma = c(220, 300) # simulate sex-specific
mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
pixelWidth = 100 # store pixelWidth

# create an array of traps, as an approach where individuals will only be detected 
# at one of the trap arrays (e.g., Furnas et al. 2018)
Xarray = array(NA, dim=c(nrow(traps),2,2))
Xarray[,,1]=traps
Xarray[,,2]=traps+4000 # shift trapping grid to new locations

# create grid and extent for 3D trap array
GridX = grid_classic(X = Xarray, crs_ = mycrs, buff = 3*max(mysigma), res = 100)

# make simple plot
par(mfrow=c(1,1))
plot(GridX$grid[,,1],xlim=c(-1600,6000),ylim=c(-1600,6000),col="darkgrey",
pch=20,ylab="Northing",xlab="Easting")
points(Xarray[,,1],col="blue",pch=20)
points(GridX$grid[,,2],pch=20,col="darkgrey")
points(Xarray[,,2],col="blue",pch=20)

# create polygon to use as a mask
library(sf)
poly = st_sfc(st_polygon(x=list(matrix(c(-1660,-1900,5730,-1050,5470,
5650,0,6050,-1800,5700,-1660,-1900),ncol=2, byrow=TRUE))), crs =  mycrs)

# add polygon to plot
plot(poly, add=TRUE)
```

<p align="center">

<img src="man/figures/Fig6.png" />

</p>

``` r

# get 3D habitat mask array for 3D grid
hab_mask = mask_polygon(poly = poly, grid = GridX$grid, crs_ = mycrs, prev_mask = NULL)

# simulate data for uniform state-space and habitat mask (N is simulated abundance per site)
data4d = sim_classic(X = Xarray, ext = GridX$ext, crs_ = mycrs, sigma_ = mysigma, prop_sex = 0.7,
N = 200, K = 4, base_encounter = 0.15, enc_dist = "binomial",hab_mask = hab_mask, setSeed = 500)

# total augmented population size 
M = 400

# augment site identifier
site = c(data4d$site,c(rep(1,((M-length(data4d$site))/2)),rep(2,((M-length(data4d$site))/2))))

# get initial activity center starting values 
s.st4d = initialize_classic(y=data4d$y, M=M, X=Xarray, buff = 3*max(mysigma), site = site, hab_mask = hab_mask)

# rescale inputs
rescale_list = rescale_classic(X = Xarray, ext = GridX$ext, s.st = s.st4d, site = site, hab_mask = hab_mask)

# store rescaled extent and convert to matrix
ext = do.call(rbind, lapply(rescale_list$ext, as.vector))

# prepare constants (note get density in activity center/100 m2 rather than activity centers/m2)
constants = list(M = M,n0 =  length(which(apply(data4d$y,1,sum)!=0)),
J=dim(data4d$y)[2], K=dim(data4d$y)[3], sigma_upper = 1000, A = (sum(hab_mask)*(pixelWidth/100)^2),
pixelWidth=pixelWidth,nSites=dim(Xarray)[3],site = site)

# prepare data
data = list(X = rescale_list$X,sex = c(data4d$sex,rep(NA,M-length(data4d$sex))),
x_lower = ext[,1],x_upper = ext[,2],y_lower = ext[,3],y_upper = ext[,4])

# store and format encounter history data
data$y = data4d$y[which(apply(data4d$y, 1, sum)!=0),,] # remove augmented records 
data$y = apply(data$y, c(1,2), sum) # covert to 2d by summing over occasions

# add z and zeros vector data for latent inclusion indicator
data$z = c(rep(1,constants$n0),rep(NA,constants$M - constants$n0))
data$zeros =  c(rep(NA,constants$n0),rep(0,constants$M - constants$n0))

# add hab_mask, proportion of available habitat, and OK for habitat check
data$hab_mask = hab_mask
data$prop.habitat=apply(hab_mask,3,mean) # need to adjust proportion of habitat available
data$OK = rep(1,constants$M)

# get initial activity center starting values
s.st = rescale_list$s.st

# define all initial values
inits = list(sigma = runif(2, 250, 350), s = s.st,psi=runif(1,0.2,0.3),
p0 = runif(dim(data$X)[3], 0.1, 0.2),sex=ifelse(is.na(data$sex),rbinom(constants$M-constants$n0,1,0.5),NA),
pOK=data$OK,z=c(rep(NA,constants$n0),rep(0,constants$M-constants$n0)),
psi_sex=runif(1,0.4,0.6))

# parameters to monitor
params = c("sigma","psi","p0","N","D","psi_sex","s","z")

# get model
scr_model = get_classic(dim_y = 2, enc_dist = "binomial",sex_sigma = TRUE,hab_mask=TRUE,trapsClustered = TRUE)

# run model
library(tictoc)
tic() # track time elapsed
out = run_classic(model = scr_model, data=data, constants=constants,
inits=inits, params = params,niter = 10000, nburnin=1000, thin=1, nchains=2, 
parallel=TRUE, RNGseed = 500)
toc()
#> 171.27 sec elapsed

# summary table of MCMC output (exclude "s" and "z" parameters)
nimSummary(out, exclude_params = c("s","z"), trace = TRUE, plot_all = FALSE)
#>          post.mean post.sd    q2.5     q50   q97.5 f0   n.eff  Rhat
#> D            0.104   0.017   0.076   0.102   0.143  1 173.373 1.012
#> N          232.925  38.702 170.000 228.000 321.000  1 173.373 1.012
#> p0[1]        0.112   0.027   0.067   0.110   0.175  1 434.862 1.004
#> p0[2]        0.106   0.023   0.068   0.104   0.157  1 545.243 1.023
#> psi          0.592   0.100   0.425   0.582   0.815  1 179.821 1.010
#> psi_sex      0.509   0.080   0.354   0.509   0.664  1 266.020 1.001
#> sigma[1]   234.625  34.141 177.724 231.260 311.807  1 274.116 1.004
#> sigma[2]   344.204  38.158 282.668 339.885 430.743  1 263.555 1.027
```

<p align="center">

<img src="man/figures/Fig7.png" />

</p>

``` r

# generate realized density surface
r = realized_density(samples=out, grid=GridX$grid, crs_=mycrs,
                      site=constants$site, hab_mask=hab_mask)

# load needed packages for multiplot
library(viridis) 
library(grid)
library(cowplot)
library(ggpubr) 
#> Loading required package: ggplot2
#> 
#> Attaching package: 'ggpubr'
#> The following object is masked from 'package:cowplot':
#> 
#>     get_legend
#> The following object is masked from 'package:raster':
#> 
#>     rotate
library(rasterVis)
#> Loading required package: lattice

# plot raster from site 1
p1<-gplot(r[[1]]) + geom_raster(aes(fill = value)) +
          scale_fill_viridis(na.value = NA, name="Density",
          limits=c(0,0.3),breaks=seq(0,0.3,by=0.1)) +
          xlab("") + ylab("") + theme_classic() +
          scale_x_continuous(expand=c(0, 0)) + 
          scale_y_continuous(expand=c(0, 0)) + 
           theme(axis.text = element_text(size=18))

# plot raster from site 2
p2<-gplot(r[[2]]) + geom_raster(aes(fill = value)) +
          scale_fill_viridis(na.value = NA, name="Density",
          limits=c(0,0.3),breaks=seq(0,0.3,by=0.1)) +
          xlab("") + ylab("") + theme_classic() + 
          scale_x_continuous(expand=c(0, 0)) + 
          scale_y_continuous(expand=c(0, 0)) + 
          theme(axis.text = element_text(size=18))

# arrange the two plots in a single row
prow <- plot_grid(p1 + theme(legend.position="none"),
           p2 + theme(legend.position="none"),
           align = 'vh',
           labels = NULL,
           hjust = -1,
           nrow = 1
           )
#> Warning: Removed 13 rows containing missing values (geom_raster).
#> Warning: Removed 54 rows containing missing values (geom_raster).

# extract the legend from one of the plots
legend_t <- get_legend(p1 + theme(legend.position = "top",legend.direction = "horizontal",legend.text = element_text(size=14),legend.title = element_text(size=16)))
#> Warning: Removed 13 rows containing missing values (geom_raster).

# add the legend above the row we made earlier. Give it 20% of the height
# of one plot (via rel_heights).
pcomb <- plot_grid(legend_t, prow, ncol = 1, rel_heights = c(.2, 1))

# add x and y axis labels
pcomb <-annotate_figure(pcomb, bottom = textGrob("Easting", gp=gpar(fontsize=18), vjust = -1, hjust = 0),
    left = textGrob("Northing", rot=90, gp=gpar(fontsize=18),
         vjust = 1, hjust = 0.5))
pcomb
```

<p align="center">

<img src="man/figures/Fig8.png" />

</p>

## Literature Cited

Bischof R., D. Turek, C. Milleret, T. Ergon, P. Dupont, S. Dey, W. Zhang
and P. de Valpine. 2021. nimbleSCR: Spatial Capture-Recapture (SCR)
Methods Using ‘nimble’. R package version 0.1.3.
<https://CRAN.R-project.org/package=nimbleSCR>.

Chandler, R. B. 2018. Speeding up data augmentation in BUGS.
<https://groups.google.com/forum/#!topic/hmecology/o6cWDqHHgOE>.

de Valpine P, C. Paciorek, D. Turek, N. Michaud, C. Anderson-Bergman, F.
Obermeyer, C. C. Wehrhahn, A. Rodrìguez, L. D. Temple, and S. Paganin.
2022. *NIMBLE: MCMC, Particle Filtering, and Programmable Hierarchical
Modeling*. doi: 10.5281/zenodo.1211190 (URL:
<https://doi.org/10.5281/zenodo.1211190>), R package version 0.12.2,
URL:<https://cran.r-project.org/package=nimble>.

Furnas, B. J., R. H. Landers, S. Hill, S. S. Itoga, and B. N. Sacks.
2018. Integrated modeling to estimate population size and composition of
mule deer. Journal of Wildlife Management 82:1429–1441.

Kellner, K. 2018. jagsUI: a wrapper around ‘rjags’ to streamline ‘JAGS’
analyses. R package version 1.5.0.
<https://CRAN.R-project.org/package=jagsUI>.

Milleret, C., P. Dupont, C. Bonenfant, H. Henrik Brøseth, Ø. Flagstad,
C. Sutherland, and R. Bischof. 2019. A local evaluation of the
individual state‐space to scale up Bayesian spatial capture‐recapture.
Ecology and Evolution 9:352–363.

Plummer, M. 2017. JAGS version 4.3.0. user manual.
<https://people.stat.sc.edu/hansont/stat740/jags_user_manual.pdf>.

Woodruff, S., D. R. Eacker, and L. Waits. 2020. Estimating coyote
density in local, discrete Bayesian capture-recapture models. Journal of
Wildlife Management 10.1002/jwmg.21967
