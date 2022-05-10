## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
# load 'localSCR' package
library(localSCR)

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # simulate a single trap array with random positional noise
#  x <- seq(-800, 800, length.out = 5)
#  y <- seq(-800, 800, length.out = 5)
#  traps <- as.matrix(expand.grid(x = x, y = y))
#  set.seed(200)
#  traps <- traps + runif(prod(dim(traps)),-20,20)
#  
#  mysigma = c(300) # simulate sex-specific scaling parameter
#  mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#  pixelWidth = 200 # store pixelWidth or grid resolution
#  
#  # create state-space
#  Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma,
#                      res = pixelWidth)
#  
#  # create polygon for mask
#  library(sf)
#  poly = st_sfc(st_polygon(x=list(matrix(c(-1765,-1765,1730,-1650,1600,
#                1650,0,1350,-800,1700,-1850,1000,-1765,-1765)
#              ,ncol=2, byrow=TRUE))), crs =  mycrs)
#  
#  # create habitat mask
#  hab_mask = mask_polygon(poly = poly, grid = Grid$grid, crs_ = mycrs,
#  prev_mask = NULL)
#  
#  # Simulated abundance
#  Nsim = 200
#  
#  # simulate data for uniform state-space and habitat mask
#  data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs,
#                       sigma_ = mysigma,  prop_sex = 1, N = Nsim, K = 4,
#                       base_encounter = 0.15,enc_dist = "binomial",
#                       hab_mask = hab_mask, setSeed = 200)
#  
#  # total augmented population size
#  M = 400
#  
#  # get initial activity center starting values
#  s.st = initialize_classic(y=data3d$y, M=M, X=traps, ext = Grid$ext,
#                              hab_mask = hab_mask)
#  
#  # convert traps and starting locations to discrete format
#  d_list = discretize_classic(X=traps, grid = Grid$grid, s.st = s.st,
#                              crs_ = mycrs, hab_mask = hab_mask)
#  
#  # inspect discrete data list
#  str(d_list)
#  #> List of 4
#  #>  $ grid: num [1:270, 1:2] -808 -608 1592 -1208 -1008 ...
#  #>   ..- attr(*, "dimnames")=List of 2
#  #>   .. ..$ : NULL
#  #>   .. ..$ : chr [1:2] "x" "y"
#  #>  $ nPix: int 270
#  #>  $ X   : num [1:25, 1:2] -807.71 -407.71 -7.71 392.29 792.29 ...
#  #>   ..- attr(*, "dimnames")=List of 2
#  #>   .. ..$ : NULL
#  #>   .. ..$ : chr [1:2] "x" "y"
#  #>  $ s.st: int [1:400] 174 194 54 178 130 91 156 163 190 193 ...
#  
#  # make ggplot of grid, and discretized trap locations
#  # and starting activity center locations
#  library(ggplot2)
#  ggplot() + geom_point(data=as.data.frame(d_list$grid),aes(x=x,y=y),
#                        color="grey60",size=2) +
#    geom_point(data=as.data.frame(d_list$X),aes(x=x,y=y),color="blue",size=3) +
#      geom_point(data=as.data.frame(d_list$grid[d_list$s.st,]),
#                 aes(x=x,y=y),color="orangered",size=0.75) +
#      theme_classic() + ylab("Northing") + xlab("Easting") +
#      theme(axis.text = element_text(size=12),
#            axis.title = element_text(size=16))

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # prepare data
#  data = list(y=data3d$y)
#  data$y = data$y[which(apply(data$y, 1, sum)!=0),,] # remove augmented records
#  # covert to 2d by summing over individuals and traps
#  data$y = apply(data$y, c(1,2), sum)
#  
#  # add discretized traps
#  data$X = d_list$X/1000 # convert to km units
#  
#  # add grid
#  data$grid = d_list$grid/1000 # convert to km units
#  
#  # prepare constants (note that density is now in activity centers/km2
#  # and each cell is now 0.01 km2 in area
#  constants = list(M = M,n0 = nrow(data$y),J=dim(data$y)[2],
#                   K=dim(data3d$y)[3],nPix=d_list$nPix,pixArea = (pixelWidth/1000)^2,
#                   sigma_upper = 1, A = sum(hab_mask)*((pixelWidth/1000)^2))
#  
#  # add z and zeros vector data for latent inclusion indicator
#  data$z = c(rep(1,constants$n0),rep(NA,constants$M - constants$n0))
#  data$zeros =  c(rep(NA,constants$n0),rep(0,constants$M - constants$n0))
#  
#  # define all initial values
#  inits = list(sigma = runif(1, 0.250, 0.350), s = d_list$s.st,
#              alpha0 = 2.8, p0 = runif(1, 0.05, 0.15),
#              z=c(rep(NA,constants$n0),rep(0,constants$M-constants$n0)))
#  
#  # parameters to monitor
#  params = c("sigma","psi","p0","N","D","EN","alpha0","s","z")
#  
#  # get model
#  discrete_model = get_discrete(type="marked",dim_y = 2,
#                enc_dist = "binomial",sex_sigma = FALSE,
#                trapsClustered=FALSE)
#  
#  # show model (not run)
#  # discrete_model
#  
#  # run model (note this was run on a Mac with 16 GB 2667 MHz DDR4
#  # and 2.3 GHz 8-Core Intel Core i9)
#  library(tictoc)
#  tic() # track time elapsed
#  out = run_discrete(model = discrete_model, data=data, constants=constants,
#          inits=inits, params = params,niter = 5000, nburnin=1000,
#          thin=1, nchains=2, parallel=TRUE,  RNGseed = 500)
#  toc()
#  #> 1025.869 sec elapsed
#  
#  # histogram of posterior samples for N (abundance)
#  samples = do.call(rbind, out)
#  par(mfrow=c(1,1))
#  hist(samples[,which(dimnames(out[[1]])[[2]]=="N")], xlab = "Abundance",
#       xlim = c(0,400), main="")
#  abline(v=Nsim, col="red") # add line for simulated abundance

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # summarize MCMC samples (exclude parameters and don't plot)
#  nimSummary(out, exclude = c("s","z"), trace=FALSE)
#  #>        post.mean post.sd    q2.5     q50   q97.5 f0   n.eff  Rhat
#  #> D         18.222   1.768  15.093  18.148  21.944  1 314.957 1.006
#  #> EN       196.003  21.492 157.391 194.404 241.130  1 307.866 1.006
#  #> N        196.798  19.091 163.000 196.000 237.000  1 314.957 1.006
#  #> alpha0     2.893   0.110   2.679   2.890   3.106  1 316.798 1.006
#  #> p0         0.167   0.023   0.126   0.167   0.214  1 316.707 1.016
#  #> psi        0.490   0.054   0.393   0.486   0.603  1 307.866 1.006
#  #> sigma      0.310   0.021   0.274   0.308   0.358  1 201.439 1.043
#  
#  # make realized density plot (we don't use the habitat mask here for
#  # discrete model)
#  r = realized_density(samples = out, grid = d_list$grid, crs_ = mycrs,
#                       site = NULL, hab_mask = FALSE,discrete=TRUE)
#  
#  # load virdiis color palette and raster libraries
#  library(viridis)
#  library(raster)
#  
#  # make simple raster plot
#  plot(r, col=viridis(100),
#       main=expression("Realized density (activity centers/0.2 km"^2*")"),
#       ylab="Northing",xlab="Easting")

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # create an array of traps, as an approach where individuals will only be
#  # detected at one of the trap arrays (e.g., Furnas et al. 2018)
#  Xarray = array(NA, dim=c(nrow(traps),2,2))
#  Xarray[,,1]=traps
#  Xarray[,,2]=traps+6000 # shift trapping grid to new locations
#  
#  # create grid and extent for 3D trap array
#  GridX = grid_classic(X = Xarray, crs_ = mycrs, buff = 3*max(mysigma),
#                       res = pixelWidth)
#  
#  # create polygon to use as a mask
#  library(sf)
#  poly = st_sfc(st_polygon(x=list(matrix(c(-1660,-1750,8730,-1050,7470,
#  7550,0,7950,-1800,8200,-1660,-1750),ncol=2, byrow=TRUE))), crs =  mycrs)
#  
#  # make ggplot
#  g1=ggplot() + geom_point(data=as.data.frame(GridX$grid[,,1]),
#                        aes(x=V1,y=V2),color="grey60",size=1.25) +
#      geom_point(data=as.data.frame(Xarray[,,1]),
#                 aes(x=V1,y=V2),color="blue",size=2) +
#      geom_point(data=as.data.frame(GridX$grid[,,2]),
#                 aes(x=V1,y=V2),color="grey60",size=1.25) +
#      geom_point(data=as.data.frame(Xarray[,,2]),
#                 aes(x=V1,y=V2),color="blue",size=2) +
#      geom_sf(data=poly, fill = NA) + coord_sf(datum=st_crs(mycrs)) +
#      theme_classic() + ylab("Northing") + xlab("Easting") +
#      scale_x_continuous(limits=c(-2000,9000)) +
#      scale_y_continuous(limits=c(-2000,9000)) +
#      theme(axis.text = element_text(size=12),
#                 axis.title = element_text(size=16))

## ---- fig.show='hide',fig.width = 8,fig.height=14,eval=FALSE------------------
#  # get 3D habitat mask array for 3D grid
#  hab_mask = mask_polygon(poly = poly, grid = GridX$grid, crs_ = mycrs,
#  prev_mask = NULL)
#  
#  # simuated population size
#  Nsim = 300
#  
#  # augmented population size
#  M=600
#  
#  # simulate data for uniform state-space and habitat mask
#  # (N is simulated abundance)
#  data4d = sim_classic(X = Xarray, ext = GridX$ext, crs_ = mycrs,
#                   sigma_ = mysigma, prop_sex = 1,N = Nsim,
#                   K = 4, base_encounter = 0.15,
#                   enc_dist = "binomial",hab_mask = hab_mask,
#                   setSeed = 300)
#  
#  # augment site identifier
#  site = c(data4d$site,c(rep(1,((M-length(data4d$site))/2)),
#                         rep(2,((M-length(data4d$site))/2))))
#  
#  # get initial activity center starting values
#  s.st = initialize_classic(y=data4d$y, M=M, X=Xarray, ext = GridX$ext,
#                              site = site, hab_mask = hab_mask)
#  
#  # convert traps and starting locations to discrete format
#  d_list = discretize_classic(X=Xarray, grid = GridX$grid, s.st = s.st,
#                              crs_ = mycrs,site=site, hab_mask = hab_mask)
#  
#  # inspect discrete data list
#  str(d_list)
#  #> List of 4
#  #>  $ grid: num [1:284, 1:2, 1:2] -1608 -1408 -1208 -1008 -808 ...
#  #>  $ nPix: int [1:2] 284 278
#  #>  $ X   : num [1:25, 1:2, 1:2] -807.71 -407.71 -7.71 392.29 792.29 ...
#  #>  $ s.st: num [1:600] 98 115 115 165 77 94 77 96 75 109 ...
#  
#  # prepare data
#  data = list(y=data4d$y)
#  data$y = data$y[which(apply(data$y, 1, sum)!=0),,] # remove augmented records
#  # covert to 2d by summing over individuals and traps
#  data$y = apply(data$y, c(1,2), sum)
#  
#  # add discretized traps
#  data$X = d_list$X #/1000 # convert to km units
#  
#  # add grid
#  data$grid = d_list$grid #/1000 # convert to km units
#  
#  # prepare constants (note that density is now in activity centers/km2
#  # and each cell is now 0.02 km2 in area
#  constants = list(M = M,n0 = nrow(data$y),J=dim(data$y)[2],site=site,
#                   K=dim(data4d$y)[3],nPix=sum(d_list$nPix),
#                   pixArea = (pixelWidth^2),
#                   sigma_upper = 1000,
#                   A = (sum(hab_mask)*((pixelWidth/1000)^2)),
#                   nSites = dim(d_list$X)[3])
#  
#  constants$npixSite = matrix(c(1,d_list$nPix[1],
#                                d_list$nPix[1]+1,
#                                sum(d_list$nPix)),
#                              ncol=2,byrow=TRUE)
#  
#  # add z and zeros vector data for latent inclusion indicator
#  data$z = c(rep(1,constants$n0),rep(NA,constants$M - constants$n0))
#  data$zeros =  c(rep(NA,constants$n0),rep(0,constants$M - constants$n0))
#  
#  # define all initial values
#  inits = list(sigma = runif(1, 250, 350), s = d_list$s.st,
#              alpha0 = -11.25, p0 = runif(dim(d_list$X)[3], 0.05, 0.15),
#              z=c(rep(NA,constants$n0),rep(0,constants$M-constants$n0)))
#  
#  # parameters to monitor
#  params = c("sigma","psi","p0","N","D","EN","alpha0","s","z")
#  
#  # get model
#  discrete_model = get_discrete(type="marked",dim_y = 2,
#                enc_dist = "binomial",sex_sigma = FALSE,
#                trapsClustered=TRUE)
#  
#  # show model (not run)
#  # discrete_model
#  
#  # run model (note this was run on a Mac with 16 GB 2667 MHz DDR4
#  # and 2.3 GHz 8-Core Intel Core i9)
#  library(tictoc)
#  tic() # track time elapsed
#  out = run_discrete(model = discrete_model, data=data, constants=constants,
#          inits=inits, params = params,niter = 5000, nburnin=1000,
#          thin=1, nchains=2, parallel=TRUE,  RNGseed = 500)
#  toc()
#  #> 1612.656 sec elapsed
#  
#  # histogram of posterior samples for N (abundance)
#  samples = do.call(rbind, out)
#  par(mfrow=c(1,1))
#  hist(samples[,which(dimnames(out[[1]])[[2]]=="EN")], xlab = "Expected abundance",
#       xlim = c(0,600), main="")
#  abline(v=Nsim, col="red") # add line for simulated abundance

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # summary table of MCMC output (exclude "s" and "z" parameters)
#  nimSummary(out, exclude = c("s","z"))
#  #>        post.mean post.sd    q2.5     q50   q97.5 f0   n.eff  Rhat
#  #> D         13.131   1.397  10.676  12.989  16.192  1 175.165 1.011
#  #> EN       294.564  33.831 232.157 291.304 367.164  1 178.972 1.011
#  #> N        295.180  31.395 240.000 292.000 364.000  1 175.165 1.011
#  #> alpha0   -11.249   0.114 -11.481 -11.254 -11.022  0 188.347 1.013
#  #> p0[1]      0.148   0.025   0.104   0.147   0.199  1 294.255 1.033
#  #> p0[2]      0.128   0.023   0.088   0.127   0.179  1 324.642 1.032
#  #> psi        0.491   0.056   0.387   0.486   0.612  1 178.972 1.011
#  #> sigma    287.695  23.052 250.820 284.677 339.804  1 139.887 1.103
#  
#  # generate realized density surface
#  r = realized_density(samples=out, grid=d_list$grid, ext = GridX$ext,
#              crs_=mycrs,site=constants$site, discrete=TRUE)
#  
#  # load needed packages for multiplot
#  library(viridis)
#  library(grid)
#  library(cowplot)
#  library(ggpubr)
#  library(rasterVis)
#  
#  # plot raster from site 1
#  p1<-gplot(r[[1]]) + geom_raster(aes(fill = value)) +
#            scale_fill_viridis(na.value = NA, name="Density",
#            limits=c(0,1.5),breaks=seq(0,1.5,by=0.5)) +
#            xlab("") + ylab("") + theme_classic() +
#            scale_x_continuous(expand=c(0, 0)) +
#            scale_y_continuous(expand=c(0, 0)) +
#             theme(axis.text = element_text(size=18))
#  
#  # plot raster from site 2
#  p2<-gplot(r[[2]]) + geom_raster(aes(fill = value)) +
#            scale_fill_viridis(na.value = NA, name="Density",
#            limits=c(0,1.5),breaks=seq(0,1.5,by=0.5)) +
#            xlab("") + ylab("") + theme_classic() +
#            scale_x_continuous(expand=c(0, 0)) +
#            scale_y_continuous(expand=c(0, 0)) +
#            theme(axis.text = element_text(size=18))
#  
#  # arrange the two plots in a single row
#  prow <- plot_grid(p1 + theme(legend.position="none"),
#             p2 + theme(legend.position="none"),
#             align = 'vh',
#             labels = NULL,
#             hjust = -1,
#             nrow = 1
#             )
#  
#  # extract the legend from one of the plots
#  legend_t <- get_legend(p1 + theme(legend.position = "top",
#                          legend.direction = "horizontal",
#                          legend.text = element_text(size=14),
#                          legend.title = element_text(size=16)))
#  
#  # add the legend above the row we made earlier. Give it 20% of the height
#  # of one plot (via rel_heights).
#  pcomb <- plot_grid(legend_t, prow, ncol = 1, rel_heights = c(.2, 1))
#  
#  # add x and y axis labels
#  pcomb <-annotate_figure(pcomb, bottom = textGrob("Easting",
#                gp=gpar(fontsize=18), vjust = -1, hjust = 0),
#                left = textGrob("Northing", rot=90, gp=gpar(fontsize=18),
#                vjust = 1, hjust = 0.5))
#  pcomb

