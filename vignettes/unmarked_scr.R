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
#  traps <- traps + runif(prod(dim(traps)),-20,20)
#  
#  mysigma = 300 # simulate sigma of 300 m
#  mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#  
#  # create state-space
#  Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = 100)
#  
#  # make ggplot of grid and trap locations
#  library(ggplot2)
#  ggplot() + geom_point(data=as.data.frame(Grid$grid),aes(x=x,y=y),
#                        color="grey60", size=1.25) +
#      geom_point(data=as.data.frame(traps),aes(x=x,y=y),color="blue",size=2) +
#      theme_classic() + ylab("Northing") + xlab("Easting") +
#      scale_x_continuous(expand=c(-0.1, 0.1)) +
#      scale_y_continuous(expand=c(-0.1, 0.1)) +
#      theme(axis.text = element_text(size=12),axis.title = element_text(size=16))

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # simulate SCR data
#  data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs,
#                       sigma_ = mysigma, prop_sex = 1,N = 200, K = 4,
#                       base_encounter = 0.10, enc_dist = "poisson",
#                       hab_mask = FALSE, setSeed = 100)
#  
#  # inspect simulated data
#  str(data3d)
#  #> List of 3
#  #>  $ y  : int [1:200, 1:25, 1:4] 0 0 0 0 0 0 0 0 0 0 ...
#  #>  $ sex: int [1:200] 1 1 1 1 1 1 1 1 1 1 ...
#  #>  $ s  : num [1:200, 1:2] -662 -834 177 -1526 -110 ...
#  #>   ..- attr(*, "dimnames")=List of 2
#  #>   .. ..$ : NULL
#  #>   .. ..$ : chr [1:2] "sx" "sy"
#  
#  # We sum over traps and occasions to produce a 2-dimensional spatial count data set
#  n = apply(data3d$y, c(2,3), sum)
#  
#  # inspect n[j,k]
#  str(n)
#  #> int [1:25, 1:4] 0 3 2 0 1 1 1 1 1 1 ...

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # simulate a dense, spatially autocorrelated trap array with random positional noise
#  x <- seq(-800, 800, length.out = 5)
#  y <- seq(-800, 800, length.out = 5)
#  traps <- as.matrix(expand.grid(x = x, y = y))
#  traps <- traps + runif(prod(dim(traps)),-20,20) # add some random noise to locations
#  
#  mysigma = 300 # simulate a single scaling parameter
#  mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#  pixelWidth = 100 # store pixelWidth or grid resolution
#  
#  # create state-space grid and extent
#  Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*max(mysigma), res = pixelWidth)
#  
#  # create polygon for mask
#  library(sf)
#  poly = st_sfc(st_polygon(x=list(matrix(c(-1765,-1765,1730,-1650,1600,1650,0,1350,-800,1700,
#  -1850,1000,-1765,-1765),ncol=2, byrow=TRUE))), crs =  mycrs)
#  
#  # create habitat mask
#  hab_mask = mask_polygon(poly = poly, grid = Grid$grid, crs_ = mycrs,
#                          prev_mask = NULL)
#  
#  # simulate data for uniform state-space and habitat mask
#  data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs,
#              sigma_ = mysigma, prop_sex = 1,N = 200, K = 4,
#              base_encounter = 0.3, enc_dist = "poisson",
#              hab_mask = hab_mask, setSeed = 100)
#  
#  
#  # total augmented population size
#  M = 500
#  
#  # get initial activity center starting values
#  s.st3d = initialize_classic(y=NULL, M=M, X=traps, buff = 3*max(mysigma),
#                              hab_mask = hab_mask, all_random=TRUE)
#  
#  # make ggplot
#  ggplot() + geom_point(data=as.data.frame(Grid$grid),aes(x=x,y=y),color="grey60",
#                        size=1.25) +
#      geom_point(data=as.data.frame(traps),aes(x=x,y=y),color="blue",size=3) +
#      geom_point(data=as.data.frame(s.st3d),aes(x=V1,y=V2),color = "orangered",size=2.5,alpha=0.5) +
#      geom_sf(data=poly, fill = NA) + coord_sf(datum=st_crs(mycrs)) +
#      theme_classic() + ylab("Northing") + xlab("Easting") +
#      scale_x_continuous(expand=c(0.025, 0.025)) +
#      scale_y_continuous(expand=c(0.025, 0.025)) +
#      theme(axis.text = element_text(size=12),axis.title = element_text(size=16))

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # rescale inputs
#  rescale_list = rescale_classic(X = traps, ext = Grid$ext, s.st = s.st3d,
#                                 hab_mask = hab_mask)
#  
#  # store rescaled extent
#  ext = rescale_list$ext
#  
#  # Prepare data by summing over traps and occasions and add to list
#  data = list(n = apply(data3d$y, c(2,3), sum))
#  
#  # add rescaled traps
#  data$X = rescale_list$X
#  
#  # prepare constants (note get density in activity center/100 m2 rather than activity centers/m2)
#  constants = list(M = M,J=dim(data3d$y)[2], K=dim(data3d$y)[3],
#  x_lower = ext[1], x_upper = ext[2], y_lower = ext[3], y_upper = ext[4],
#  lam0_upper = 1,sigma_upper = 1000,
#  A = (sum(hab_mask)*(pixelWidth/100)^2),pixelWidth=pixelWidth)
#  
#  # add hab_mask and OK for habitat check
#  data$hab_mask = hab_mask
#  data$OK = rep(1,constants$M)
#  
#  # get initial activity center starting values
#  s.st = rescale_list$s.st
#  
#  # define all initial values
#  inits = list(sigma = runif(1, 250, 350), su = s.st,psi=runif(1,0.4,0.6),
#  lam0 = runif(1, 0.05, 0.15),pOK=data$OK,zu=rbinom(constants$M,1,0.5))
#  
#  # parameters to monitor
#  params = c("sigma","psi","lam0","N","D","su","zu")
#  
#  # get spatial count model
#  sc_model = get_unmarked(occ_specific = FALSE,
#                           hab_mask=TRUE,trapsClustered=FALSE)
#  
#  # run model (note we set s_alias to "su" for spatial count model)
#  library(tictoc)
#  tic() # track time elapsed
#  out = run_classic(model = sc_model, data=data, constants=constants,
#  inits=inits, params = params,niter = 10000, nburnin=1000, thin=1, nchains=2, parallel=TRUE,
#  RNGseed = 500, s_alias="su")
#  toc()
#  #> 511.62 sec elapsed
#  
#  # summarize output
#  samples = do.call(rbind, out)
#  
#  par(mfrow=c(1,1))
#  hist(samples[,which(dimnames(out[[1]])[[2]]=="N")], xlab = "Abundance", xlim = c(0,500), main="")
#  abline(v=200, col="red") # add line for simulated abundance

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # summarize MCMC samples (exclude parameters and don't plot)
#  nimSummary(out, exclude_params = c("su","zu"), trace=FALSE)
#  #>       post.mean post.sd    q2.5     q50   q97.5 f0  n.eff  Rhat
#  #> D         0.234   0.121   0.040   0.221   0.453  1 33.817 1.121
#  #> N       250.204 129.023  43.000 237.000 485.000  1 33.817 1.121
#  #> lam0      0.385   0.238   0.065   0.357   0.906  1 68.708 1.020
#  #> psi       0.500   0.258   0.086   0.475   0.968  1 33.594 1.119
#  #> sigma   344.133 202.596 146.564 270.782 932.274  1 31.582 1.282

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  
#  # get from Github and load
#  #install_github("austinnam/modeltools",force=TRUE)
#  library(modeltools)
#  
#  # get 'alpha' and 'beta' parameters of Gamma distribution
#  gparam = estGammaParam(mu = mysigma, sigma = 30)
#  
#  # view prior distribution for sigma (scaling parameter)
#  hist(rgamma(100000,shape=gparam$alpha,rate=1/gparam$beta),main="",xlab="Gamma prior")

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # define new code to replace line 3 in 'sc_model'
#  library(nimble)
#  new_prior = nimbleCode({
#     sigma ~ dgamma(alpha, 1/beta) # note that dgamma takes rate (1/beta)
#  })
#  
#  # delete old prior on line 3 and replace with new prior
#  sc_model_inf = customize_model(model = sc_model, append_code = new_prior,
#                                 line_append=3,line_remove=3)
#  
#  # inspect model (not run)
#  # sc_model_inf
#  
#  # add 'alpha' and 'beta' to list of constants
#  constants$alpha = gparam$alpha
#  constants$beta = gparam$beta
#  
#  # run model (note we set s_alias to "su" for spatial count model)
#  library(tictoc)
#  tic() # track time elapsed
#  out = run_classic(model = sc_model_inf, data=data, constants=constants,
#  inits=inits, params = params,niter = 10000, nburnin=1000, thin=1, nchains=2, parallel=TRUE,
#  RNGseed = 500,s_alias="su")
#  toc()
#  #> 836.56 sec elapsed
#  
#  # summarize output
#  samples = do.call(rbind, out)
#  par(mfrow=c(1,1))
#  hist(samples[,which(dimnames(out[[1]])[[2]]=="N")], xlab = "Abundance", xlim = c(0,500), main="")
#  abline(v=200, col="red") # add line for simulated abundance

## ---- fig.show='hide',eval=FALSE,fig.width = 8,fig.height=14------------------
#  # define new code to replace line 3 in 'sc_model'
#  new_prior = nimbleCode({
#     lam0 ~ dgamma(alpha_lam0, 1/beta_lam0) # note that dgamma takes rate (1/beta)
#  })
#  
#  # delete old prior on line 3 and replace with new prior
#  sc_model_inf2 = customize_model(model = sc_model_inf, append_code = new_prior,
#                                 line_append=2,line_remove=2)
#  
#  # inspect model (not run)
#  # sc_model_inf2
#  
#  # get 'alpha' and 'beta' parameters of Gamma distribution
#  lam0_param = estGammaParam(mu = 0.30, sigma = 0.03)
#  
#  # add 'alpha' and 'beta' to list of constants
#  constants$alpha_lam0 = lam0_param$alpha
#  constants$beta_lam0 = lam0_param$beta
#  
#  # run model (note we set s_alias to "su" for spatial count model)
#  library(tictoc)
#  tic() # track time elapsed
#  out = run_classic(model = sc_model_inf2, data=data, constants=constants,
#  inits=inits, params = params,niter = 10000, nburnin=1000, thin=1, nchains=2, parallel=TRUE,
#  RNGseed = 500,s_alias="su")
#  toc()
#  #> 1193.95 sec elapsed
#  
#  # summarize output (exclude "su" and "zu" from table and make posterior/trace plots)
#  
#  nimSummary(out, exclude_params = c("su","zu"), trace=TRUE, plot_all=FALSE)
#  #>       post.mean post.sd    q2.5     q50   q97.5 f0    n.eff  Rhat
#  #> D         0.223   0.053   0.137   0.219   0.343  1  275.831 1.001
#  #> N       238.946  56.449 147.000 234.000 367.000  1  275.831 1.001
#  #> lam0      0.298   0.029   0.244   0.297   0.358  1 1374.996 1.000
#  #> psi       0.478   0.114   0.288   0.468   0.737  1  287.766 1.001
#  #> sigma   287.307  28.393 234.682 286.116 346.612  1  301.998 1.013

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # make realized density plot (need to specify s_alias and z_alias)
#  r = realized_density(samples = out, grid = Grid$grid, crs_ = mycrs, site = NULL,
#                       hab_mask = hab_mask, s_alias = "su", z_alias = "zu")
#  
#  # load virdiis color palette and raster libraries
#  library(viridis)
#  library(raster)
#  
#  # make simple raster plot
#  plot(r, col=viridis(100),
#       main=expression("Realized density (activity centers/100 m"^2*")"),
#       ylab="Northing",xlab="Easting")

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # simulate a single trap array with random positional noise
#  x <- seq(-800, 800, length.out = 5)
#  y <- seq(-800, 800, length.out = 5)
#  traps <- as.matrix(expand.grid(x = x, y = y))
#  traps <- traps + runif(prod(dim(traps)),-20,20)
#  
#  mysigma = 300 # simulate single scaling parameter
#  mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#  pixelWidth = 100 # store pixelWidth
#  
#  # create an array of traps, as an approach where individuals will only be detected
#  # at one of the trap arrays (e.g., Furnas et al. 2018)
#  Xarray = array(NA, dim=c(nrow(traps),2,2))
#  Xarray[,,1]=traps
#  Xarray[,,2]=traps+4000 # shift trapping grid to new locations
#  
#  # create grid and extent for 3D trap array
#  GridX = grid_classic(X = Xarray, crs_ = mycrs, buff = 3*max(mysigma), res = 100)
#  
#  # create polygon to use as a mask
#  library(sf)
#  poly = st_sfc(st_polygon(x=list(matrix(c(-1660,-1900,5730,-1050,5470,
#  5650,0,6050,-1800,5700,-1660,-1900),ncol=2, byrow=TRUE))), crs =  mycrs)
#  
#  # make ggplot
#  ggplot() + geom_point(data=as.data.frame(GridX$grid[,,1]),aes(x=V1,y=V2),
#                        color="grey60",size=1.25) +
#      geom_point(data=as.data.frame(Xarray[,,1]),aes(x=V1,y=V2),color="blue",size=2) +
#      geom_point(data=as.data.frame(GridX$grid[,,2]),aes(x=V1,y=V2),color="grey60",
#                 size=1.25) +
#      geom_point(data=as.data.frame(Xarray[,,2]),aes(x=V1,y=V2),color="blue",size=2) +
#      geom_sf(data=poly, fill = NA) + coord_sf(datum=st_crs(mycrs)) +
#      theme_classic() + ylab("Northing") + xlab("Easting") +
#      scale_x_continuous(limits=c(-2000,6000)) +
#      scale_y_continuous(limits=c(-2000,6000)) +
#      theme(axis.text = element_text(size=12),axis.title = element_text(size=16))

## ---- fig.show='hide',eval=FALSE,fig.width = 8,fig.height=14------------------
#  # get 3D habitat mask array for 3D grid
#  hab_mask = mask_polygon(poly = poly, grid = GridX$grid, crs_ = mycrs, prev_mask = NULL)
#  
#  # simulate data for uniform state-space and habitat mask (N is simulated abundance)
#  data4d = sim_classic(X = Xarray, ext = GridX$ext, crs_ = mycrs,
#                       sigma_ = mysigma,prop_sex = 1,N = 200, K = 4,
#                     base_encounter = 0.3, enc_dist = "poisson",
#                     hab_mask = hab_mask, setSeed = 100)
#  
#  # organize by site and bind into an array
#  library(abind) # load abind package
#  y = abind(data4d$y[which(data4d$site==1),,],
#            data4d$y[which(data4d$site==2),,], along = 4)
#  
#  # total augmented population size
#  M = 400
#  
#  # augment site identifier
#  site = c(rep(1,200),rep(2,200))
#  
#  # get initial activity center starting values
#  s.st4d = initialize_classic(y=NULL, M=M, X=Xarray, buff = 3*max(mysigma),
#                        site = site, hab_mask = hab_mask,all_random = TRUE)
#  
#  # rescale inputs
#  rescale_list = rescale_classic(X = Xarray, ext = GridX$ext, s.st = s.st4d,
#                                 site = site, hab_mask = hab_mask)
#  
#  # store rescaled extent and convert to matrix
#  ext = do.call(rbind, lapply(rescale_list$ext, as.vector))
#  
#  # Prepare data by summing over traps and occasions and add to list
#  data = list(n = apply(y, c(2,3,4), sum),x_lower = ext[,1],
#              x_upper = ext[,2],y_lower = ext[,3]
#              ,y_upper = ext[,4],X = rescale_list$X)
#  
#  # add hab_mask, proportion of available habitat, and OK for habitat check
#  data$hab_mask = hab_mask
#  # need to adjust proportion of habitat available
#  data$prop.habitat=apply(hab_mask,3,mean)
#  data$OK = rep(1,constants$M)
#  
#  # prepare constants (note get density in activity center/100 m2)
#  constants = list(M = M,J=dim(data4d$y)[2],
#   K=dim(data4d$y)[3], sigma_upper = 1000, A = (sum(hab_mask)*(pixelWidth/100)^2),
#  pixelWidth=pixelWidth,nSites=dim(Xarray)[3],site = site)
#  
#  # add indexes for sites and individuals
#  constants$site_indexL = seq(1,M,200)
#  constants$site_indexU = seq(200,M,200)
#  
#  # priors for sigma: 'alpha' and 'beta'
#  constants$alpha = gparam$alpha
#  constants$beta = gparam$beta
#  
#  # priors for lam0: 'alpha' and 'beta'
#  constants$alpha_lam0 = lam0_param$alpha
#  constants$beta_lam0 = lam0_param$beta
#  
#  # get initial activity center starting values
#  s.st = rescale_list$s.st
#  
#  # define all initial values
#  inits = list(sigma = runif(1, 250, 350), su = s.st,psi=runif(1,0.4,0.6),
#  lam0 = runif(1, 0.1, 0.3),pOK=data$OK,zu=rbinom(constants$M,1,0.5))
#  
#  # get initial activity center starting values
#  s.st = rescale_list$s.st
#  
#  # define all initial values
#  inits = list(sigma = runif(1, 250, 350), su = s.st, psi=runif(1,0.2,0.3),
#  lam0 = runif(dim(data$X)[3], 0.1, 0.2),
#  pOK=data$OK,z=rbinom(M,1,0.5))
#  
#  # parameters to monitor
#  params = c("sigma","psi","lam0","N","D","su","zu")
#  
#  # get model
#  sc_model = get_unmarked(occ_specific = FALSE, hab_mask = TRUE,
#                          trapsClustered = TRUE)
#  
#  # model code to replace old code
#  add_model = nimbleCode({
#    lam0[g] ~ dgamma(alpha_lam0,1/beta_lam0)
#    sigma ~ dgamma(alpha, 1/beta)
#  })
#  
#  # now create new model
#  sc_model_inf = customize_model(sc_model, add_model, line_append = c(3,4),
#                                 line_remove = c(3,5))
#  
#  # inspect model (not run)
#  # sc_model_inf
#  
#  # run model (need to set s_alias)
#  library(tictoc)
#  tic() # track time elapsed
#  out = run_classic(model = sc_model_inf, data=data, constants=constants,
#  inits=inits, params = params,niter = 10000, nburnin=1000, thin=1, nchains=2,
#  parallel=TRUE, RNGseed = 500, s_alias = "su")
#  toc()
#  #> 968.93 sec elapsed
#  
#  # summary table of MCMC output (exclude "su" and "zu" parameters)
#  nimSummary(out, exclude_params = c("su","zu"))
#  #>         post.mean post.sd    q2.5     q50   q97.5 f0    n.eff  Rhat
#  #> D           0.102   0.022   0.066   0.099   0.152  1  270.831 1.008
#  #> N         228.569  49.067 149.000 223.000 342.000  1  270.831 1.008
#  #> lam0[1]     0.295   0.027   0.247   0.294   0.352  1 3595.823 1.002
#  #> lam0[2]     0.299   0.027   0.249   0.298   0.355  1 2805.852 1.001
#  #> psi         0.581   0.125   0.376   0.569   0.865  1  305.124 1.007
#  #> sigma     295.497  27.822 243.024 295.005 352.028  1  319.216 1.002

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # generate realized density surface (note setting z_alias and s_alias)
#  r = realized_density(samples=out, grid=GridX$grid, crs_=mycrs,
#                        site=constants$site, hab_mask=hab_mask,
#                         z_alias = "zu", s_alias = "su")
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
#            limits=c(0,0.3),breaks=seq(0,0.3,by=0.1)) +
#            xlab("") + ylab("") + theme_classic() +
#            scale_x_continuous(expand=c(0, 0)) +
#            scale_y_continuous(expand=c(0, 0)) +
#             theme(axis.text = element_text(size=18))
#  
#  # plot raster from site 2
#  p2<-gplot(r[[2]]) + geom_raster(aes(fill = value)) +
#            scale_fill_viridis(na.value = NA, name="Density",
#            limits=c(0,0.3),breaks=seq(0,0.3,by=0.1)) +
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

