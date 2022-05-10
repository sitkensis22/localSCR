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
#  x <- seq(-1600, 1600, length.out = 6)
#  y <- seq(-1600, 1600, length.out = 6)
#  traps <- as.matrix(expand.grid(x = x, y = y))
#  # add some random noise to locations
#  set.seed(100)
#  traps <- traps + runif(prod(dim(traps)),-20,20)
#  mysigma = 300 # simulate sigma of 300 m
#  mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#  pixelWidth = 100 # grid resolution
#  
#  # Simulated abundance
#  Nsim = 250
#  
#  # manually create state-space grid and extent (note that I ran the code through
#  # the localize_classic() using the same settings as other vignettes for
#  # grid_classic() without using a habitat mask to get the extent
#  # for the scaled-up state-space grid)
#  
#  # make raster layer of grid
#  r_grid = raster::raster(xmn=-2698.706,xmx=2701.294,ymn=-2700.908,ymx=2699.092,
#                          res = pixelWidth, crs = mycrs)
#  Grid = list() # create Grid list
#  Grid$grid = raster::coordinates(r_grid)
#  Grid$ext = raster::extent(r_grid)
#  
#  # create polygon to use as a mask
#  library(sf)
#  poly = st_sfc(st_polygon(x=list(matrix(c(-2465,-2465,2530,-2550,2650,2550,
#  0,2550,-800,2500,-2350,2300,-2465,-2465),ncol=2, byrow=TRUE))), crs =  mycrs)
#  
#  # make simple plot
#  par(mfrow=c(1,1))
#  plot(Grid$grid, pch=20, col="gray60")
#  points(traps, col="blue",pch=20)
#  plot(poly, add=TRUE)

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # create habitat mask from polygon
#  hab_mask = mask_polygon(poly = poly, grid = Grid$grid, crs_ = mycrs,
#  prev_mask = NULL)
#  
#  # simulate SCR data
#  data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, sigma_ = mysigma,
#                    prop_sex = 1, N = Nsim, K = 4, base_encounter = 0.15,
#                    enc_dist = "binomial", hab_mask = hab_mask, setSeed = 100)
#  
#  # augmented population size (detected + augmented individuals)
#  M=500
#  
#  # generate initial activity center coordinates for 2D trap array without
#  # habitat mask
#  s.st = initialize_classic(y=data3d$y, M=M, X=traps, ext = Grid$ext,
#  hab_mask = hab_mask)
#  
#  # now use grid_classic to create an individual-level state-space (with origin 0, 0)
#  Grid_ind = grid_classic(X = matrix(c(0,0),nrow=1), crs_ = mycrs, buff = 3*mysigma, res = 100)
#  
#  # now localize the data components created above (a bit time consuming ~ 20 sec)
#  # set layers to equal to evenly augment state-space
#  library(tictoc)
#  tic()
#  local_list = localize_classic(y = data3d$y, grid_ind = Grid_ind$grid, X=traps,
#                              crs_ = mycrs, sigma_ = mysigma, s.st = s.st,
#                              hab_mask = hab_mask)
#  toc()
#  #> 40.25 sec elapsed
#  
#  # inspect local_list
#  str(local_list)
#  #> List of 8
#  #>  $ y           : int [1:99, 1:36, 1:4] 0 0 0 0 0 0 0 0 0 0 ...
#  #>  $ X           : num [1:500, 1:36, 1:2] -1608 -965 -970 965 -1588 ...
#  #>  $ grid        : num [1:2916, 1:2] -2649 -2549 -2449 -2349 -2249 ...
#  #>   ..- attr(*, "dimnames")=List of 2
#  #>   .. ..$ : NULL
#  #>   .. ..$ : chr [1:2] "x" "y"
#  #>  $ prop_habitat: num [1:500] 1 1 1 1 0.981 ...
#  #>  $ ext_mat     : num [1:500, 1:4] -573 -1226 -572 716 -1843 ...
#  #>  $ ext         :Formal class 'Extent' [package "raster"] with 4 slots
#  #>   .. ..@ xmin: num -2699
#  #>   .. ..@ xmax: num 2701
#  #>   .. ..@ ymin: num -2701
#  #>   .. ..@ ymax: num 2699
#  #>  $ Jind        : num [1:500] 35 26 35 19 24 24 24 26 35 25 ...
#  #>  $ s.st        : num [1:500, 1:2] 327 -326 328 1616 -943 ...

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # rescale inputs
#  rescale_list = rescale_local(X = local_list$X, ext = local_list$ext,
#                               ext_mat = local_list$ext_mat,
#                               s.st = local_list$s.st, hab_mask = hab_mask)
#  
#  # prepare encounter data data
#  data = list(y=local_list$y)
#  data$y = apply(data$y, c(1,2), sum) # covert to 2d by summing over individuals and traps
#  
#  # add rescaled traps
#  data$X = rescale_list$X
#  
#  # prepare constants (rescale area to activity centers/km 2)
#  constants = list(M = M,n0 = nrow(data$y),Jind = local_list$Jind, K = 4,
#                   xy_bounds = rescale_list$ext_mat, sigma_upper = 1000,
#                   pixelWidth=pixelWidth,
#                   A = prod(c(abs(diff(local_list$ext[1:2]))/1000,
#                              abs(diff(local_list$ext[3:4]))/1000)))
#  
#  # add z and zeros vector data for latent inclusion indicator
#  data$z = c(rep(1,constants$n0),rep(NA,constants$M - constants$n0))
#  data$zeros =  c(rep(NA,constants$n0),rep(0,constants$M - constants$n0))
#  
#  # add hab_mask,OK for habitat check, and proportion of available habitat
#  data$hab_mask = hab_mask
#  data$OK = rep(1,constants$M)
#  data$prop.habitat = local_list$prop_habitat
#  
#  # define all initial values
#  inits = list(sigma = runif(1, 250, 350), s = rescale_list$s.st,psi=runif(1,0.2,0.3),
#            p0 = runif(1, 0.05, 0.15), z=c(rep(NA,constants$n0),rep(0,constants$M-constants$n0)),
#            pOK=data$OK)
#  
#  # parameters to monitor
#  params = c("sigma","psi","p0","N","D","s","z")
#  
#  # get model
#  scr_model = get_classic(dim_y = 2, enc_dist = "binomial",sex_sigma = FALSE,
#  hab_mask=TRUE,trapsClustered = FALSE)
#  
#  # show lines of model to help with editing
#  print_model(scr_model)
#  
#  # add code to inject localization lines
#  add_model = nimble::nimbleCode({
#    z[i] ~ dbern(psim[i])
#    psim[i] <- (1 - (1 - psi)^prop.habitat[i])
#    s[i,1] ~ dunif(xy_bounds[i,1], xy_bounds[i,2])
#    s[i,2] ~ dunif(xy_bounds[i,3], xy_bounds[i,4])
#    dist[i, 1:Jind[i]] <- sqrt((s[i,1] - X[i,1:Jind[i],1])^2 + (s[i, 2] - X[i,1:Jind[i], 2])^2)
#    p[i,1:Jind[i]] <- p0 * exp(-dist[i, 1:Jind[i]]^2/(2 * sigma.pixel^2))
#    for (i in 1:n0) {
#       for (j in 1:Jind[i]) {
#           y[i, j] ~ dbin(p[i, j], K)
#       }
#    }
#    for (i in (n0 + 1):M) {
#       zeros[i] ~ dbern((1 - prod(1 - p[i, 1:Jind[i]])^K) * z[i])
#    }
#  })
#  
#  # edit model (injecting two new lines on line 7, otherwise 1-1 line replacement)
#  local_model = customize_model(model = scr_model,append_code = add_model,
#                                remove_line = c(7,8:9,12:13,15:22),
#                                append_line = c(7,7,8:9,12:13,15:22))
#  
#  # run model
#  library(tictoc)
#  tic() # track time elapsed
#  out = run_classic(model = local_model, data=data, constants=constants,
#        inits=inits, params = params, niter = 10000, nburnin=1000,
#        thin=1, nchains=2, parallel=TRUE, RNGseed = 500)
#  toc()
#  #> 125.08 sec elapsed
#  
#  # summarize output
#  samples = do.call(rbind, out)
#  par(mfrow=c(1,1))
#  hist(samples[,which(dimnames(out[[1]])[[2]]=="N")], xlab = "Abundance", xlim = c(0,500), main="")
#  abline(v=Nsim, col="red") # add line for simulated abundance

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # summarize MCMC samples (exclude parameters and don't plot)
#  nimSummary(out, exclude = c("s","z"), trace=FALSE)
#  #>       post.mean post.sd    q2.5     q50   q97.5 f0   n.eff  Rhat
#  #> D         8.901   1.099   7.064   8.779  11.351  1 346.433 1.013
#  #> N       259.560  32.057 206.000 256.000 331.000  1 346.433 1.013
#  #> p0        0.194   0.039   0.127   0.191   0.280  1 361.349 1.014
#  #> psi       0.567   0.069   0.446   0.561   0.716  1 330.496 1.014
#  #> sigma   272.652  22.015 233.584 271.202 322.137  1 294.354 1.004
#  
#  # make realized density plot
#  r = realized_density(samples = out, grid = local_list$grid,
#                       crs_ = mycrs, hab_mask = hab_mask)
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
#  x <- seq(-1600, 1600, length.out = 6)
#  y <- seq(-1600, 1600, length.out = 6)
#  traps <- as.matrix(expand.grid(x = x, y = y))
#  # add some random noise to locations
#  set.seed(100)
#  traps <- traps + runif(prod(dim(traps)),-20,20)
#  mysigma = 300 # simulate sigma of 300 m
#  mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#  pixelWidth = 100 # grid resolution
#  
#  # Simulated abundance
#  Nsim = 400
#  
#  # create an array of traps, as an approach where individuals will only be detected
#  # at one of the trap arrays (e.g., Furnas et al. 2018)
#  Xarray = array(NA, dim=c(nrow(traps),2,2))
#  Xarray[,,1]=traps
#  Xarray[,,2]=traps+10000 # shift trapping grid to new locations
#  
#  # create grid and extent for 3D trap array
#  GridX = grid_classic(X = Xarray, crs_ = mycrs, buff = 3*max(mysigma), res = 100)
#  
#  # create polygon to use as a mask
#  library(sf)
#  poly = st_sfc(st_polygon(x=list(matrix(c(-2160,-2900,14430,-1550,12270,
#  12400,0,13050,-2800,13100,-2160,-2900),ncol=2, byrow=TRUE))), crs =  mycrs)
#  
#  # make a simple plot
#  plot(Xarray[,,1],xlim=c(-3000,16000),ylim=c(-3000,16000),
#       xlab="Easting",ylab="Northing")
#  points(GridX$grid[,,2],col="gray60",pch=20)
#  points(GridX$grid[,,1],col="gray60",pch=20)
#  points(Xarray[,,1],pch=20,col="blue")
#  points(Xarray[,,2],pch=20,col="blue")
#  plot(poly, add=TRUE)

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # get 3D habitat mask array for 3D grid
#  hab_mask = mask_polygon(poly = poly, grid = GridX$grid, crs_ = mycrs,
#  prev_mask = NULL)
#  
#  # simulate data for uniform state-space and habitat mask (N is simulated abundance)
#  data4d = sim_classic(X = Xarray, ext = GridX$ext, crs_ = mycrs, sigma_ = mysigma,
#                  prop_sex = 1,N = Nsim, K = 4, base_encounter = 0.15,
#                  enc_dist = "binomial",hab_mask = hab_mask, setSeed = 500)
#  
#  M=800
#  
#  # augment site identifier
#  site = c(data4d$site,c(rep(1,((M-length(data4d$site))/2)),
#                         rep(2,((M-length(data4d$site))/2))))
#  
#  # get initial activity center starting values
#  s.st = initialize_classic(y=data4d$y, M=M, X=Xarray, ext=GridX$ext,
#                              site = site, hab_mask = hab_mask)
#  
#  # now use grid_classic to create an individual-level state-space (with origin 0, 0)
#  Grid_ind = grid_classic(X = matrix(c(0,0),nrow=1), crs_ = mycrs, buff = 3*mysigma, res = 100)
#  
#  # now localize the data components created above (a bit time consuming ~ 40 sec)
#  # set layers to equal to evenly augment state-space
#  library(tictoc)
#  tic()
#  local_list = localize_classic(y = data4d$y, grid_ind = Grid_ind$grid, X=Xarray,
#                              crs_ = mycrs, sigma_ = mysigma, s.st = s.st,
#                              site = site,hab_mask = hab_mask)
#  toc()
#  
#  # inspect local list
#  str(local_list)
#  #> List of 8
#  #>  $ y           : int [1:150, 1:36, 1:4] 0 0 0 0 0 0 0 0 0 0 ...
#  #>  $ X           : num [1:800, 1:36, 1:2] 965 -970 959 -1608 -1588 ...
#  #>  $ grid        : num [1:2916, 1:2, 1:2] -2649 -2549 -2449 -2349 -2249 ...
#  #>  $ prop_habitat: num [1:800] 1 1 1 1 1 ...
#  #>  $ ext_mat     : num [1:800, 1:4] 716 707 691 -593 -1843 ...
#  #>  $ ext         :List of 2
#  #>   ..$ :Formal class 'Extent' [package "raster"] with 4 slots
#  #>   .. .. ..@ xmin: num -2699
#  #>   .. .. ..@ xmax: num 2701
#  #>   .. .. ..@ ymin: num -2701
#  #>   .. .. ..@ ymax: num 2699
#  #>   ..$ :Formal class 'Extent' [package "raster"] with 4 slots
#  #>   .. .. ..@ xmin: num 7301
#  #>   .. .. ..@ xmax: num 12701
#  #>   .. .. ..@ ymin: num 7299
#  #>   .. .. ..@ ymax: num 12699
#  #>  $ Jind        : num [1:800] 19 21 24 32 24 26 30 25 29 24 ...
#  #>  $ s.st        : num [1:800, 1:2] 1616 1607 1591 307 -943 ...
#  
#  # rescale inputs
#  rescale_list = rescale_local(X = local_list$X, ext = local_list$ext, ext_mat = local_list$ext_mat,
#                               s.st = local_list$s.st, site = site, hab_mask = hab_mask)
#  
#  # prepare encounter data data
#  data = list(y=local_list$y)
#  data$y = apply(data$y, c(1,2), sum) # covert to 2d by summing over individuals and traps
#  
#  # add rescaled traps
#  data$X = rescale_list$X
#  
#  # prepare constants (rescale area to activity centers/km 2)
#  constants = list(M = M,n0 = nrow(data$y),Jind = local_list$Jind, K = 4,
#                   xy_bounds = rescale_list$ext_mat, sigma_upper = 1000,
#                   pixelWidth=pixelWidth,
#                   A = (sum(hab_mask)*(pixelWidth/1000)^2),
#                   nSites=dim(Xarray)[3],site = site)
#  
#  # add z and zeros vector data for latent inclusion indicator
#  data$z = c(rep(1,constants$n0),rep(NA,constants$M - constants$n0))
#  data$zeros =  c(rep(NA,constants$n0),rep(0,constants$M - constants$n0))
#  
#  # add hab_mask,OK for habitat check, and proportion of available habitat
#  data$hab_mask = hab_mask
#  data$OK = rep(1,constants$M)
#  data$prop.habitat = local_list$prop_habitat
#  
#  # define all initial values
#  inits = list(sigma = runif(1, 250, 350), s = rescale_list$s.st,psi=runif(1,0.2,0.3),
#            p0 = runif(constants$nSites, 0.05, 0.15),
#            z=c(rep(NA,constants$n0),rep(0,constants$M-constants$n0)),
#            pOK=data$OK)
#  
#  # parameters to monitor
#  params = c("sigma","psi","p0","N","D","s","z")
#  
#  # get model
#  scr_model = get_classic(dim_y = 2, enc_dist = "binomial",sex_sigma = FALSE,
#  hab_mask=TRUE,trapsClustered = TRUE)
#  
#  # show lines of model to help with editing
#  print_model(scr_model)
#  
#  # create code to inject localization lines
#  add_model = nimble::nimbleCode({
#    psim[i] <- (1 - (1 - psi)^prop.habitat[i])
#    s[i,1] ~ dunif(xy_bounds[i,1], xy_bounds[i,2])
#    s[i,2] ~ dunif(xy_bounds[i,3], xy_bounds[i,4])
#    dist[i, 1:Jind[i]] <- sqrt((s[i,1] - X[i,1:Jind[i],1])^2 + (s[i, 2] - X[i,1:Jind[i], 2])^2)
#    p[i,1:Jind[i]] <- p0[site[i]] * exp(-dist[i, 1:Jind[i]]^2/(2 * sigma.pixel^2))
#    for (i in 1:n0) {
#       for (j in 1:Jind[i]) {
#           y[i, j] ~ dbin(p[i, j], K)
#       }
#    }
#    for (i in (n0 + 1):M) {
#       zeros[i] ~ dbern((1 - prod(1 - p[i, 1:Jind[i]])^K) * z[i])
#    }
#  })
#  
#  # edit model (a straight 1-1 line replacement here)
#  local_model = customize_model(model = scr_model,append_code = add_model,
#                                remove_line = c(10:12,15:16,18:25),
#                                append_line = c(10:12,15:16,18:25))
#  
#  # run model
#  library(tictoc)
#  tic() # track time elapsed
#  out = run_classic(model = local_model, data=data, constants=constants,
#        inits=inits, params = params, niter = 10000, nburnin=1000,
#        thin=1, nchains=2, parallel=TRUE, RNGseed = 500)
#  toc()
#  #> 237.82 sec elapsed
#  
#  # summarize MCMC output
#  nimSummary(out, exclude = c("s","z"))
#  #>       post.mean post.sd    q2.5     q50   q97.5 f0   n.eff  Rhat
#  #> D         8.754   0.823   7.314   8.694  10.527  1 361.800 1.004
#  #> N       424.939  39.928 355.000 422.000 511.000  1 361.800 1.004
#  #> p0[1]     0.167   0.029   0.116   0.165   0.230  1 476.028 1.010
#  #> p0[2]     0.176   0.030   0.126   0.174   0.242  1 514.860 1.013
#  #> psi       0.540   0.053   0.447   0.537   0.650  1 348.219 1.004
#  #> sigma   276.489  16.936 244.493 276.152 310.937  1 328.285 1.010
#  
#  # plot posterior abundance with line for simulated abundance
#  samples = do.call(rbind, out)
#  par(mfrow=c(1,1))
#  hist(samples[,which(dimnames(out[[1]])[[2]]=="N")], xlab = "Abundance",
#       xlim = c(0,900), main="")
#  abline(v=Nsim, col="red") # add line for simulated abundance

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # make realized density plot (must trace "s" and "z" above)
#  r = realized_density(samples=out, grid=local_list$grid, ext = local_list$ext,
#                    crs_=mycrs, site=constants$site, hab_mask=hab_mask)
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

