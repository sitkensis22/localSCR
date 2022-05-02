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
#  mysigma = 300 # simulate sigma of 300 m
#  mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#  
#  # create state-space
#  Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = 100)
#  
#  # make ggplot of grid and trap locations
#  library(ggplot2)
#  ggplot() + geom_point(data=as.data.frame(Grid$grid),aes(x=x,y=y),color="grey60",
#                        size=1.25) +
#      geom_point(data=as.data.frame(traps),aes(x=x,y=y),color="blue",size=2) +
#      theme_classic() + ylab("Northing") + xlab("Easting") +
#      scale_x_continuous(expand=c(-0.1, 0.1)) +
#      scale_y_continuous(expand=c(-0.1, 0.1)) +
#      theme(axis.text = element_text(size=12),axis.title = element_text(size=16))

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # simulate SCR data
#  data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs,
#                sigma_ = mysigma, prop_sex = 1,N = 200, K = 4,
#                base_encounter = 0.10, enc_dist = "binomial",
#                hab_mask = FALSE, setSeed = 100)
#  
#  # inspect simulated data
#  str(data3d)
#  #> List of 3
#  #>  $ y  : int [1:200, 1:25, 1:4] 0 0 0 0 0 0 0 0 0 0 ...
#  #>  $ sex: int [1:200] 1 1 1 1 1 1 1 1 1 1 ...
#  #>  $ s  : num [1:200, 1:2] -658 -830 182 -1521 -106 ...
#  #>   ..- attr(*, "dimnames")=List of 2
#  #>   .. ..$ : NULL
#  #>   .. ..$ : chr [1:2] "sx" "sy"
#  
#  # bind simulated activity centers with vector of the number of individual detections
#  sdata = as.data.frame(cbind(data3d$s,Ndet=apply(data3d$y,1,sum)))
#  
#  # make ggplot
#  ggplot() + geom_point(data=as.data.frame(Grid$grid),aes(x=x,y=y),color="grey60",
#                        size=1.25) +
#      geom_point(data=as.data.frame(traps),aes(x=x,y=y),color="blue",size=2) +
#      geom_point(data=sdata,aes(x=sx,y=sy,size=Ndet),color="orangered",alpha=0.75) +
#      theme_classic() + ylab("Northing") + xlab("Easting") +
#      scale_x_continuous(expand=c(0.025, 0.025)) +  scale_size_continuous(range = c(2, 7)) +
#      scale_y_continuous(expand=c(0.025, 0.025)) +
#      theme(axis.text = element_text(size=12),axis.title = element_text(size=16))

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # simulate a single trap array with random positional noise
#  x <- seq(-800, 800, length.out = 5)
#  y <- seq(-800, 800, length.out = 5)
#  traps <- as.matrix(expand.grid(x = x, y = y))
#  set.seed(200)
#  traps <- traps + runif(prod(dim(traps)),-20,20)
#  
#  mysigma = c(220, 300) # simulate sex-specific
#  mycrs = 32608 # EPSG for WGS 84 / UTM zone 8N
#  pixelWidth = 100 # store pixelWidth or grid resolution
#  
#  # create state-space grid and extent
#  Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*max(mysigma), res = pixelWidth)
#  # create polygon for mask
#  library(sf)
#  poly = st_sfc(st_polygon(x=list(matrix(c(-1765,-1765,1730,-1650,1600,1650,0,1350,
#              -800,1700,-1850,1000,-1765,-1765),ncol=2, byrow=TRUE))), crs =  mycrs)
#  
#  # create habitat mask
#  hab_mask = mask_polygon(poly = poly, grid = Grid$grid, crs_ = mycrs,
#  prev_mask = NULL)
#  
#  # simulate data for uniform state-space and habitat mask
#  data3d = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs, sigma_ = mysigma, prop_sex = 0.7,
#  N = 200, K = 4, base_encounter = 0.15, enc_dist = "binomial",
#  hab_mask = hab_mask, setSeed = 100)
#  
#  # total augmented population size
#  M = 400
#  
#  # get initial activity center starting values
#  s.st3d = initialize_classic(y=data3d$y, M=M, X=traps, buff = ext=Grid$ext,
#                              hab_mask = hab_mask)
#  
#  # bind simulated activity centers with vector of either 1 (detected) or 0 (not detected)
#  s.stdata = as.data.frame(cbind(s.st3d,Det=c(as.numeric(apply(data3d$y,1,sum)>0),
#  rep(0,M-nrow(data3d$y)))))
#  
#  # convert Det to factor
#  s.stdata$Det = as.factor(s.stdata$Det)
#  
#  # make 1 as the baseline level for plotting
#  s.stdata$Det = factor(s.stdata$Det, levels = c("1","0"))
#  
#  # make ggplot
#  ggplot() + geom_point(data=as.data.frame(Grid$grid),aes(x=x,y=y),color="grey60",
#                        size=1.5) +
#      geom_point(data=as.data.frame(traps),aes(x=x,y=y),color="blue",size=3) +
#      geom_point(data=s.stdata,aes(x=V1,y=V2,fill=Det,color=Det),size = 3.5) +
#      geom_point(data=s.stdata[s.stdata$Det==1,],aes(x=V1,y=V2),size = 3.5,color="lawngreen") +
#      geom_sf(data=poly, fill = NA) + coord_sf(datum=st_crs(mycrs)) +
#      theme_classic() + ylab("Northing") + xlab("Easting") +
#      scale_x_continuous(expand=c(0.025, 0.025)) +
#      scale_color_manual(values = c("lawngreen","orangered")) +
#      scale_y_continuous(expand=c(0.025, 0.025)) +
#      theme(axis.text = element_text(size=12),axis.title = element_text(size=16))

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # rescale inputs
#  rescale_list = rescale_classic(X = traps, ext = Grid$ext, s.st = s.st3d,
#  hab_mask = hab_mask)
#  
#  # store rescaled extent
#  ext = rescale_list$ext
#  
#  # prepare data
#  data = list(y=data3d$y)
#  data$y = data$y[which(apply(data$y, 1, sum)!=0),,] # remove augmented records
#  data$y = apply(data$y, c(1,2), sum) # covert to 2d by summing over individuals and traps
#  
#  # add rescaled traps
#  data$X = rescale_list$X
#  
#  # prepare constants (note get density in activity center/100 m2 rather than activity centers/m2)
#  constants = list(M = M,n0 = nrow(data$y),J=dim(data$y)[2], K=dim(data3d$y)[3],
#  x_lower = ext[1], x_upper = ext[2], y_lower = ext[3], y_upper = ext[4],
#  sigma_upper = 1000, A = (sum(hab_mask)*(pixelWidth/100)^2),pixelWidth=pixelWidth)
#  
#  # augment sex
#  data$sex = c(data3d$sex,rep(NA,constants$M-length(data3d$sex)))
#  
#  # add z and zeros vector data for latent inclusion indicator
#  data$z = c(rep(1,constants$n0),rep(NA,constants$M - constants$n0))
#  data$zeros =  c(rep(NA,constants$n0),rep(0,constants$M - constants$n0))
#  
#  # add hab_mask and OK for habitat check
#  data$hab_mask = hab_mask
#  data$OK = rep(1,constants$M)
#  
#  # get initial activity center starting values
#  s.st3d = rescale_list$s.st
#  
#  # define all initial values
#  inits = list(sigma = runif(2, 250, 350), s = s.st3d,psi=runif(1,0.2,0.3),
#  p0 = runif(1, 0.05, 0.15),pOK=data$OK,z=c(rep(NA,constants$n0),
#  rep(0,constants$M-constants$n0)))
#  
#  # parameters to monitor
#  params = c("sigma","psi","p0","N","D","psi_sex","s","z")
#  
#  # get model
#  scr_model = get_classic(dim_y = 2, enc_dist = "binomial",sex_sigma = TRUE,hab_mask=TRUE,
#  trapsClustered=FALSE)
#  
#  # run model
#  library(tictoc)
#  tic() # track time elapsed
#  out = run_classic(model = scr_model, data=data, constants=constants,
#  inits=inits, params = params,niter = 10000, nburnin=1000, thin=1, nchains=2, parallel=TRUE,
#  RNGseed = 500)
#  toc()
#  #> 106.658 sec elapsed
#  
#  # summarize output
#  samples = do.call(rbind, out)
#  par(mfrow=c(1,1))
#  hist(samples[,which(dimnames(out[[1]])[[2]]=="N")], xlab = "Abundance", xlim = c(0,500), main="")
#  abline(v=200, col="red") # add line for simulated abundance

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # summarize MCMC samples (exclude parameters and don't plot)
#  nimSummary(out, exclude = c("s","z"), trace=FALSE)
#  #>          post.mean post.sd    q2.5     q50   q97.5 f0   n.eff  Rhat
#  #> D            0.217   0.027   0.174   0.214   0.279  1 287.723 1.003
#  #> N          232.708  28.761 186.000 229.000 298.000  1 287.723 1.003
#  #> p0           0.133   0.024   0.091   0.131   0.183  1 405.761 1.012
#  #> psi          0.581   0.076   0.456   0.574   0.749  1 298.192 1.003
#  #> psi_sex      0.758   0.054   0.648   0.761   0.855  1 465.875 1.003
#  #> sigma[1]   306.583  45.277 236.604 299.853 411.239  1 254.652 1.000
#  #> sigma[2]   278.029  25.147 232.727 276.585 331.032  1 273.259 1.011
#  
#  # make realized density plot
#  r = realized_density(samples = out, grid = Grid$grid, crs_ = mycrs,
#                       site = NULL, hab_mask = hab_mask)
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
#  set.seed(200)
#  traps <- traps + runif(prod(dim(traps)),-20,20) # add some random noise to locations
#  
#  mysigma = c(220, 300) # simulate sex-specific
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
#  ggplot() + geom_point(data=as.data.frame(GridX$grid[,,1]),aes(x=V1,y=V2),color="grey60",size=1.25) +
#      geom_point(data=as.data.frame(Xarray[,,1]),aes(x=V1,y=V2),color="blue",size=2) +
#      geom_point(data=as.data.frame(GridX$grid[,,2]),aes(x=V1,y=V2),color="grey60",size=1.25) +
#      geom_point(data=as.data.frame(Xarray[,,2]),aes(x=V1,y=V2),color="blue",size=2) +
#      geom_sf(data=poly, fill = NA) + coord_sf(datum=st_crs(mycrs)) +
#      theme_classic() + ylab("Northing") + xlab("Easting") +
#      scale_x_continuous(limits=c(-2000,6000)) +
#      scale_y_continuous(limits=c(-2000,6000)) +
#      theme(axis.text = element_text(size=12),axis.title = element_text(size=16))

## ---- fig.show='hide',fig.width = 8,fig.height=14,eval=FALSE------------------
#  # get 3D habitat mask array for 3D grid
#  hab_mask = mask_polygon(poly = poly, grid = GridX$grid, crs_ = mycrs,
#  prev_mask = NULL)
#  
#  # simulate data for uniform state-space and habitat mask (N is simulated abundance)
#  data4d = sim_classic(X = Xarray, ext = GridX$ext, crs_ = mycrs, sigma_ = mysigma,
#                  prop_sex = 0.7,N = 200, K = 4, base_encounter = 0.15,
#                  enc_dist = "binomial",hab_mask = hab_mask, setSeed = 500)
#  
#  # total augmented population size
#  M = 400
#  
#  # augment site identifier
#  site = c(data4d$site,c(rep(1,((M-length(data4d$site))/2)),
#                         rep(2,((M-length(data4d$site))/2))))
#  
#  # get initial activity center starting values
#  s.st4d = initialize_classic(y=data4d$y, M=M, X=Xarray, ext = GridX$ext,
#                              site = site, hab_mask = hab_mask)
#  
#  # rescale inputs
#  rescale_list = rescale_classic(X = Xarray, ext = GridX$ext, s.st = s.st4d,
#                                 site = site, hab_mask = hab_mask)
#  
#  # store rescaled extent and convert to matrix
#  ext = do.call(rbind, lapply(rescale_list$ext, as.vector))
#  
#  # prepare constants (note get density in activity center/100 m2)
#  constants = list(M = M,n0 =  length(which(apply(data4d$y,1,sum)!=0)),
#              J=dim(data4d$y)[2], K=dim(data4d$y)[3], sigma_upper = 1000,
#              A = (sum(hab_mask)*(pixelWidth/100)^2),
#              pixelWidth=pixelWidth,nSites=dim(Xarray)[3],site = site)
#  
#  # prepare data
#  data = list(X = rescale_list$X,sex = c(data4d$sex,rep(NA,M-length(data4d$sex))),
#  x_lower = ext[,1],x_upper = ext[,2],y_lower = ext[,3],y_upper = ext[,4])
#  
#  # store and format encounter history data
#  data$y = data4d$y[which(apply(data4d$y, 1, sum)!=0),,] # remove augmented records
#  # covert to 2d by summing over individuals and traps
#  data$y = apply(data$y, c(1,2), sum)
#  
#  # add z and zeros vector data for latent inclusion indicator
#  data$z = c(rep(1,constants$n0),rep(NA,constants$M - constants$n0))
#  data$zeros =  c(rep(NA,constants$n0),rep(0,constants$M - constants$n0))
#  
#  # add hab_mask, proportion of available habitat, and OK for habitat check
#  data$hab_mask = hab_mask
#  # need to adjust proportion of habitat available
#  data$prop.habitat=apply(hab_mask,3,mean)
#  data$OK = rep(1,constants$M)
#  
#  # get initial activity center starting values
#  s.st = rescale_list$s.st
#  
#  # define all initial values
#  inits = list(sigma = runif(2, 250, 350), s = s.st,psi=runif(1,0.2,0.3),
#  p0 = runif(dim(data$X)[3], 0.1, 0.2),sex=ifelse(is.na(data$sex),
#  rbinom(constants$M-constants$n0,1,0.5),NA),
#  pOK=data$OK,z=c(rep(NA,constants$n0),rep(0,constants$M-constants$n0)),
#  psi_sex=runif(1,0.4,0.6))
#  
#  # parameters to monitor
#  params = c("sigma","psi","p0","N","D","psi_sex","s","z")
#  
#  # get model
#  scr_model = get_classic(dim_y = 2, enc_dist = "binomial",
#                  sex_sigma = TRUE,hab_mask=TRUE,trapsClustered = TRUE)
#  
#  # run model
#  library(tictoc)
#  tic() # track time elapsed
#  out = run_classic(model = scr_model, data=data, constants=constants,
#  inits=inits, params = params,niter = 10000, nburnin=1000, thin=1, nchains=2,
#  parallel=TRUE, RNGseed = 500)
#  toc()
#  #> 110.978 sec elapsed
#  
#  # summary table of MCMC output (exclude "s" and "z" parameters)
#  nimSummary(out, exclude = c("s","z"), trace = TRUE, plot_all = FALSE)
#  #>          post.mean post.sd    q2.5     q50   q97.5 f0   n.eff  Rhat
#  #> D            0.100   0.015   0.074   0.098   0.136  1 221.187 1.063
#  #> N          223.906  34.525 166.000 221.000 306.000  1 221.187 1.063
#  #> p0[1]        0.112   0.027   0.066   0.109   0.173  1 531.737 1.005
#  #> p0[2]        0.108   0.024   0.069   0.106   0.162  1 534.815 1.006
#  #> psi          0.569   0.090   0.409   0.562   0.776  1 229.467 1.062
#  #> psi_sex      0.553   0.079   0.394   0.553   0.703  1 304.397 1.021
#  #> sigma[1]   252.076  37.377 192.030 247.679 335.659  1 264.634 1.009
#  #> sigma[2]   330.496  35.807 270.581 326.672 411.406  1 268.864 1.015

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # generate realized density surface
#  r = realized_density(samples=out, grid=GridX$grid, crs_=mycrs,
#                        site=constants$site, hab_mask=hab_mask)
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

