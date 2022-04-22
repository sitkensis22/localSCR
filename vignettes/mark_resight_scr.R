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
#  pixelWidth=100 # define pixelWidth or grid resolution
#  
#  # create state-space
#  Grid = grid_classic(X = traps, crs_ = mycrs, buff = 3*mysigma, res = pixelWidth)
#  
#  # create polygon for habitat mask
#  library(sf)
#  poly = st_sfc(st_polygon(x=list(matrix(c(-1765,-1765,1730,-1650,1600,1650,0,1350,-800,1700,
#  -1850,1000,-1765,-1765),ncol=2, byrow=TRUE))), crs =  mycrs)
#  
#  # create habitat mask
#  hab_mask = mask_polygon(poly = poly, grid = Grid$grid, crs_ = mycrs,
#                          prev_mask = NULL)
#  
#  # make ggplot of grid/habitat mask and trap locations
#  library(ggplot2)
#  grid_data = data.frame(x=Grid$grid[,1],y=Grid$grid[,2],Suitable=as.factor(as.vector(t(apply(hab_mask,2,rev)))))
#  
#  # make 1 as the baseline level for plotting
#  grid_data$Suitable= factor(grid_data$Suitable, levels = c("1","0"))
#  
#  # ggplot
#  ggplot() + geom_point(data=grid_data,aes(x=x,y=y,fill=Suitable,color=Suitable),
#                        size=1.25) +
#      geom_point(data=as.data.frame(traps),aes(x=x,y=y),color="blue",size=2) +
#      theme_classic() + ylab("Northing") + xlab("Easting") +
#      geom_sf(data=poly, fill = NA) + coord_sf(datum=st_crs(mycrs)) +
#      scale_color_manual(values = c("orangered","gray60")) +
#      scale_x_continuous(expand=c(0.025, 0.025)) +
#      scale_y_continuous(expand=c(0.025, 0.025)) +
#      theme(axis.text = element_text(size=12),axis.title = element_text(size=16))

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # simulate marked data
#  data_m = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs,
#                       sigma_ = mysigma, prop_sex = 1,N = 100, K = 4,
#                       base_encounter = 0.15, enc_dist = "poisson",
#                       hab_mask = hab_mask, setSeed = 100)
#  
#  # inspect simulated data
#  str(data_m)
#  #> List of 3
#  #>  $ y  : int [1:100, 1:25, 1:4] 0 0 0 0 0 0 0 1 0 0 ...
#  #>  $ sex: int [1:100] 1 1 1 1 1 1 1 1 1 1 ...
#  #>  $ s  : num [1:100, 1:2] -657 -828 178 -1515 -108 ...
#  #>   ..- attr(*, "dimnames")=List of 2
#  #>   .. ..$ : NULL
#  #>   .. ..$ : chr [1:2] "sx" "sy"

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # simulate marked data (we use a different 'setSeed' to produce different data)
#  data_u = sim_classic(X = traps, ext = Grid$ext, crs_ = mycrs,
#                       sigma_ = mysigma, prop_sex = 1,N = 100, K = 4,
#                       base_encounter = 0.15, enc_dist = "poisson",
#                       hab_mask = hab_mask, setSeed = 500)
#  
#  # Covert to unmarked data
#  # sum over traps and occasions to produce a 2-dimensional spatial count data set
#  n = apply(data_u$y, c(2,3), sum)
#  
#  # inspect n[j,k]
#  str(n)
#  #> int [1:25, 1:4] 1 3 2 0 0 1 0 0 0 1 ...

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # we'll use the traps, state-space, habitat mask, and
#  # data sets from above rather than recreating them here
#  
#  # Prep marked data
#  # marked augmented population size
#  M = 200
#  
#  # get initial activity center starting values for marked individuals
#  s.st_m = initialize_classic(y=data_m$y, M=M, X=traps, buff = 3*max(mysigma),
#                              hab_mask = hab_mask)
#  
#  # rescale inputs
#  rescale_list_m = rescale_classic(X = traps, ext = Grid$ext, s.st = s.st_m,
#                                 hab_mask = hab_mask)
#  
#  # store rescaled starting activity center coordinates
#  s.st_m = rescale_list_m$s.st
#  
#  # store rescaled extent
#  ext = rescale_list_m$ext
#  
#  # Prep marked data
#  # unmarked augmented population size
#  m = 200
#  
#  # get initial activity center starting values for unmarked individuals
#  s.st_u = initialize_classic(y=NULL, M=m, X=traps, buff = 3*max(mysigma),
#                              hab_mask = hab_mask, all_random = TRUE)
#  
#  # rescale inputs (note we'll use the traps and extent from above since
#  # they are the same for all individuals, and we'll just grab s.st)
#  s.st_u = rescale_classic(X = traps, ext = Grid$ext, s.st = s.st_u,
#                                 hab_mask = hab_mask)$s.st
#  
#  # create constants list
#  constants = list(m=m, M = M,J=dim(data_m$y)[2], K=dim(data_m$y)[3],
#        n0 = length(which(apply(data_m$y,c(1),sum)!=0)),
#        x_lower = ext[1], x_upper = ext[2], y_lower = ext[3], y_upper = ext[4],
#        lam0_upper = 1,sigma_upper = 1000,
#        A = (sum(hab_mask)*(pixelWidth/100)^2),pixelWidth=pixelWidth)
#  
#  # remove rows with no detections from marked data y
#  y = data_m$y[which(apply(data_m$y,1,sum)!=0),,]
#  # convert y to 2D encounter history
#  y = apply(y, c(1,2), sum)
#  
#  # create combined data list
#  data = list(y = y,
#              n = n, X = rescale_list_m$X)
#  data$hab_mask = hab_mask
#  data$OK = rep(1,constants$M)
#  data$OKu = rep(1,constants$m)
#  
#  # add z and zeros vector data for latent inclusion indicator for marked individuals
#  data$z = c(rep(1,constants$n0),rep(NA,constants$M - constants$n0))
#  data$zeros =  c(rep(NA,constants$n0),rep(0,constants$M - constants$n0))
#  
#  # define all initial values
#  inits = list(sigma = runif(1, 250, 350),s=s.st_m,su = s.st_u,
#           psi=runif(1,0.4,0.6),psiu=runif(1,0.4,0.6),
#          lam0 = runif(1, 0.05, 0.15),pOK=data$OK,pOKu=data$OKu,
#          z=c(rep(NA,constants$n0),rep(0,constants$M-constants$n0)),
#          zu=rbinom(constants$m,1,0.5))
#  
#  # parameters to monitor
#  params = c("sigma","psi","psiu","lam0","N","Nm","Nu","D","s","z","su","zu")
#  
#  # get 'classic' SCR model for marked data
#  scr_model = get_classic(dim_y = 2, enc_dist = "poisson",
#                                 sex_sigma = FALSE,hab_mask=TRUE,
#                                  trapsClustered=FALSE)
#  
#  # get spatial count model for unmarked data
#  sc_model = get_unmarked(occ_specific = FALSE,
#                          hab_mask=TRUE,trapsClustered=FALSE)
#  
#  # now combine models to form spatial mark-resight
#  smr_model = customize_model(scr_model, sc_model)
#  
#  # now write new model code to update smr_model lines
#  library(nimble)
#  add_model = nimbleCode({
#    Nm <- sum(z[1:M])
#    Nu <- sum(zu[1:m])
#    N <- Nm + Nu
#  })
#  
#  # erase rows and insert add_model
#  # add new lines all to the same line (line 44 of the smr_model)
#  smr_model2 = customize_model(smr_model, append_code = add_model,
#              append_line=rep(44,3), remove_line = c(2:4,23:24,44))
#  
#  # inspect mark-resight model (not run)
#  # smr_model2
#  
#  library(tictoc)
#  tic() # track time elapsed
#  out = run_classic(model = smr_model2, data=data, constants=constants,
#  inits=inits, params = params,niter = 10000, nburnin=1000, thin=1, nchains=2, parallel=TRUE,
#  RNGseed = 500, s_alias=c("s","su"))
#  toc()
#  #> 175.6 sec elapsed
#  
#  nimSummary(out, exclude_params = c("s","z","su","zu"), trace=TRUE,
#             plot_all = FALSE)
#  #>       post.mean post.sd    q2.5     q50   q97.5 f0   n.eff  Rhat
#  #> D         0.206   0.041   0.141   0.202   0.300  1 123.795 1.057
#  #> N       220.077  43.200 150.000 215.000 320.000  1 123.795 1.057
#  #> Nm      110.097  20.580  77.000 108.000 158.000  1 217.391 1.018
#  #> Nu      109.980  27.535  65.000 106.000 172.000  1 156.183 1.078
#  #> lam0      0.143   0.039   0.081   0.139   0.221  1 236.030 1.056
#  #> psi       0.551   0.107   0.371   0.540   0.799  1 213.259 1.015
#  #> psiu      0.549   0.140   0.312   0.533   0.862  1 170.283 1.082
#  #> sigma   285.441  32.722 227.818 282.206 356.657  1 206.288 1.002

## ---- fig.show='hide',eval=FALSE----------------------------------------------
#  # generate realized density surface marked individuals
#  rm = realized_density(samples=out, grid=Grid$grid, crs_=mycrs,
#                    hab_mask=hab_mask, z_alias = "z", s_alias = "s")
#  
#  # generate realized density surface unmarked individuals
#  ru = realized_density(samples=out, grid=Grid$grid, crs_=mycrs,
#                    hab_mask=hab_mask, z_alias = "zu", s_alias = "su")
#  
#  # load needed packages for multiplot
#  library(viridis)
#  library(grid)
#  library(cowplot)
#  library(ggpubr)
#  library(rasterVis)
#  
#  # plot raster for marked individuals
#  p1<-gplot(rm) + geom_raster(aes(fill = value)) +
#            scale_fill_viridis(na.value = NA, name="Density",
#            limits=c(0,0.3),breaks=seq(0,0.3,by=0.1)) +
#            xlab("") + ylab("") + theme_classic() +
#            scale_x_continuous(expand=c(0, 0)) +
#            scale_y_continuous(expand=c(0, 0)) +
#             theme(axis.text = element_text(size=18))
#  
#  # plot raster for unmarked individuals
#  p2<-gplot(ru) + geom_raster(aes(fill = value)) +
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

