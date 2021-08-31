#load packages
library(maptools)
library(MuMIn)  
library(raster) 
library(rgdal) 
library(rgeos) 
library(sp) 
library(spatstat)
library(here)

#load data from working directory
pare10_points <- readOGR('.', "Pare10")
pare20_points <- readOGR('.', "Pare20")
shore_line <- readOGR('.', "shore limit")
DTM <- raster("tingrid_Copy.tif")
path_distance <- raster("euc_pd_clip.tif") #path distance cost path calculations based on tobler time 
ridge_distance <- raster("euc_rdg_clip.tif")
tot_vis <- raster("5m_tot_vis2_Copy.tif")
water_distance <- raster("euc_wtr_clip.tif") #points based off of comparing GIS hydrology stream order against published map from NAE project. 



#convert to spatstat format 
shore <- as.owin(shore_line) #convert to window format
pare10 <- ppp(pare10_points$POINT_X, pare10_points$POINT_Y, window=shore) #convert pare to a ppp object/ error? maybe use X_Cor
pare20 <- ppp(pare20_points$POINT_X, pare20_points$POINT_Y, window=shore)
elev <- as.im(DTM) #convert DTM to a pixel image
pd_costpath <- as.im(path_distance) #convert g visibility to a pixel image
ridge <- as.im(ridge_distance)
tot_vis <- as.im(tot_vis) #convert total visibility to a pixel image
water <- as.im(water_distance)

###########################
### Exploratory Analyses ##
###########################

#Compute nearest neighbor distances for forts 
pare10_nn <- nndist(pare10)  
#mean nn
mean(pare10_nn)   #mean=1003
#median_nn
median(pare10_nn)    #median=1066.6
#20 pare stats Nearest neighbor
pare20_nn <- nndist(pare20)
#mean nn
mean(pare20_nn)      #mean=762.9
#median_nn
median(pare20_nn)    #median=723.4

#Perform L function test againt 39 realizations of CSR with fixed number of points
set.seed(1234) #set random seed to get reproducible result
pare10_L <- envelope(pare10, fun=Lest,  fix.n=T,nsim=39) 
#check it
plot(pare10_L)

set.seed(1234) #set random seed to get reproducible result
pare20_L <- envelope(pare20, fun=Lest,  fix.n=T,nsim=39) 
#check it
plot(pare20_L)

#Examine form of possible covariate effects of the sample of 10 large pare
#using relative distribution (nonparametric regression (rhohat))
elev_rh10 <- rhohat(pare10, elev, confidence = 0) #intensity as a function of elevation 
pd_costpath_rh10 <- rhohat(pare10, pd_costpath, confidence = 0) #intensity as a function of distance from path distance cost path
ridge_rh10 <- rhohat(pare10, ridge, confidence = 0) #intensity as a function of distance from ridgeline 
tot_vis_rh10 <- rhohat(pare10, tot_vis, confidence = 0) #intensity as a function of total visibility for the island 
water_rh10 <- rhohat(pare10, water, confidence = 0) #intensity as a function of distance to water

#Examine form of possible covariate effects of the sample of 10 large pare
#using relative distribution (nonparametric regression (rhohat))
elev_rh20 <- rhohat(pare20, elev, confidence = 0) #intensity as a function of elevation 
pd_costpath_rh20 <- rhohat(pare20, pd_costpath, confidence = 0) #intensity as a function of distance from path distance cost path
ridge_rh20 <- rhohat(pare20, ridge, confidence = 0) #intensity as a function of distance from ridgeline 
tot_vis_rh20 <- rhohat(pare20, tot_vis, confidence = 0) #intensity as a function of total visibility for the island 
water_rh20 <- rhohat(pare20, water, confidence = 0) #intensity as a function of to water



############
## Models ##
############

# Beginning of point process models for 10 forts             
#correction="none" so models do not calculate unobserved points outside the window, 
#since this is an island poinst could not exist in the ocean  
#the default is a "border correction" which uses the fixed threshold value and has 
#an erroded border (not really sure what it means) but either could be valid. 

# Modeling with Gibbs  Hard-core process: forbidden to come close to oneanother within a specified distance   
#Inhomogeneous models include models that involve covariate effects (like elevation)  - 
# Worth considering a hybrid model as well that combines two different processes
# Need to add an edge correction to make sure it does not account for points "outside" of the window. Appropriate when borders are "real" borders like an island          correction="none"
#possibly use a hybrid of Strauss and Area Interaction? 
#Stick with Strauss, hybrid is better for multiple interactions, but this is really just one interaction between forts. 

#Strauss. range of pairwise interaction determined by looking at Kr throguh mean(K_sim$r) = 942  and the mean of nearest neighbor = 1003 
ppm100 <- ppm(pare10, ~ 1, Strauss(1000), correction = 'none')
ppm101 <- ppm(pare10, ~ pd_costpath, Strauss(1000), correction = 'none')
ppm102 <- ppm(pare10, ~ ridge, Strauss(1000), correction = 'none')
ppm103 <- ppm(pare10, ~ tot_vis, Strauss(1000), correction = 'none')
ppm104 <- ppm(pare10, ~ elev, Strauss(1000), correction = 'none')
ppm105 <- ppm(pare10, ~ water, Strauss(1000), correction = 'none')
ppm106 <- ppm(pare10, ~ pd_costpath+ridge, Strauss(1000), correction= 'none') 
ppm107 <- ppm(pare10, ~ pd_costpath+tot_vis, Strauss(1000), correction= 'none')
ppm108 <- ppm(pare10, ~ pd_costpath+elev, Strauss(1000), correction= 'none')
ppm109 <- ppm(pare10, ~ pd_costpath+water, Strauss(1000), correction= 'none')
ppm110 <- ppm(pare10, ~ ridge+tot_vis, Strauss(1000), correction= 'none')
ppm111 <- ppm(pare10, ~ ridge+elev, Strauss(1000), correction= 'none')
ppm112 <- ppm(pare10, ~ ridge+water, Strauss(1000), correction= 'none')
ppm113 <- ppm(pare10, ~ tot_vis+elev, Strauss(1000), correction= 'none')
ppm114 <- ppm(pare10, ~ tot_vis+water, Strauss(1000), correction= 'none')
ppm115 <- ppm(pare10, ~ elev+water, Strauss(1000), correction= 'none')
ppm116 <- ppm(pare10, ~ pd_costpath+ridge+tot_vis, Strauss(1000), correction= 'none')
ppm117 <- ppm(pare10, ~ pd_costpath+ridge+elev, Strauss(1000), correction= 'none')
ppm118 <- ppm(pare10, ~ pd_costpath+ridge+water, Strauss(1000), correction= 'none')
ppm119 <- ppm(pare10, ~ pd_costpath+tot_vis+elev, Strauss(1000), correction= 'none')
ppm120 <- ppm(pare10, ~ pd_costpath+tot_vis+water, Strauss(1000), correction= 'none')
ppm121 <- ppm(pare10, ~ pd_costpath+elev+water, Strauss(1000), correction= 'none')
ppm122 <- ppm(pare10, ~ ridge+tot_vis+elev, Strauss(1000), correction= 'none')
ppm123 <- ppm(pare10, ~ ridge+tot_vis+water, Strauss(1000), correction= 'none')
ppm124 <- ppm(pare10, ~ ridge+elev+water, Strauss(1000), correction= 'none')
ppm125 <- ppm(pare10, ~ tot_vis+elev+water, Strauss(1000), correction= 'none')
ppm126 <- ppm(pare10, ~ pd_costpath+ridge+tot_vis+elev, Strauss(1000), correction= 'none')
ppm127 <- ppm(pare10, ~ pd_costpath+ridge+elev+water, Strauss(1000), correction= 'none')
ppm128 <- ppm(pare10, ~ pd_costpath+ridge+tot_vis+water, Strauss(1000), correction= 'none')
ppm129 <- ppm(pare10, ~ pd_costpath+tot_vis+elev+water, Strauss(1000), correction= 'none')
ppm130 <- ppm(pare10, ~ ridge+tot_vis+elev+water, Strauss(1000), correction= 'none')
ppm131 <- ppm(pare10, ~ pd_costpath+ridge+tot_vis+elev+water, Strauss(1000), correction= 'none')


#compare Gibbs models
ppm_AICc10 <- model.sel(ppm100, ppm101, ppm102, ppm103, ppm104, ppm105, ppm106, ppm107, ppm108, ppm109, ppm110, ppm111, ppm112, ppm113, ppm114, ppm115, ppm116, ppm117, ppm118, ppm119, ppm120, ppm121, ppm122, ppm123, ppm124, ppm125, ppm126, ppm127, ppm128, ppm129, ppm130, ppm131,rank=AICc)
ppm_AICc10

#print results of best fitting model
summary(ppm110)

#check residual K
set.seed(1234)
K_sim10 <- envelope(ppm110, Kres, nsim=99, fix.n=T)
#check fit
plot(K_sim10, lwd=3, legend='F')

#plot predicted first-order intensity of best-fitting model as a geographic map
plot(intensity.ppm(ppm110), 
     col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, rev = T),
     main="", riblab="Fitted intensity")
plot(pare10, pch=16, col='red', cex=0.25, add=T)

#partial residual plot for elevation
par_res_vis10 <- parres(ppm110, "tot_vis")
#check fit
plot(par_res_vis10)

#partial residual plot for elevation
par_res_ridge10 <- parres(ppm110, "ridge")
#check fit
plot(par_res_ridge10)

#PPM for 20 pare

#HybridModel       Strauss may not be the best interaction for this since sattelite communities can be more close. Perhaps Area Interaction or a hybrid model of Strauss and Area Interaction. Hybrids are good at modeling interaction at multiple scales (i.e. the large pare and the small forts)  hybrid <- hybrid(Strauss(1000), AreaInter(700))

## ALTERNATIVE HYBRID MODEL
#Strauss. range of pairwise interaction determined by looking at a plot of pare20_nn  = 700       Strauss may not be the best interaction for this since sattelite communities can be more close. Perhaps Area Interaction or a hybrid model of Strauss and Area Interaction. Hybrids are good at modeling interaction at multiple scales (i.e. the large pare and the small forts)  hybrid <- hybrid(Strauss(1000), AreaInter(700))

hybrid20 <- Hybrid(Strauss(1000), AreaInter(740))

ppm200 <- ppm(pare20, ~ 1, hybrid20, correction = 'none')
ppm201 <- ppm(pare20, ~ pd_costpath, hybrid20, correction = 'none')
ppm202 <- ppm(pare20, ~ ridge, hybrid20, correction = 'none')
ppm203 <- ppm(pare20, ~ tot_vis, hybrid20, correction = 'none')
ppm204 <- ppm(pare20, ~ elev, hybrid20, correction = 'none')
ppm205 <- ppm(pare20, ~ water, hybrid20, correction = 'none')
ppm206 <- ppm(pare20, ~ pd_costpath+ridge, hybrid20, correction= 'none') 
ppm207 <- ppm(pare20, ~ pd_costpath+tot_vis, hybrid20, correction= 'none')
ppm208 <- ppm(pare20, ~ pd_costpath+elev, hybrid20, correction= 'none')
ppm209 <- ppm(pare20, ~ pd_costpath+water, hybrid20, correction= 'none')
ppm210 <- ppm(pare20, ~ ridge+tot_vis, hybrid20, correction= 'none')
ppm211 <- ppm(pare20, ~ ridge+elev, hybrid20, correction= 'none')
ppm212 <- ppm(pare20, ~ ridge+water, hybrid20, correction= 'none')
ppm213 <- ppm(pare20, ~ tot_vis+elev, hybrid20, correction= 'none')
ppm214 <- ppm(pare20, ~ tot_vis+water, hybrid20, correction= 'none')
ppm215 <- ppm(pare20, ~ elev+water, hybrid20, correction= 'none')
ppm216 <- ppm(pare20, ~ pd_costpath+ridge+tot_vis, hybrid20, correction= 'none')
ppm217 <- ppm(pare20, ~ pd_costpath+ridge+elev, hybrid20, correction= 'none')
ppm218 <- ppm(pare20, ~ pd_costpath+ridge+water, hybrid20, correction= 'none')
ppm219 <- ppm(pare20, ~ pd_costpath+tot_vis+elev, hybrid20, correction= 'none')
ppm220 <- ppm(pare20, ~ pd_costpath+tot_vis+water, hybrid20, correction= 'none')
ppm221 <- ppm(pare20, ~ pd_costpath+elev+water, hybrid20, correction= 'none')
ppm222 <- ppm(pare20, ~ ridge+tot_vis+elev, hybrid20, correction= 'none')
ppm223 <- ppm(pare20, ~ ridge+tot_vis+water, hybrid20, correction= 'none')
ppm224 <- ppm(pare20, ~ ridge+elev+water, hybrid20, correction= 'none')
ppm225 <- ppm(pare20, ~ tot_vis+elev+water, hybrid20, correction= 'none')
ppm226 <- ppm(pare20, ~ pd_costpath+ridge+tot_vis+elev, hybrid20, correction= 'none')
ppm227 <- ppm(pare20, ~ pd_costpath+ridge+elev+water, hybrid20, correction= 'none')
ppm228 <- ppm(pare20, ~ pd_costpath+ridge+tot_vis+water, hybrid20, correction= 'none')
ppm229 <- ppm(pare20, ~ pd_costpath+tot_vis+elev+water, hybrid20, correction= 'none')
ppm230 <- ppm(pare20, ~ ridge+tot_vis+elev+water, hybrid20, correction= 'none')
ppm231 <- ppm(pare20, ~ pd_costpath+ridge+tot_vis+elev+water, hybrid20, correction= 'none')


#compare Gibbs models
ppm_AICc20 <- model.sel(ppm200, ppm201, ppm202, ppm203, ppm204, ppm205, ppm206, ppm207, ppm208, ppm209, ppm210, ppm211, ppm212, ppm213, ppm214, ppm215, ppm216, ppm217, ppm218, ppm219, ppm220, ppm221, ppm222, ppm223, ppm224, ppm225, ppm226, ppm227, ppm228, ppm229, ppm230, ppm231,rank=AICc)
ppm_AICc20


#print results of best fitting model
summary(ppm210)

#check residual K
set.seed(1234)
K_sim20 <- envelope(ppm210, Kres, nsim=99, fix.n=T)
#check fit
plot(K_sim20, lwd=3, legend='F')

#plot predicted first-order intensity of best-fitting model
plot(intensity.ppm(ppm210), 
     col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, rev = T),
     main="", riblab="Fitted intensity")
plot(pare20, pch=16, col='red', cex=0.25, add=T)

#partial residual plot for visibility
par_res_vis20 <- parres(ppm210, "tot_vis")
#check fit
plot(par_res_vis20)

#partial residual plot for ridge
par_res_ridge20 <- parres(ppm110, "ridge")
#check fit
plot(par_res_ridge20)


##########################
##########################
#### Figures Elements ####
##########################
##########################
#Nearest neighbor and L function plots

par(mfrow=c(1,2))
hist(pare10_nn, xlab="meters", breaks=8, col="grey", xlim=c(0,1500), ylim=c(0,5), main="") #creates histogram
mtext(side=3, line=1, at=-1, adj=0, cex=1, "a)")  # creates the textline above the figure
hist(pare20_nn, xlab="meters", breaks=10, col="grey", xlim=c(0,1500), ylim=c(0,5), main="") #creates histogram
mtext(side=3, line=1, at=-1, adj=0, cex=1, "b)")  # creates the textline above the figure
par(mfrow=c(1,1))


par(mfrow=c(1,2))
plot(pare10_L, lwd=3, xlab="r (meters)", main="", legend=F)
mtext(side=3, line=1, at=-1, adj=0, cex=1, "a)")
plot(pare20_L, lwd=3, xlab="r (meters)", main="", legend=F)
mtext(side=3, line=1, at=-1, adj=0, cex=1, "b)")
par(mfrow=c(1,1))


#plot relative distributions of pare10        SOmething needs fixed here with the first line in order to create the file

par(mfrow=c(2,3))
plot(elev_rh10, legend=F, xlab=("Elevation ASL"), main="", xlim=c(0,500), lwd=3)
mtext(side=3, line=1, at=0, adj=0, cex=0.9, "a)")
plot(water_rh10, legend=F, xlab=("Water"), main="", lwd=3)
mtext(side=3, line=1, at=0, adj=0, cex=0.9, "b)")
plot(tot_vis_rh10, legend=F, xlab=("Visibility"), main="", lwd=3)
mtext(side=3, line=1, at=2550, adj=0, cex=0.9, "c)")
plot(pd_costpath_rh10, legend=F, xlab=("Cost Path"), main="", lwd=3)
mtext(side=3, line=1, at=0, adj=0, cex=0.9, "d)")
plot(ridge_rh10, legend=F, xlab=("Ridges"), main="", lwd=3)
mtext(side=3, line=1, at=0, adj=0, cex=0.9, "e)")
par(mfrow=c(1,1))


#plot relative distributions of pare20 

par(mfrow=c(2,3))
plot(elev_rh20, legend=F, xlab=("Elevation ASL"), main="", xlim=c(0,500), lwd=3)
mtext(side=3, line=1, at=0, adj=0, cex=0.9, "a)")
plot(water_rh20, legend=F, xlab=("Water"), main="", lwd=3)
mtext(side=3, line=1, at=0, adj=0, cex=0.9, "b)")
plot(tot_vis_rh20, legend=F, xlab=("Visibility"), main="", lwd=3)
mtext(side=3, line=1, at=2550, adj=0, cex=0.9, "c)")
plot(pd_costpath_rh20, legend=F, xlab=("Cost Path"), main="", lwd=3)
mtext(side=3, line=1, at=0, adj=0, cex=0.9, "d)")
plot(ridge_rh20, legend=F, xlab=("Cost Path"), main="", lwd=3)
mtext(side=3, line=1, at=0, adj=0, cex=0.9, "e)")
par(mfrow=c(1,1))


#residual diagnostic plots  c & d don't really fit here. Maybe witch them for something else. 

par(mfrow=c(2,2))
plot(K_sim10, lwd=3, main='', legend=F, xlab="r (meters)") #model 110
mtext(side=3, line=1, at=0, adj=0, cex=1, "a)")
plot(K_sim20, lwd=3, main='', legend=F, xlab="r (meters)") #model 210
mtext(side=3, line=1, at=0, adj=0, cex=1, "b)")
par(mfrow=c(1,1))
plot(par_res_ridge10, legend=F, main='', xlab="Visibility")  
mtext(side=3, line=1, at=0, adj=0, cex=1, "c)")
par(mar=c(2.25,2.25,2.25,2.25))
plot(intensity.ppm(ppm110), 
     col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, rev = T),
     main="") #riblab="Fitted intensity")
plot(pare10, pch=16, col='red', cex=0.5, add=T)
mtext(side=3, line=1, at=759002, adj=0, cex=1, "d)")
par(mfrow=c(1,1))


# Partial Residual plots for ridges and visibility 
par(mfrow=c(1,2))
plot(par_res_ridge10, lwd=3, xlab="r (meters)", main="", legend=F)
mtext(side=3, line=1, at=-1, adj=0, cex=1, "a)")
plot(par_res_vis10, lwd=3, xlab="r (meters)", main="", legend=F)
mtext(side=3, line=1, at=-1, adj=0, cex=1, "b)")
par(mfrow=c(1,1))


par(mfrow=c(1,2))
plot(par_res_ridge20, lwd=3, xlab="r (meters)", main="", legend=F)
mtext(side=3, line=1, at=-1, adj=0, cex=1, "a)")
plot(par_res_vis20, lwd=3, xlab="r (meters)", main="", legend=F)
mtext(side=3, line=1, at=-1, adj=0, cex=1, "b)")
par(mfrow=c(1,1))


###################
## Basic Maps ##### need work to use
###################

#plots maps with pare and the covariate
par(mfrow=c(2,2))
plot(elev, col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, rev = FALSE), riblab="Elevation")
plot(pare20, pch=16, col='green', add=T)
plot(pare10, pch=16, col='red', add=T)
plot(shore, add=T)
plot(water, col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, rev = FALSE), riblab="Water")
plot(pare20, pch=16, col='green', add=T)
plot(pare10, pch=16, col='red', add=T)
plot(shore, add=T)
plot(tot_vis, col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, rev = FALSE), riblab="Visibility")
plot(pare20, pch=16, col='green', add=T)
plot(pare10, pch=16, col='red', add=T)
plot(shore, add=T)
plot(pd_costpath, col=gray.colors(10, start = 0.3, end = 0.9, gamma = 2.2, rev = FALSE), riblab="Cost Path")
plot(pare20, pch=16, col='green', add=T)
plot(pare10, pch=16, col='red', add=T)
plot(shore, add=T)
par(mfrow=c(1,1))