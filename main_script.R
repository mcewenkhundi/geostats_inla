#Task: Understand how to safely predict estimates
#Script adapted from https://www.paulamoraga.com/book-geospatial/sec-geostatisticaldataexamplespatial.html
#Date 27May2024

#Using Malaria data in the Gambia

library(geoR) #to get the malaria data
data(gambia)

head(gambia)
dim(gambia)

#unique village location
dim(unique(gambia[, c("x", "y")]))

#data with malaria prevelance by village
library(dplyr)
d <- group_by(gambia, x, y) %>%
  summarize(
    total = n(),
    positive = sum(pos),
    prev = positive / total
  )
head(d)

#Plot the malaria prevalence
library(sp)
library(rgdal)
sps <- SpatialPoints(d[, c("x", "y")],
                     proj4string = CRS("+proj=utm +zone=28")
)

#Transform the UTM coordinates in sps to geographic coordinates using spTransform() where we set CRS to
spst <- spTransform(sps, CRS("+proj=longlat +datum=WGS84"))

#Finally, we add the longitude and latitude variables to the data frame d.
d[, c("long", "lat")] <- coordinates(spst)
head(d)

#Mapping prevalence
library(leaflet)
library(viridis)

pal <- colorBin("viridis", bins = c(0, 0.25, 0.5, 0.75, 1))
leaflet(d) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(lng = ~long, lat = ~lat, color = ~ pal(prev)) %>%
  addLegend("bottomright",
            pal = pal, values = ~prev,
            title = "Prev."
  ) %>%
  addScaleBar(position = c("bottomleft"))

#Get environmental covariates, altitude
library(raster)
#r <- getData(name = "alt", country = "GMB", mask = TRUE)
#Add the altitude to the data.frame d
library(geodata)
elevation_30s(country="GMB", path=here::here())
r <- raster(here::here("GMB_elv_msk.tif"))

#Plot the elevation
pal <- colorNumeric("viridis", values(r),
                    na.color = "transparent"
)

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addRasterImage(r, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright",
            pal = pal, values = values(r),
            title = "Altitude"
  ) %>%
  addScaleBar(position = c("bottomleft"))

#get data into the data.frame
d$alt <- raster::extract(r, d[, c("long", "lat")])

head(d)

#Modelling
library(INLA)
coo <- cbind(d$long, d$lat)
mesh <- inla.mesh.2d(
  loc = coo, max.edge = c(0.1, 5),
  cutoff = 0.01
)

#The number of mesh vertices given by mesh$n
mesh$n
plot(mesh)
points(coo, col = "red")

#Build the spde model on the mesh
spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)
indexs <- inla.spde.make.index("s", spde$n.spde)
lengths(indexs)

#Projection matrix
A <- inla.spde.make.A(mesh = mesh, loc = coo)

#Prediction data
ra <- aggregate(r, fact = 5, fun = mean)
dp <- rasterToPoints(ra)
dim(dp)

#prediction projection
coop <- dp[, c("x", "y")]
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)

#Stack with data fro estimation and prediction
# stack for estimation stk.e
stk.e <- inla.stack(
  tag = "est",
  data = list(y = d$positive, numtrials = d$total),
  A = list(1, A),
  effects = list(data.frame(b0 = 1, altitude = d$alt), s = indexs)
)

# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA, numtrials = NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = 1, altitude = dp[, 3]),
                 s = indexs
  )
)

# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)

#Model formula
formula <- y ~ 0 + b0 + altitude + f(s, model = spde)

#Inla() call
# In control.predictor we set compute = TRUE to compute the posteriors of the 
#predictions. We set link=1 to compute the fitted values (res$summary.fitted.values
#and res$marginals.fitted.values) with the same link function as the family 
#specified in the model. We also add control.compute = list(return.marginals.predictor = TRUE)
#to obtain the marginals.
res <- inla(formula,
            family = "binomial", Ntrials = numtrials,
            control.family = list(link = "logit"),
            data = inla.stack.data(stk.full),
            control.predictor = list(compute = TRUE, link = 1,
                                     A = inla.stack.A(stk.full)),
            control.compute = list(return.marginals.predictor = TRUE)
)

#Mapping prevalence
#The mean prevalence and lower and upper limits of 95% credible intervals are in
#the data frame res$summary.fitted.values. The rows of res$summary.fitted.values 
#that correspond to the prediction locations can be obtained by selecting the indices
#of the stack stk.full that are tagged with tag = "pred".
#We can obtain these indices by using inla.stack.index() passing stk.full and tag = "pred".
index <- inla.stack.index(stack = stk.full, tag = "pred")$data

prev_mean <- res$summary.fitted.values[index, "mean"]
prev_ll <- res$summary.fitted.values[index, "0.025quant"]
prev_ul <- res$summary.fitted.values[index, "0.975quant"]

#Some of the fitted values are negative
#How can we make sure that all the fitted.values are non-negative
#Since they are supposed to be on the probability scale
min(res$summary.fitted.values$mean)
View(res$summary.fitted.values)

#Graph of the prevalence 
pal <- colorNumeric("viridis", c(0, 1), na.color = "transparent")

leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addCircles(
    lng = coop[, 1], lat = coop[, 2],
    color = pal(prev_mean)
  ) %>%
  addLegend("bottomright",
            pal = pal, values = prev_mean,
            title = "Prev."
  ) %>%
  addScaleBar(position = c("bottomleft"))

