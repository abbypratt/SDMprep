#' @title Occurrence Preparation
#' @description This function prepares occurrence data for species distribution modeling
#' @details This function deletes duplicates and NAs from occurrence data. It also rarefies data by a user-specified distance.
#' @param param1 x = data frame containing occurrence data with column names 'name', 'latitude', 'longitude'.
#' @param param1 y = the distance (in km) at which you would like to rarefy points
#' @returns The output is a dataframe.
#' @export
occ.prep = function(x,y) {
  occs.dups <- duplicated(x[c('longitude', 'latitude')])
  x <- x[!occs.dups,]
  # remove NAs
  occs.full <- x[complete.cases(x$longitude, x$latitude), ]
  # give all records a unique ID
  occs.full$occID <- row.names(occs.full)

  #Spatially Rarefy
  output <- spThin::thin(occs.full, 'latitude', 'longitude', 'name', thin.par = y, reps = 100,
                         locs.thinned.list.return = TRUE, write.files = FALSE, verbose = FALSE)
  # find the iteration that returns the max number of occurrences
  maxThin <- which(sapply(output, nrow) == max(sapply(output, nrow)))
  # if there's more than one max, pick the first one
  maxThin <- output[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]
  # subset occs to match only thinned occs
  occs <- occs.full[as.numeric(rownames(maxThin)),]
}


#' @title Background Bounding Box Formation
#' @description This function creates the spatial extent from which background points will be selected.
#' @details This function finds the minimum and maximum lat and long values from the occurrence points and uses that to create the extent from which background points will be selected.
#' @param param1 x = output from the occ.prep function, or a data frame with occurrence points and column names 'name', 'latitude', 'longitude'.
#' @returns The output is a Formal Class SpatialPolygon.
#' @export
bg.bb = function(x) {
  ##Background selection technique chosen as Bounding Box.
  xmin <- min(x$longitude)
  xmax <- max(x$longitude)
  ymin <- min(x$latitude)
  ymax <- max(x$latitude)
  bb <- matrix(c(xmin, xmin, xmax, xmax, xmin, ymin, ymax, ymax, ymin, ymin), ncol=2)
  bgExt <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(bb)), 1)))
  #Buffer size of the study extent polygon defined as 1 degrees.
  bgExt <- rgeos::gBuffer(bgExt, width = 1)
  assign("bgExt", bgExt, envir = .GlobalEnv)
}

#' @title Background Point Selection
#' @description This function selects background points from the area specified in the bg.bb function.
#' @details This function randomly selects 10000 background points from the area specified in the bg.bb function.
#' @param param1 x = raster stack of environmental variable.
#' @param param2 y = output from the bg.bb function.
#' @param param3 z = number of background points you would like to select.
#' @returns The output is a dataframe with columns 'x', 'y' corresponding to longitude and latitude.
#' @export
# x = env raster stack; y = output from bg.bb function
bg.sel = function(x,y,z) {
  #Mask environmental variables by Bounding Box, and take a random sample of background values from the study extent.
  # crop the environmental rasters by the background extent shape
  envsBgCrop <- raster::crop(x, y)
  assign("envsBgCrop", envsBgCrop, envir = .GlobalEnv)
  # mask the background extent shape from the cropped raster
  envsBgMsk <- raster::mask(envsBgCrop, y)
  assign("envsBgMsk", envsBgMsk, envir = .GlobalEnv)
  # sample random background points
  bg.xy <- dismo::randomPoints(envsBgMsk, z)
  # convert matrix output to data frame
  bg.xy <- as.data.frame(bg.xy)
}

#' @title Data partitioning
#' @description This function partitions both occurrence and background points into groups.
#' @details This function partitions occurrence data and background points into 4 groups to be used in k-fold jackknifing.
#' @param param1 x = output from the occ.prep function, or a data frame with occurrence points and column names 'name', 'latitude', 'longitude'.
#' @param param2 y = output from the bg.sel function.
#' @returns The output is two numeric vectors with numbers 1-4, one vector each for occurrences and background.
#' @export
partition = function(x,y) {
  occs.xy <- x[c('longitude', 'latitude')]
  assign("occs.xy", occs.xy, envir = .GlobalEnv)
  group.data <- ENMeval::get.block(occ = occs.xy, bg = y)
  # pull out the occurrence and background partition group numbers from the list
  occs.grp <- group.data[[1]]
  assign("occs.grp", occs.grp, envir = .GlobalEnv)
  bg.grp <- group.data[[2]]
  assign("bg.grp", bg.grp, envir = .GlobalEnv)
}




