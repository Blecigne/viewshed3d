
#' Computes the visibility from a single location in a 3D point cloud
#'
#' @description  Computes visibility from a user-defined location with
#' user-defined sightline angles and returns the visibility as function of
#' distance and, optionally, the 3D point cloud classified as visible
#' and non-visible points.
#'
#' @param data LAS class object containing the xyz coordinates of a 3D point
#' cloud
#' @param position vector of length 3 containing the xyz coordinates of the
#' animal location. Default = c(0,0,0).
#' @param angular_res numeric. The angular resolution of a single sightline.
#' Default = 1.
#' @param scene_radius (optional) numeric. Defines the radius of the scene
#' relative to the animal position. Can be used to apply a cut-off distance to
#' visibility analyses.
#' @param store_points logical. If \code{TRUE}, the 3D point cloud is returned
#' with visible and not visible points classified (see details).
#'
#' @note In most cases, a ground reconstruction should be performed before
#' visibility computation to avoid sightlines passing through the ground. This
#' can be done with the \code{\link{reconstruct_ground}} function.
#'
#' @return If \code{store_points = FALSE}, a data.table of the visibility
#' (Visibility) as a function of distance to the animal location (r) is
#' returned. If \code{store_points = TRUE}, a list containing two objects is
#' returned. The first object is similar to the data.table returned when
#' \code{store_points = FALSE}. The second object is a LAS class object
#' containing the coordinates of the point cloud (X, Y, Z), the distance of
#' each point to the animal position (r) and the class of each point: visible
#' or not visible from the animal position (Visibility = 2 or 1, respectively).
#'
#' @details Sightline directions are computed from the method described by
#' Malkin (2016). This ensures that every sightline explores
#' a similar portion of the 3d space.
#'
#' @references Malkin, Z. (2016). A new method to subdivide a spherical surface
#' into equal-area cells. arXiv:1612.03467.
#'
#' @importFrom data.table := .N
#'
#' @export
#'
#' @examples
#' \donttest{
#' #- import the tree_line_plot dataset
#' file <- system.file("extdata", "tree_line_plot.laz", package="viewshed3d")
#' tls <- lidR::readLAS(file)
#'
#' center <- c(0,0,2) # defines the scene center for the entire process
#' angle <- 1 # defines the angular resolution for the entire process
#'
#' #- remove noise to avoid visibility estimates error
#' tls_clean <- viewshed3d::denoise_scene(tls,method="sd",
#'                                        filter=6)
#'
#'
#' #- class ground and vegetation points
#' class <- lidR::classify_ground(tls_clean, lidR::csf(rigidness = 1L,
#'                                               class_threshold = 0.2,
#'                                               sloop_smooth = FALSE))
#'
#' #- reconstruct the ground with the optimal resolution
#' recons <- viewshed3d::reconstruct_ground(data=class,
#'                                          position = center,
#'                                          ground_res = 0.05,
#'                                          angular_res = angle,
#'                                          method="knnidw",
#'                                          full_raster = TRUE)
#'
#' #- compute the visibility and store the output point cloud.
#' #- As the input file is a LAS object, the output
#' #- point cloud is also stored in a LAS file.
#' view.data <- viewshed3d::visibility(data = recons,
#'                                     position = center,
#'                                     angular_res = angle,
#'                                     scene_radius = 17, # apply cut_oof distance
#'                                     store_points = TRUE)
#'
#' #- 3D plot with visible points in white and non-visible points in darkgrey
#' x=lidR::plot(view.data$points,color="Visibility",colorPalette = c("grey24","white"))
#'
#' #- add animal position to the plot
#' position=data.frame(X=center[1],Y=center[2],Z=center[3])
#' lidR::add_treetops3d(x,sp::SpatialPointsDataFrame(position,position),
#'                      radius=0.2,col="red")
#'
#' #- plot the visibility as function of distance
#' plot(view.data$visibility$r,view.data$visibility$visibility,
#'      type="l",ylim=c(0,100),lwd=4)
#'}
#'
visibility <- function(data,position,angular_res,scene_radius,store_points){

  #- declare variables to pass CRAN check as suggested by data.table mainaitners
  N=theta=X=Y=Z=phi=r=.=Visibility=NULL

  #- default parameters
  if(missing(store_points)) store_points <- F
  if(missing(position)) position <- c(0,0,0)
  if(missing(angular_res)) angular_res <- 1

  #- check for possible problems in arguments
  if(missing(data)) stop("no input data provided.")
  if(class(data)[1]!="LAS"){
    stop("data must be a LAS. You can use lidr::LAS(data) to convert a
         data.frame or data.table to LAS.")
  }
  if(is.numeric(angular_res)==F) stop("angular_res must be numeric.")
  if(is.numeric(position)==F) stop("position must be numeric.")
  if(length(position)!=3) stop("position must be a vector of length 3.")
  if(is.logical(store_points)==F) stop("store_coord must be logical")

  #- BUILD PARAMETER FILE
  #- with theta = the elevation angle, and N = the number of cells
  if(180%%angular_res==0){ # build the entire theta range
    param <- data.table::data.table(theta=seq(-90,90,angular_res),N=NA)
  }else{
    param <- data.table::data.table(theta=seq(-90,90+angular_res,angular_res),
                                    N=NA)
  }
  # compute the number of cells in each latitudinal band following the
  # method described by Malkin (2016)
  L=pi/(360/angular_res) # cell side length at the equator
  N_cell=round(pi/L) # number of cells at the equator
  # the number of decreases with increasing the elevation angle
  param[,N:=round(N_cell*(cos(pracma::deg2rad(theta))))]

  # N = 1 for the two extremes and transform theta range from -90, 90 to 0, 180
  param[theta==90 | theta==-90,N:=1]
  param[,theta:=theta+90]

  # transform cumulated number of cells into the number of cells in one row
  n_dir <- sum(param$N)

  #- FIND THE NEAREST POINT IN EACH DIRECTION
  data <- data@data[,.(X,Y,Z)] # transform data into a data.table

  #- translate the data so the center is 0,0,0
  data[, ':=' (X = X - position[1],
               Y = Y - position[2],
               Z = Z - position[3])]

  #- compute spherical coordinates: elevation angle (theta), azimuth (phi) and
  #- radius (r)
  data[, ':=' (theta = round((acos(Z * 100/(sqrt(X^2 + Y^2 + Z^2) * 100)) *
                                (180/pi))/angular_res)*angular_res, #-
               phi = acos(Y * 100/(sqrt(Y^2 + X^2) * 100)) * (180/pi),
               r = sqrt(X^2+Y^2+Z^2))]
  data[ X < 0 , phi := (180-phi)+180] # transform phi so its range is now 0-360

  if(missing(scene_radius)==F){ #- apply the scene radius if defined
    #- ensure scene.radius is valid
    if(is.numeric(scene_radius)==F) stop("scene.radius must be numeric.")
    data=data[r <= scene_radius]
    if(nrow(data)==0) stop("scene_radius is too small, the scene contains 0 point.")
  }

  #- match param and data
  data.table::setkey(param,theta)
  data.table::setkey(data,theta)
  data[, N:= param[data,N]]

  #- compute the phi for each cell center
  data[, phi := round(phi/(360/N))*(360/N)]
  data[ phi>=360 , phi:=0]

  #- find the nearest point in each direction and store xyz coordinates
  near <- data[, .(r = min(r)), keyby = .(theta,phi)]
  near <- unique(near)

  #- store coordinates for export (needed for the cumview function)
  if(store_points==T){
    data[,Visibility := 1]
    data.table::setkeyv(data,c("theta","phi","r"))
    data.table::setkeyv(near,c("theta","phi","r"))
    data[near, Visibility := 2]
    data[,c("theta","phi","N"):=NULL]
    data[, ':=' (X = X + position[1],Y = Y + position[2], Z = Z + position[3])]
  }

  #- COMPUTE VISIBILITY
  #- number of nearest points at each distance
  near[, r := round(r,2)]
  near <- near[, .(N = .N), keyby = r]

  #- add distances from 0 to the minimum recorded
  start <- data.table::data.table(r = seq(0,min(near$r),0.1), N = 0)
  near <- rbind(start, near)

  #- compute visibility
  near[, visibility := (1-cumsum(N)/n_dir)*100]
  near[,N:=NULL]

  if(store_points==T){
    data <- lidR::LAS(data) # export a LAS
    ret <- list(visibility=near,points=data)
    return(ret)
  }else{
    return(near)
  }
}
