#' Removes the slope from a TLS scene
#'
#' @description Removes the slope from a TLS scene by fitting a plan to the ground
#' points and rotating the scene according to the plan orientation. Note that the
#' animal position is also rotated so the relative position remain unchanged.
#'
#' @param data LAS class object containing the xyz coordinates of a 3D point
#' cloud with segmented ground and non-ground points.
#' @param position vector of length 3 containing the xyz coordinates of the
#' animal location. Default = c(0,0,0).
#'
#' @return The LAS of the TLS scene without slope and the new animal position. Note
#' that in the LAS, the original coordinates are stored.
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
#'                                                  class_threshold = 0.2,
#'                                                  sloop_smooth = FALSE))
#'
#'
#' # remove the slope from the classified scene
#' no_slope <- viewshed3d::remove_slope(data = class,c(0,0,0))
#'
#' # plot the data with slope
#' lidR::plot(class,axis = TRUE)
#'
#' # plot the data vithout slope
#' lidR::plot(no_slope$data,axis = TRUE)
#'
#' # compute visibility from the data without slope
#' view.data <- viewshed3d::visibility(data = no_slope$data,
#'                                     position = no_slope$position,
#'                                     angular_res = angle,
#'                                     scene_radius = 17, # apply cut_oof distance
#'                                     store_points = TRUE)
#'}

remove_slope = function(data,position){

  #- declare variables to pass CRAN check as suggested by data.table mainaitners
  X=Y=Z=Classification=dist=index=.=x_round=y_round=NULL

  #- default parameters
  if(missing(position)) position <- c(0,0,0)

  #- check for possible problems in arguments
  if(missing(data)) stop("no input data provided.")
  if(class(data)[1]!="LAS"){
    stop("data must be a LAS. You can use lidr::LAS(data) to convert a
         data.frame or data.table to LAS.")
  }
  if(is.numeric(position)==F) stop("position must be numeric.")
  if(length(position)!=3) stop("position must be a vector of length 3.")
  if(!"Classification" %in% colnames(data@data)) stop("Ground points must be classified.")

  # store original coordinates in three new columns
  orig_coords <- data@data[,.(X , Y , Z)]
  data.table::setnames(orig_coords,c("X_orig","Y_orig","Z_orig"))

  data@data <- cbind(data@data,orig_coords)

  # select ground points
  ground <- data@data[Classification == 2]

  # find the ground centre elevation
  ground[,dist := sqrt((X-position[1])^2 + (Y-position[2])^2)]
  gr_Z <- mean(ground[dist == min(dist)]$Z)

  # translate the ground so the origin is on the ground
  ground[,Z:=Z-gr_Z]

  # homogenize the ground point distribution to speed up the plane fitting
  ground[,':='(x_round = round(X/0.2)*0.2,y_round = round(Y/0.2)*0.2,index=1:nrow(ground))]
  ground <- ground[ground[,max(index),by=.(x_round,y_round)]$V1]

  # fit a plan to the ground and store the angles of the plane with X and Y axes
  fit_plan <- hyper.fit::hyper.fit(ground[,1:3],coord.type = "theta")$parm

  rm(ground)

  ######- rotate the scene to remove the slopes
  #- moove the scene to the same position than the ground
  data_mat <- as.matrix(data@data[,.(X,Y,Z)])
  data_mat[,3] <- data_mat[,3]-gr_Z

  #- apply rotations
  data_mat <- rgl::rotate3d(data_mat,-pracma::deg2rad(fit_plan[2]),1,0,0)
  data_mat <- rgl::rotate3d(data_mat,-pracma::deg2rad(fit_plan[1]),0,1,0)

  #- add new positions in the original data and reset Z to the original
  data@data[,':='(X=data_mat[,1],Y=data_mat[,2],Z=data_mat[,3]+gr_Z)]

  data@data[,':='(X_orig = orig_coords$X , Y_orig = orig_coords$Y , Z_orig = orig_coords$Z)]
  rm(data_mat)

  #- rotation of the position
  position[3] <- position[3]-gr_Z
  position <- rgl::rotate3d(position,-pracma::deg2rad(fit_plan[2]),1,0,0)
  position <- rgl::rotate3d(position,-pracma::deg2rad(fit_plan[1]),0,1,0)
  position[3] <- position[3]+gr_Z

  #data@data[,':='(Xrot = fit_plan[2] , Yrot = fit_plan[1] , ground_Z = gr_Z)]

  return(list(data=data,position=position))
}
