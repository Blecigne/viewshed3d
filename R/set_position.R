#' Set animal position elevation
#'
#' @description Updates the elevation of the user defined animal position based on the nearest ground point.
#' The updated elevation correspond to the user defined elevation + the elevation
#' of the nearest ground point.
#'
#' @param data @param data LAS class object containing the xyz coordinates of a 3D point
#' cloud classified as ground and non-ground points.
#' @param position vector of length 3 containing the xyz coordinates of the
#' animal location. Default = c(0,0,0).
#'
#' @return A vector of length 3 containing the updated position.
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
#' new_position <- viewshed3d::set_position(class,center)
#'
#' # input position
#' center
#' # updated position
#' new_position
#' }
set_position = function(data,position){

  #- declare variables to pass CRAN check as suggested by data.table mainaitners
  Classification=NULL

  # select ground points
  ground <- data@data[Classification == 2]

  # find the Z value for the closest ground point
  elevation <- ground$Z[which.min(
    sqrt( (ground$X - position[1])^2 + (ground$Y - position[2])^2 )
  )]

  return(c(position[1:2], # original X and Y
           elevation + position[3])) # modified Z
}
