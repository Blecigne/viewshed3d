#' Put slope back into data
#'
#' @description Reverse the result of the \code{\link[viewshed3d]{remove_slope}} function.
#'
#' @param data a LAS of a TLS scene from which the slope was previously removed.
#'
#' @return a LAS with slope.
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
#' # add the slope
#' slope_back <- viewshed3d::add_slope(no_slope$data)
#'}
add_slope = function(data){

  #- declare variables to pass CRAN check as suggested by data.table mainaitners
  X_orig=Y_orig=Z_orig=NULL

  if(missing(data)) stop("no input data provided.")
  if(class(data)[1]!="LAS"){
    stop("data must be a LAS. You can use lidr::LAS(data) to convert a
         data.frame or data.table to LAS.")
  }

  if(!all(c("X_orig","Y_orig","Z_orig") %in% colnames(data@data))){
    stop("The data does not contain original coordinates.")
  }

  #- add new positions in the original data and reset Z to the original
  data@data[,':='(X=X_orig,Y=Y_orig,Z=Z_orig)]
  data@data[,':='(X_orig=NULL,Y_orig=NULL,Z_orig=NULL)]

  return(data)
}
