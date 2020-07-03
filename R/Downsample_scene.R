#' Reduces the point cloud density
#'
#' @param data LAS file of a 3D point cloud.
#' @param method character string. Defines the method to use for downsampling.
#' Can be "space" or "random". See details. Default = "space".
#' @param filter numeric. The intensity of the filter that depends on the
#' method. See details.
#'
#' @details
#' \subsection{\code{method = "space"}}{
#' a single point is conserved within a voxel of
#' \code{filter} size.}
#' \subsection{\code{method = "random"}}{
#' randomly select a user defined proportion of the point cloud.
#' Here, \code{filter} is the proportion of points to keep in the point cloud.}
#'
#' @return The downsampled data.
#'
#' @importFrom data.table :=
#'
#' @export
#'
#' @examples
#' \donttest{
#' #- import the tree_line_plot dataset
#' file <- system.file("extdata", "tree_line_plot.laz", package="viewshed3d")
#' tls <- lidR::readLAS(file,select="xyz")
#'
#' #- reduce the point cloud density: keep one point in a voxel of 4cm.
#' sub = viewshed3d::downsample_scene(tls,filter=0.04)
#'
#' #- plot the downsampled point cloud
#' lidR::plot(sub)
#'}


downsample_scene <- function(data,method,filter){

  #- declare variables to pass CRAN check as suggested by data.table mainaitners
  X=Y=Z=X_vox=Y_vox=Z_vox=temp=.=NULL

  #- default parameters
  if(missing(method)) method <- "space"
  if(method=="space" & missing(filter)) filter <- 0.02
  if(method=="random" & missing(filter)) filter <- 0.5

  #- checking for potential errors in arguments
  if(missing(data)){stop("no input data provided.")}
  if(class(data)[1]!="LAS"){
    stop("data must be a LAS. You can use lidr::LAS(data) to convert a
         data.frame or data.table to LAS.")
  }
  if(method %in% c("space","random")==F){
    stop("unknow method, method can be 'space' or 'random'.")
  }
  if(is.numeric(filter)==F) stop("filter must be numeric")
  if(method=="random") if(filter<0 | filter>1) stop("filter must be between 0 and 1")
  if(length(filter)!=1) stop("filter must be of length 1.")

  #- data as a data.table
  data <- data@data

  if(method == "random"){
    #- randomly sample the rows to keep
    data <- data[sample(1:nrow(data),round(nrow(data)*filter)),]
  }

  if(method == "space"){
    #- create a voxel cloud with an index
    data[,':='(X_vox = round(X/filter)*filter,
               Y_vox = round(Y/filter)*filter,
               Z_vox = round(Z/filter)*filter,
               row = 1:nrow(data))]

    #- for each voxel keep the point with the greater index
    temp <- data[,.(to_keep = max(row)), by = 'X_vox,Y_vox,Z_vox']

    #- keep the selected points from the original dataset
    data <- data[temp$to_keep,]

    #- remove useless dataframe and columns
    rm(temp)
    data[,':='(X_vox = NULL,Y_vox = NULL,Z_vox = NULL,row = NULL)]
  }

  data <- lidR::LAS(data) # export a LAS
  return(data)
}
