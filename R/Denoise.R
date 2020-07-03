
#' Filters isolated points from a point cloud
#'
#' @param data LAS file of a 3D point cloud.
#' @param method character string. Defines the method to use for noise filtering.
#' Can be "quantile", "sd" or "voxel". See details. Default = "sd".
#' @param filter numeric. The intensity of the filter that depends on the
#' method. See details.
#' @param k numeric. The number of nearest neighbours to use to compute the mean
#' nearest neighbour distance. Required only if
#' method = "quantile" or "sd". Default = 5.
#' @param store_noise logical. If TRUE, the surveyed points considered as noise are not
#' removed from the data and a column "Noise" is added with a value of 1
#' indicating non-noisy points and a value of 2 indicating noisy points.
#' Default = FALSE.
#'
#' @details
#' \subsection{\code{method = "quantile"}}{
#' the quantile-based method computes the distance of the k
#' nearest neighbours for each surveyed point and considers points that fall in
#' the last user defined quantile as noise. If quantile is used as the filtering
#' method, the default is set to = 0.999.}
#' \subsection{\code{method = "sd"}}{
#' the standard deviation-based method computes the average
#' distance of the k nearest neighbours of each surveyed point and considers
#' points as noise if they are more than the average distance plus a number of
#' times the standard deviation away from other surveyed points. The filter
#' parameter sets the standard deviation multiplier. Default = 4. This filter is
#' similar to the "SOR filter" available in
#' \href{https://www.cloudcompare.org/doc/wiki/index.php?title=SOR_filter}{CloudCompare}.}
#' \subsection{\code{method = "voxel"}}{
#' the voxel-based method considers surveyed points as noise if
#' they are the only surveyed point within a user defined voxel volume. The
#' \code{filter} parameter sets the voxel size (i.e., voxel side length).
#' Default = 0.5.}
#'
#' @return The filtered data (if \code{store_noise = FALSE}) or the classified
#' data (if \code{store_noise = TRUE}) with noisy points labeled as 2.
#'
#' @importFrom data.table := .N
#'
#' @export
#'
#' @examples
#' \donttest{
#' #- import the tree_line_plot dataset
#' file <- system.file("extdata", "tree_line_plot.laz", package="viewshed3d")
#' tls <- lidR::readLAS(file,select="xyz")
#'
#' #- remove duplicated points
#' tls <- lidR::filter_duplicates(tls)
#'
#' #- filter noise with the quantile base method
#' data <- viewshed3d::denoise_scene(tls,
#'                                   method="quantile",
#'                                   filter=0.999,
#'                                   k=5,
#'                                   store_noise = TRUE)
#'
#' lidR::plot(data,color="Noise",colorPalette=c("white","red")) # plot
#'
#' #- filter noise with the standard deviation based method
#' data <- viewshed3d::denoise_scene(tls,
#'                                   method="sd",
#'                                   filter=4,
#'                                   k=5,
#'                                   store_noise = TRUE)
#'
#' lidR::plot(data,color="Noise",colorPalette=c("white","red")) # plot
#'
#' #- filter noise with the voxel based method
#' data <- viewshed3d::denoise_scene(tls,
#'                                   method="voxel",
#'                                   filter=0.5,
#'                                   store_noise = TRUE)
#' lidR::plot(data,color="Noise",colorPalette=c("white","red")) # plot
#'}
denoise_scene = function(data,method,filter,k,store_noise){

  #- declare variables to pass CRAN check as suggested by data.table mainaitners
  Noise=dist=X=Y=Z=X_vox=Y_vox=Z_vox=N=.=.N=NULL

  #- default parameters
  if(missing(store_noise)) store_noise=FALSE
  if(missing(method)) method <- "sd"
  if(method=="quantile"){
    if(missing(k)) k <- 5
    if(missing(filter)) filter <- 0.999
  }
  if(method=="sd"){
    if(missing(k)) k <- 5
    if(missing(filter)) filter <- 4
  }
  if(method=="voxel" & missing(filter)) filter <- 0.5

  #- checking for potential errors in arguments
  if(missing(data)){stop("no input data provided.")}
  if(is.logical(store_noise)==F) stop("store noise must be logical.")
  if(class(data)[1]!="LAS"){
    stop("data must be a LAS. You can use lidr::LAS(data) to convert a
         data.frame or data.table to LAS.")
  }
  if(method %in% c("sd","quantile","voxel")==F){
    stop("unknow method, method can be 'sd', 'quantile' or 'voxel'.")
  }
  if(method=="quantile" | method == "sd"){
    if(length(k)==F) stop("k must be of length 1.")
    if(is.numeric(k)==F) stop("k must be numeric.")
  }
  if(is.numeric(filter)==F) stop("filter must be numeric")
  if(method=="quantile") if(filter<0 | filter>1) stop("filter must be between 0 and 1")
  if(length(filter)==F) stop("filter must be of length 1.")

  data <- data@data
  data[,Noise := 1] #- create a field for noise classification

  if(method == "quantile" | method == "sd"){
    #- computes k nearest neighbours distance
    #- k = k+1 because the first column is the point itself.
    data[,dist:=rowSums(nabor::knn(data=data, k=(k+1))$nn.dists)/k]

    # points above the threshold are noise
    if(method == "quantile") data[dist>stats::quantile(dist,filter),Noise := 2]

    if(method == "sd"){
      mean_d <- mean(data[,dist]) #- mean distance among all points
      sd <- sd(data[,dist]) #- standard deviation of distance
      #- points above the threshold are noise
      data[dist > (mean_d+(filter*sd)),Noise := 2]
    }
    data[,dist := NULL]
  }

  if(method == "voxel"){
    #- transform data coordinates into voxels
    data[,":=" (X_vox=round(X/filter)*filter,
                     Y_vox=round(Y/filter)*filter,
                     Z_vox=round(Z/filter)*filter)]
    #- count the number of point per voxel
    vox_dat <- data[, .N, keyby = .(X_vox, Y_vox, Z_vox)]

    #- find isolated points
    vox_dat <- vox_dat[N==1]

    #- add Noise = 2 to isolated points
    data.table::setkeyv(data,c("X_vox", "Y_vox", "Z_vox"))
    data.table::setkeyv(vox_dat,c("X_vox", "Y_vox", "Z_vox"))
    data[vox_dat, Noise := 2]

    #- remove the voxelised coordinates from data
    data[,":=" (X_vox=NULL,
                Y_vox=NULL,
                Z_vox=NULL)]
  }

  if(store_noise==FALSE){ #- remove noise from the data
    data <- data[Noise==1]
    data <- data[,Noise:=NULL]
  }

  data <- lidR::LAS(data) # export a LAS
  return(data)
}

