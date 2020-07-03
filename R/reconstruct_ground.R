
#' Optimal ground reconstruction for visibility computation
#'
#' @description Reconstructs the ground surface with a grid resolution defined
#' by the user and adds a second grid around the animal position with an
#' optimal resolution so that no sightline can pass through the ground when
#' computing visibility with the \code{\link{visibility}} function.
#'
#' @param data LAS class object containing a 3d point cloud + a Classification
#' field that classes points as ground and non-ground, as provided by the
#' \code{\link[lidR]{classify_ground}} function from the
#' \code{\link[lidR]{lidR-package}}.
#' @param ground_res numeric. The grid resolution to reconstruct the ground on
#' the entire 3D scene. Default = 0.05. NOTE: a if needed, second grid may be added
#' with smaller (internally computed) resolution.
#' @param position vector of length 3 containing the xyz coordinates of the
#' animal position when computing the visibility with the
#' \code{\link{visibility}} function. Default = c(0,0,0).
#' @param angular_res numeric. The angular resolution of sightlines when
#' computing the visibility with the \code{\link{visibility}} function.
#' Default = 1.
#' @param method which algorithm to use for spatial interpolation. Can be
#' "knnidw", "tin" or "kriging". See documentation from the
#' \code{\link[lidR]{lidR-package}} for \code{\link[lidR]{knnidw}},
#' \code{\link[lidR]{tin}} and \code{\link[lidR]{kriging}}.
#' @param full_raster should the entire raster be interpolated for the ground
#' portion around the animal position? Parameter passed to the
#' \code{\link[lidR]{grid_terrain}} function available in the
#' \code{\link[lidR]{lidR-package}}.
#' @param ... other arguments to pass to the spatial interpolation algorithm.
#' See documentation from \code{\link[lidR]{knnidw}}, \code{\link[lidR]{tin}}
#' and \code{\link[lidR]{kriging}}
#'
#' @return A LAS class object containing the 3D point cloud coordinates
#' with the ground reconstructed to be passed directly to the
#' \code{\link{visibility}} function. Note: the Classification field is
#' preserved.
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
#' #- class ground and vegetation points
#' class <- lidR::classify_ground(tls, lidR::csf(rigidness = 1L,
#'                                         class_threshold = 0.2,
#'                                         sloop_smooth = FALSE))
#'
#' #- reconstruct the ground. Here the ground is reconstructed with the user
#' #- defined resolution only.
#' recons <- viewshed3d::reconstruct_ground(data=class,position = c(0,0,3),
#'                                          ground_res = 0.05,
#'                                          angular_res = 2,
#'                                          method="knnidw")
#'
#' lidR::plot(recons,color="Classification",
#'            colorPalette = c("darkgreen","chocolate4"))
#'
#' #- when the position is closer to the ground, the user defined resolution is
#' #- not sufficient and a second grid is added with the optimal resolution so
#' #- that no sightline can pass trough the ground when computing visibility.
#' #- In this example, full_raster = TRUE was used as a portion of the ground
#' #- near the animal location is not reconstructed because of a data
#' #- gap around a TLS scan position when using full_raster = FALSE.
#' recons <- viewshed3d::reconstruct_ground(data=class,position = c(0,0,1),
#'                                          ground_res = 0.05,
#'                                          angular_res = 2,
#'                                          method="knnidw",
#'                                          full_raster = TRUE)
#'
#' lidR::plot(recons,color="Classification",
#'            colorPalette = c("darkgreen","chocolate4"))
#'}
reconstruct_ground <- function(data,ground_res,position,angular_res,
                               method,full_raster,...){

  #- declare variables to pass CRAN check as suggested by data.table mainaitners
  X=Y=Z=x=y=Classification=.=r=NULL

  #- defaults parameters
  if(missing(angular_res)) angular_res <- 1
  if(missing(ground_res)) ground_res <- 0.05
  if(missing(position)) position <- c(0,0,0)
  if(missing(method)) method <- "knnidw"
  if(missing(full_raster)) full_raster <- F

  #- possible errors
  if(missing(data)) stop("no input data provided.")
  if(class(data)[1]!="LAS") stop("data must be a LAS. You can use
                                 lidr::LAS(data) to convert a data.frame or
                                 data.table to LAS.")
  if(is.numeric(ground_res)==F) stop("ground_res must be numeric.")
  if(is.numeric(angular_res)==F) stop("angular_res must be numeric.")
  if(is.numeric(position)==F) stop("position must be numeric.")
  if(length(position)!=3) stop("position must be a vector of length 3.")
  if(is.logical(full_raster)==F) stop("full_raster must be logical.")
  if(is.character(method)==F) stop("method should be a character string.")

  data <- data@data[,.(X,Y,Z,Classification)] # transform data into a data.table

  #- separate ground and vegetation points
  vegetation=data[Classification == 1]
  ground=data[Classification == 2]

  #- compute the ground points distance to the animal location
  ground[,r := sqrt((X-position[1])^2+(Y-position[2])^2+(Z-position[3])^2)]

  #- estimate the scene radius
  scene_radius <- sqrt((ground$X-position[1])^2+(ground$Y-position[2])^2)
  scene_radius <- max(scene_radius,na.rm=T)

  #- compute the required ground resolution based on distance to location and
  #- angular_resolution so that no sightline can pass trough the ground when
  #- using the visibility() function
  distances  <-  seq(min(ground[,r]),max(ground[,r]),0.1) #- range of distance
  #- the required ground resolution at each distance
  spacing  <-  (distances*tan(pracma::deg2rad(angular_res)))/1.5

  #- if the required ground resolution is smaller than the user defined ground
  #- resolution -> reconstruct the ground with a finest ground resolution until
  #- the user defined resolution is sufficient
  if(min(spacing)<ground_res){
    #- the portion of the ground where finest reconstruction is needed
    ground.fine.recons  <- ground[r<=distances[min(which(spacing>=ground_res))]]
  }
  #- renconstruct the ground with the user defined resolution and then,
  #- if needed, with the calculated fine resolution using the lidR package tools
  if(method == "knnidw"){
    #- reconstruct the ground for the entire plot
    dtm.all <- lidR::grid_terrain(lidR::LAS(ground[,1:4]),
                                  ground_res,
                                  lidR::knnidw(...),
                                  full_raster = TRUE)
    if(min(spacing)<ground_res){
      #- reconstruct the portion of the ground that require a fine
      #- reconstruction with the parameters calculated previously
      dtm.fine <- lidR::grid_terrain(lidR::LAS(ground.fine.recons[,1:4]),
                                     min(spacing),
                                     lidR::knnidw(...),
                                     full_raster=full_raster)
    }
  }
  if(method == "tin"){
    #- reconstruct the ground for the entire plot
    dtm.all <- lidR::grid_terrain(lidR::LAS(ground[,1:4]),
                                  ground_res,
                                  lidR::tin(...),
                                  full_raster = TRUE)
    if(min(spacing)<ground_res){
      #- reconstruct the portion of the ground that require a fine
      #- reconstruction with the parameters calculated previously
      dtm.fine = lidR::grid_terrain(lidR::LAS(ground.fine.recons[,1:4]),
                                    min(spacing),
                                    lidR::tin(...),
                                    full_raster=full_raster)
    }
  }
  if(method == "kriging"){
    #- reconstruct the ground for the entire plot
    dtm.all <- lidR::grid_terrain(lidR::LAS(ground[,1:4]),
                                  ground_res,
                                  lidR::kriging(...),
                                  full_raster = TRUE)
    if(min(spacing)<ground_res){
      #- reconstruct the portion of the ground that require a fine
      #- reconstruction with the parameters calculated previously
      dtm.fine <- lidR::grid_terrain(lidR::LAS(ground.fine.recons[,1:4]),
                                     min(spacing),
                                     lidR::kriging(...),
                                     full_raster=full_raster)
    }
  }

  #- transform the reconstructed ground into a data.table and remove missing Z
  recons.all <- data.table::data.table(raster::as.data.frame(dtm.all,xy = TRUE))
  recons.all <- recons.all[is.na(Z)==F]

  #- keep the points within the original scene ridius
  recons.all[,r := sqrt((x-position[1])^2+(y-position[2])^2)]
  recons.all <- recons.all[r <= scene_radius]
  recons.all[,r:=NULL]

  if(min(spacing)<ground_res){
    #- transform the ground reconstructed with fine resolution into a data.table
    #- and remove missing Z
    recons.fine <- data.table::data.table(raster::as.data.frame(dtm.fine,
                                                                xy = TRUE))
    recons.fine <- recons.fine[is.na(Z)==F]

    #- merge the two reconstructed ground point clouds
    data <- rbind(recons.all,recons.fine)
  }else{
    data=recons.all #- if no fine reconstruction was needed
  }
  #- add the classification as ground points as in the lidR package
  data[,Classification := 2]
  data.table::setnames(data,c("X","Y","Z","Classification"))

  #- merge the reconstructed ground and the vegetation point cloud
  data <- rbind(data,vegetation)


  data[,Classification := as.integer(Classification)]
  data <- lidR::LAS(data) # export a LAS

  return(data)
}
