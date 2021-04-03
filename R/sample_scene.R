
#' Recenters and subsets a 3D scene for visibility estimates
#'
#' @description Recenters and if needed subsets a 3D scan image  for use in the
#' \code{\link{visibility}} function. Keeps the points that fall within a user
#' defined distance from the animal location and recenters the scene so that the
#' animal location in the output point has 0,0,0 coordinates. The animal location can
#' be defined by providing xyz coordinates or can be manually selected within the scene.
#' The scene shape can be spherical or circular (see details for more information).
#'
#' @param data LAS class object containing the xyx coordinates of a 3D
#' point cloud.
#' @param scene_radius numerical. The radius of the final scene. Can refer to
#' the radius of a sphere if \code{scene_shape = "sph"} or of a circle
#' if \code{scene_shape = "circ"}.
#' @param scene_shape character string. Defines the shape of the scene: "sph" and
#' "circ" are accepted (see details for more informations). Default = "circ".
#' @param center (optional) vector of length 3 providing the xyz coordinates of
#' the user defined animal location. If not provided, the user can manually select
#' the animal location in the 3D point cloud. The average coordinates of the selected
#' region will be set as the animal location (see details).
#'
#' @param downsample numeric. Enables the user to downsample the point cloud before
#' visualizing it for scene center manual selection (if no \code{center} is
#' provided). Defines the voxel resolution within which a single point of the
#' input scene will be kept.
#' \code{downsample = 0} desable downsampling. Default is 0 if the scene contains less
#' than 5e6 surveyed points or 0.1 if the scene contains more than 5e6 surveyed points.
#' @param messages logical. Disables the messages and message box when manually
#' selecting the scene center.
#'
#' @details
#' \subsection{Scene shape}{if \code{scene_shape = "circ"} the distance
#' to scene center is computed in the xy dimension of the original scene only,
#' resulting in a circular scene. If \code{scene_shape = "sph"}
#' the distance to the scene center is similar to arguments \code{scene.radius}
#' in \code{\link{visibility}} function and \code{cut_off} in
#' \code{\link{viewsheds}}.}
#' \subsection{Manual selection of scene center}{ if no \code{center} is
#' provided, a 3D plot automatically opens. The user can navigate around the
#' surveyed points (rotate = left
#' click, pan = right click) and select points (in a rectangular region)
#' with the middle click. Once the points are selected, a message box opens
#' (if not disabled). If the user clicks "yes", the plot window closes and the
#' average point coordinates are defined as the animal location. If the user clicks
#' "no", he/she can revise the previous selection.}
#'
#' @return A LAS class object containing the coordinates of the reshaped scene.
#' @export
#'
#' @examples
#'
#' #- import the tree_line_plot dataset
#' file <- system.file("extdata", "tree_line_plot.laz", package="viewshed3d")
#' tls  <-  lidR::readLAS(file,select="xyz")
#'
#' #- define the animal location
#' center <- c(0,-6,1)
#'
#' #- reshape the TLS scene with scene_shape="circ" and the calculated center
#' reshaped <- viewshed3d::sample_scene(tls,scene_radius = 4,
#'                                      center=center,
#'                                      scene_shape = "circ")
#'
#' lidR::plot(reshaped)
#'
#' #- reshape the TLS scene with scene_shape="sph" and the calculated center
#' reshaped <- viewshed3d::sample_scene(tls,scene_radius = 4,
#'                                      center=center,
#'                                      scene_shape = "sph")
#'
#' lidR::plot(reshaped)
#'
#' if (interactive()){
#' #- manual selection of the center
#' reshaped <- viewshed3d::sample_scene(tls,scene_radius = 4,
#'                                      scene_shape = "circ")
#'
#' lidR::plot(reshaped)
#' }
#'
sample_scene  <-  function(data,scene_radius,scene_shape,center,
                           downsample,messages){

  #- declare variables to pass CRAN check as suggested by data.table mainaitners
  r=X=Y=Z=NULL

  #- defaults parameters
  if(missing(scene_radius)) stop("scene_radius is missing.") #- no default
  if(missing(scene_shape)) scene_shape <- "circ"

  #- preventing some miss uses
  if(missing(data)) stop("no input data provided.")
  if(class(data)[1]!="LAS"){
   stop("data must be a LAS. You can use lidr::LAS(data) to convert a data.frame
   or data.table to LAS.")
  }
  if(is.numeric(scene_radius)==F | length(scene_radius)!=1){
    stop("scene radius must be numeric of length 1.")
  }
  if((scene_shape=="sph" | scene_shape=="circ")==F){
    stop("scene_shape must be 'sph' or 'circ'.")
  }

  #- if no scene center is provided -> manual selection in the 3d plot
  if(missing(center)){
    if(missing(messages)) messages=T  #- default message_box parameter
    if(missing(downsample)){
      if(nrow(data@data)>5e6){downsample=0.1}else{downsample=0}
    }
    #- if down sampling is required apply it, if not display the original scene
    if(downsample>0){
      if(is.numeric(downsample)==F | length(downsample)!= 1){
        stop("donwsample must be a numeric of length 1.")
      }

      #- down sampling with TreeLS tools
      down <- downsample_scene(data,method = "space", filter = 0.03)
      x <- lidR::plot(down,mapview=F) #- display donw sampled point cloud
    }else{
      x <- lidR::plot(data,mapview=F) #- display original scene
    }

    #- while the desired region was not selected -> continue to select the
    #- point region to calculate the scene center
    run=T
    while(run==T){
      if(messages==T){
        print("Please select a region of the point cloud (middle cick)",quote=F)
      }
      #- select a point region with middle button
      center <- rgl::selectpoints3d(value = T, closest = F,
                                    multiple = F,button="middle")
      #- show selected points in red
      pts <- rgl::plot3d(center,col="red",size=10,add=T)

      if(messages==T){
        #- open a box to confirm it is the right selected region
        rep <- tcltk::tkmessageBox(message = "Is it the right location ?",
                                   icon="question", type = "yesno",
                                   default = "yes")
      }else{
        #- if message box is desabled, the selected region is considered as good
        rep <- "yes"
      }
      if(as.character(rep)=="yes"){
        rgl::rgl.close()
        center <- c(mean(center[,1]),mean(center[,2]),mean(center[,3]))
        run <- F
      }else{rgl::rgl.pop(id=pts)}
    }

    #-transfrom data into a data.table object to apply the distance based filter
    data=data@data

    #- correct the calculated point center: when applying rgl::selectpoint3d to
    #- a lidR::plot the selected points coordinates have to be corrected based
    #- on minimum x and y coordinates
    center=c(center[1]+x[1],center[2]+x[2],center[3])
  }else{
    #-transfrom data into a data.table object to apply the distance based filter
    data <- data@data
  }

  if(scene_shape=="sph"){
    #- distance to a sphere center
    data[,r := sqrt(((X-center[1])^2)+((Y-center[2])^2)+((Z-center[3])^2))]
  }
  if(scene_shape=="circ"){
    #- distance to a disc center
    data[,r := sqrt(((X-center[1])^2)+((Y-center[2])^2))]
  }

  data=data[r<=scene_radius] #- select the points
  data[,r := NULL] #- remove the radius column

  #- recenter the scene
  data[, ':=' (X=X-center[1],
               Y=Y-center[2],
               Z=Z-center[3])]

  data <- pkgcond::suppress_messages(lidR::LAS(data)) # export a LAS

  return(data)
}
