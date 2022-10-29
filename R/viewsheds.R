#' Computes cumulated viewsheds
#'
#' @description Computes cumulated viewsheds within a 3D point cloud and return
#' a voxel cloud accounting for the number of times each voxel was visible.
#'
#' @param data LAS class object containing the xyx coordinates of a 3D point
#' cloud.
#' @param positions data.frame or data.table with 3 columns containing the xyz
#' coordinates of the animal locations from which the viewsheds will be computed.
#' @param angular_res numeric. The angular resolution of a single sightline.
#' Default = 1.
#' @param vox_res numeric. The resolution of the output voxel cloud. Default =
#' 0.2.
#' @param cut_off (optional) numeric. Defines a cut-off distance for each
#' individual viewshed. Speeds up the process
#' when \code{viewsheds} is applied to big datasets.
#' @param pb logical. If \code{FALSE}, desables the progress bar.
#'
#' @return A LAS class object containing the coordinates of the voxel cloud
#' (X, Y, Z), and the number of times each voxel was visible from any position
#' (N_visible).
#'
#' @note In most cases, a ground reconstruction should be performed before
#' viewsheds computation. This can be done with
#' the \code{\link[lidR]{classify_ground}} and \code{\link[lidR]{grid_terrain}}
#' functions from the \code{\link[lidR]{lidR-package}}.
#'
#' @details Sightline directions in each viewshed are computed from the method
#' described by Malkin (2016). This ensures that every sightline explores
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
#' file = system.file("extdata", "tree_line_plot.laz", package="viewshed3d")
#' tls = lidR::readLAS(file,select="xyz")
#'
#' #- remove noise to avoid visibility estimates error
#' tls_clean = viewshed3d::denoise_scene(tls,method="sd",
#'                                        filter=6)
#'
#' #- RECONSTRUCT THE GROUND
#' #- classify ground points
#' class=lidR::classify_ground(tls_clean, lidR::csf(rigidness = 1L,
#'                                            class_threshold = 0.1,
#'                                            sloop_smooth = TRUE), FALSE)
#'
#' #- reconstruct the ground. No need for a very fine ground reconstruction.
#' ground = lidR::grid_terrain(class, 0.05, lidR::knnidw())
#'
#' #- build the final scene
#' reconstructed = na.omit(raster::as.data.frame(ground, xy = TRUE))
#' names(reconstructed)=c("X","Y","Z")
#' recons=rbind(lidR::LAS(na.omit(reconstructed)),tls_clean)
#'
#' #- CREATE THE POSITIONS WITH RANDOM POINTS
#' N_positions = 10 #- how many points ?
#' height = 2 #- points height relative to the ground
#' positions=data.table::data.table(reconstructed[runif(N_positions,
#'                                                1,nrow(reconstructed)),])
#' positions[,Z:=Z+height]
#'
#' #- compute the cumulated viewsheds from the positions
#' cumulated=viewshed3d::viewsheds(data=recons,
#'                                 positions = positions ,
#'                                 angular_res = 1,
#'                                 vox_res = 0.2)
#'
#' #- plot the result
#' x=lidR::plot(cumulated,color="N_visible",size=3,
#'              colorPalette=viridis::cividis(nrow(positions)+1),trim=6)
#'
#' #- add the positions
#' lidR::add_treetops3d(x,sp::SpatialPointsDataFrame(positions,positions),
#'                      radius=0.5,col="red",add=TRUE)
#'}
viewsheds=function(data,positions,angular_res,vox_res,cut_off,pb){

  #- declare variables to pass CRAN check as suggested by data.table mainaitners
  x=y=z=X=Y=Z=theta=phi=r=.=N=N_visible=NULL

  #- default parameters
  if(missing(pb)) pb <- T
  if(missing(angular_res)) angular_res <- 1
  if(missing(vox_res)) vox_res <- 0.2

  #- check for possible problems in arguments
  if(missing(data)) stop("no input data provided.")
  if(class(data)[1]!="LAS"){
    stop("data must be a LAS. You can use lidr::LAS(data) to convert a
         data.frame or data.table to LAS.")
  }
  if(missing(positions)){
    stop("No positions provided")
  }else{
    if(class(positions)[1] %in% c("data.frame","data.table")==F){
      stop("positions must be a data.frame or data.table.")
    }
    if(ncol(positions)!=3){stop("positions must have 3 columns: X, Y, Z.")}
    #- transform positions into a data.table
    positions <- data.table::data.table(positions)
    data.table::setnames(positions,c("X", "Y", "Z"))
  }
  if(is.numeric(angular_res)==F) stop("angular_res must be numeric.")
  if(is.numeric(vox_res)==F) stop("vox_res must be numeric.")
  if(is.logical(pb)==F) stop("pb must be logical")

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

  #- set data names
  data <- data@data[,.(X,Y,Z)] # transform data into a data.table
  data.table::setnames(data,c("x", "y", "z"))

  #- BUILD THE VOXEL CLOUD
  vox_dat <- data[,":=" (X=round(x/vox_res)*vox_res,
                         Y=round(y/vox_res)*vox_res,
                         Z=round(z/vox_res)*vox_res)]
  vox_dat <- unique(vox_dat[,.(X,Y,Z)])
  #- add the N field to store the number of time each voxel was visible
  vox_dat[,N_visible := 0]
  data.table::setkeyv(vox_dat,c("X", "Y", "Z")) #- for matchin with "points"

  data <- data[,1:3]

  #- COMPUTE VISIBILITY from each position
  if(pb==T) prog <- utils::txtProgressBar(min=0,max=nrow(positions),style = 3,width = 33)
  for(i in 1:nrow(positions)){
    if(pb==T) utils::setTxtProgressBar(prog, i)

    #- translate the data so the center is 0,0,0
    data[, ':=' (X = x - positions[i,X],
                 Y = y - positions[i,Y],
                 Z = z - positions[i,Z])]
    dat <- data[,.(X,Y,Z)] #-duplicate data to keep intact for next iteration

    dat[,r := sqrt(X^2+Y^2+Z^2)]

    if(missing(cut_off)==F){ #- apply the scene radius if defined
      #- ensure scene_radius is valid
      if(is.numeric(cut_off)==F) stop("cut_off must be numeric.")
      dat <- dat[r <= cut_off]
      if(nrow(dat)==0){
        stop(paste("scene.radius is too small, the scene contains 0 point
                   at iteration ",i,". Increase scene_radius or change position "
                   ,i,".",sep=""))
      }
    }

    #- compute elevation angle (theta) and azimuth (phi)
    dat[, ':=' (theta = round((acos(Z * 100/(sqrt(X^2 + Y^2 + Z^2) * 100)) *
                                 (180/pi))/angular_res)*angular_res,
                phi = acos(Y * 100/(sqrt(Y^2 + X^2) * 100)) * (180/pi))]

    dat[ X < 0 , phi := (180-phi)+180] # transform phi so its range is now 0-360

    #- match param and data
    data.table::setkey(param,theta)
    data.table::setkey(dat,theta)
    dat[, N:= param[dat,N]]

    #- compute the phi for each cell center
    dat[, phi := round(phi/(360/N))*(360/N)]
    dat[ phi>=360 , phi:=0]

    #- find the nearest point in each direction and store xyz coordinates
    near <- dat[, .(r = min(r)), keyby = .(theta,phi)]
    near <- unique(near)

    #- subset visible points
    data.table::setkeyv(dat,c("theta","phi","r"))
    data.table::setkeyv(near,c("theta","phi","r"))
    dat <- dat[near, .(X,Y,Z)]

    #- create the voxel cloud of visible points
    dat[, ':=' (X = X + positions[i,X],
                Y = Y + positions[i,Y],
                Z = Z + positions[i,Z])]

    dat <- unique(dat[,":=" (X=round(X/vox_res)*vox_res,
                             Y=round(Y/vox_res)*vox_res,
                             Z=round(Z/vox_res)*vox_res)])

    #- add 1 in the "N" field of the voxel cloud of the scene
    data.table::setkeyv(dat,c("X","Y","Z")) #- to match with vox_dat
    vox_dat[dat,N_visible := N_visible+1] #- add 1 to each joint voxels
  }

  vox_dat <- pkgcond::suppress_messages(lidR::LAS(vox_dat, check = FALSE)) # export a LAS
  return(vox_dat)
}
