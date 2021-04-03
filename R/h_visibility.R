
#' Compute horizontal visibility
#'
#' @description Computes the horizontal visibility in a viewshed defined as a disk.
#'
#' @param data LAS class object containing the xyz coordinates of a 3D point
#' cloud.
#' @param position vector of length 3 containing the xyz coordinates of the
#' animal location. Default = c(0,0,0).
#' @param layer_tickness numeric. The thickness of the disk that defines the viewshed.
#' @param angular_res numeric. The angular resolution of a single sightline.
#' Default = 1.
#' @param scene_radius (optional) numeric. Defines the radius of the scene
#' relative to the animal position. Can be used to apply a cut-off distance to
#' visibility analyses.
#'
#' @return A list containing a data.table of the visibility
#' (\code{$visibility}) as a function of distance to the animal location (r), a data.table of
#' the viewshed (\code{$viewshed}) defined as the radius (r) for each azimuth (phi) with non occluded
#' sightlines distance set to scene_radius, a data.table (\code{$vegetation_distance}) of the vegetation distance (r)
#' in each azimuth (phi) for occluded sightlines only and different viewshed metrics (\code{$metrics}).
#' The metrics are: the viewshed fractal dimension (computed from Chen, 2020), the
#' viewshed area (i.e. the visible area), the proportion visible (the area visible
#'  / potential area visible), the viewshed coefficient (the area under the curve of
#' visibility as a function of distance) and the vegetation fractal dimension.
#'
#' @references
#' Chen, Y. (2020). Two Sets of Simple Formulae to Estimating Fractal Dimension
#' of Irregular Boundaries. Mathematical Problems in Engineering, 2020.
#' https://doi.org/10.1155/2020/7528703
#'
#' @export
#'
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
#' #- compute the horizontal visibility.
#' view.data <- viewshed3d::h_visibility(data = tls_clean,
#'                                       position = center,
#'                                       angular_res = angle,
#'                                       scene_radius = 17)
#'
#' # viewshed metrics
#' view.data$metrics
#'
#' # plot the view sheds in a radial plot
#' plotrix::radial.plot(view.data$viewshed$r,rp.type = "p",poly.col = "blue",
#'                      radial.lim = c(0,max(view.data$viewshed$r)))
#' }
h_visibility = function(data,position = c(0,0,0),layer_tickness = 0.1,angular_res = 1,scene_radius){

  #- declare variables to pass CRAN check as suggested by data.table mainaitners
  X=Y=Z=phi=r=.=Visibility=NULL

  data=data@data

  #- translate the data so the center is 0,0,0
  data[, ':=' (X = X - position[1],
               Y = Y - position[2],
               Z = Z - position[3])]

  # select the layer of interest
  data <- data[Z >= -layer_tickness/2 & Z <= layer_tickness/2]

  # compute azimuth and distance from center
  data[, ':=' (phi = round((acos(Y * 100/(sqrt(Y^2 + X^2) * 100)) * (180/pi))/angular_res)*angular_res,
               r = sqrt(X^2+Y^2+Z^2))]
  data[ X < 0 , phi := (180-phi)+180] # transform phi so its range is now 0-360

  # keep the data within the scene radius
  data <- data[r <= scene_radius]

  # build a table with all directions and the scene radius as distance
  radius_tab <- data.table::data.table(phi = seq(0,360,angular_res),r = scene_radius)

  if(nrow(data)>0){
    # compute the minimum distance in each phi
    data_radius <- data[,min(r),by = phi]
    data.table::setnames(data_radius,old="V1",new="r")

    # update radius_temp with the measured distances
    data.table::setkey(radius_tab,phi)
    data.table::setkey(data_radius,phi)
    radius_tab[data_radius,r := data_radius$r]
  }else{
    data_radius <- data.table::data.table(phi = radius_tab$phi,r = 0)
  }

  ########### computing area and fractal dimension of the viewshed
  # transform the polar coordinates into x y coordinates
  xy <- data.table::data.table(
    p = pracma::deg2rad(radius_tab$phi),
    t = 0,
    r = radius_tab$r
  )
  xy <- pracma::sph2cart(as.matrix(xy))

  # compute the area and perimeters of the polygon
  area_vs <- pracma::polyarea(x=xy[,1],y=xy[,2])
  perimeter <- pracma::poly_length(x=xy[,1],y=xy[,2])

  # fractal dimension from the formula presented in eq. 20 in Chen 2020 Chen, Y.
  # (2020). Two Sets of Simple Formulae to Estimating Fractal Dimension of
  #Irregular Boundaries. Mathematical Problems in Engineering.
  # Here I assume that the shape deviates from a circle
  FD_vs <- (2*log(perimeter/(2*pi)))/log(area_vs/pi)

  ########### computing area and fractal dimension of the vegetation
  # transform the polar coordinates into x y coordinates
  xy <- data.table::data.table(
    p = pracma::deg2rad(data_radius$phi),
    t = 0,
    r = data_radius$r
  )
  xy <- pracma::sph2cart(as.matrix(xy))

  # compute the area and perimeters of the polygon
  area <- pracma::polyarea(x=xy[,1],y=xy[,2])
  perimeter <- pracma::poly_length(x=xy[,1],y=xy[,2])

  # fractal dimension
  FD_veg <- (2*log(perimeter/(2*pi)))/log(area/pi)

  # compute visibility as a function of distance
  radius_tab[,visibility := 1]
  visibility <- radius_tab[order(r),.(r,visibility)]
  visibility[, visibility := (1-cumsum(visibility)/(nrow(radius_tab)))*100]

  # add the range 0 to minimum distance
  if(min(visibility$r) > 0){
    visibility <- data.table::rbindlist( list(
      data.table::data.table(r= seq(0, min(visibility$r),0.1),visibility = 100),
      visibility
    ))
  }

  # add the range maximum distance to scene radius
  if(max(visibility$r) < scene_radius){
    visibility <- data.table::rbindlist( list(
      visibility,
      data.table::data.table(r= seq(max(visibility$r), scene_radius,0.1),visibility = min(visibility$visibility))
    ))
  }

  radius_tab[,visibility := NULL]
  # return the datasets and metrics
  return(list(
    visibility = visibility,
    viewshed = radius_tab,
    vegetation_distance = data_radius,
    metrics = list(
      viewshed_fractal_dimension = FD_vs,
      viewshed_area = area_vs,
      proportion_visible = area_vs/(pi*scene_radius^2),
      viewshed_coeffecient = pracma::trapz(x = visibility$r, y = visibility$visibility),
      vegetation_fractal_dimension = FD_veg
      )
  ))
}
