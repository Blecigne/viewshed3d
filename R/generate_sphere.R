#' Generate a spherical point cloud
#'
#' @param angular_res the angular resolution that sets the spacing between points.
#' @param r numeric. The sphere radius.
#' @param sph logical. If \code{TRUE} the spherical coordinates of the points
#' are returned, if \code{FALSE} the XYZ coordinates are returned.
#'
#' @return a data.table containing points coordinates.
#' @export
#'
#' @examples
#' # generate a sphere
#' sph <- viewshed3d::generate_sphere()
#'
#' # plot the sphere
#' rgl::plot3d(sph)
#'
generate_sphere = function(angular_res = 1,r = 1,sph = FALSE){

  if(!is.numeric(angular_res)) stop("angular_res must be numeric")
  if(!is.numeric(r)) stop("r must be numeric")
  if(!is.logical(sph)) stop("sph must be logical")

  #- declare variables to pass CRAN check as suggested by data.table mainaitners
  N=theta=phi=.=NULL

  # convert angular_res to radients
  angular_res <- pracma::deg2rad(angular_res)

  #- with theta = the elevation angle, and N = the number of cells
  param <- data.table::data.table(theta=seq(-pi/2,pi/2,angular_res),N=NA)

  # compute the number of cells in each latitudinal band following the
  # method described by Malkin (2016)
  L <- pi/(pi/angular_res) # cell side length at the equator
  N_cell <- round(pi/L) # number of cells at the equator
  # the number of decreases with increasing the elevation angle
  param[,N:=round(N_cell*(cos(theta)))]

  # N = 1 for the two poles
  #param[theta==pi/2 | theta==-pi/2,N:=1]
  param[,N := N+1]

  # compute the azimuth (phi) for each theta
  # replicate N times each theta
  param <- param[rep(1:nrow(param),param$N)]

  # compute phi in each row
  param[,phi := 1]
  param[,phi := (cumsum(phi)*(2*(pi+angular_res)/N))-(2*angular_res),by=theta]

  # convert polar coordinates into xyz point cloud

  if(sph){
    return(param[,.(phi,theta,r)])
  }else{
    xyz <- pracma::sph2cart(
      as.matrix(
        data.table::data.table(
          p = param$phi,
          t = param$theta,
          r = r
        )
      )
    )

    return(
      data.table::data.table(
        X = xyz[,1],
        Y = xyz[,2],
        Z = xyz[,3]
      )
    )
  }
}
