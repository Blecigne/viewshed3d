#'Tools to compute visibility in 3D point clouds of ecosystems
#'
#'For many animals, the ability to visually assess the environment and detect approaching predators is an important
#'part of anti-predator strategies. Because this can occur across spatial scales, estimation of the viewshed can help to quantify
#'visibility as a continuous variable around animal locations and facilitate studies of habitat selection and predator-prey interactions.
#'
#'\subsection{Visibility and cumulated viewsheds}{
#'visibility within a spherical viewshed is calculated using the \code{\link{visibility}}
#'function. This function is designed to sample the point cloud in every direction of the 3D space from
#'a single user-defined location and to record the distance to the nearest point in each direction. Each direction is thus considered as a
#'sightline - of a user defined angle - that is assumed to end when an object is encountered. The \code{\link{d_visibility}} function (for directional visibility) provides
#'similar capabilities than the \code{\link{visibility}} function but enable to segment the viewshed to analyse the visibility in a particular direction.
#'The \code{\link{h_visibility}} function (for horizontal visibility) computes the visibility in a viewshed defined as a disk. This enable to compute the extend of animal visibility
#'parallel to the ground.
#'The \code{\link{viewsheds}} function computes the overlap between viewsheds calculated from different locations and returns a voxel cloud
#'quantifying for each voxel (i.e. each portion of the 3D scene) the number of times it was visible from any location.}
#'
#'\subsection{Ground reconstruction}{
#'in the point clouds, some portions of the ground is frequently not sampled by the sensor (especially in the case of a TLS).
#'That would result in infinite sightlines that continue
#'below the ground surface. To correct for this effect, the \code{\link{reconstruct_ground}} function can be used to reconstruct the ground
#'before using the \code{\link{visibility}} function. The \code{\link{reconstruct_ground}} function computes the optimal resolution to
#'reconstruct the ground based on user-defined parameters for visibility calculation.}
#'
#'\subsection{3D scene reshaping}{
#'because 3D scenes might cover a large area, but the visibility analyses might be computed for smaller areas, the \code{\link{sample_scene}} function can be used
#'at the beginning of the data preparation process to segment a scene, with the appropriate properties in terms of size and shape for visibility
#'calculation. This might be usefull to reduce computation time during the ground reconstruction process.}
#'
#'  \subsection{Noise filters}{
#'the \code{\link{denoise_scene}} function provides three different methods to filter isolated points from 3D point clouds.
#'The \code{\link{remove_slope}} function enable to remove the slope from a 3D scene and the \code{\link{add_slope}} function reverse this operation. }
#'
#' \subsection{Other functions}{
#'the \code{\link{generate_sphere}} function can be used to generate a spherical point cloud with a user defined angular resolution.
#'The \code{\link{set_position}} function computes the correct elevation for the animal position based on the ground elevation. }
#'
#'\subsection{Dataset}{the \code{viewshed3d} package provides a TLS scene of a circular forest plot located at
#'northern treeline sites in Alaska (\code{tree_line_plot.laz}). This dataset has the following specifications :
#'\itemize{
#' \item Format: LAS
#' \item 2513044 points
#' \item radius: 17 m
#' \item center coordinates: 0,0,0
#' \item duplicated points removed
#' \item dowsampled by keeping one point within a 2 cm voxel
#'}}
#'
#' @docType package
#' @name viewshed3d
NULL
