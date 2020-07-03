# viewshed3d
## Tools to compute visibility in 3D point clouds of ecosystems

For many animals, the ability to visually assess the environment and detect approaching predators is an important
part of anti-predator strategies. Because this can occur across spatial scales, estimation of the viewshed can help to quantify
visibility as a continuous variable around animal locations and facilitate studies of habitat selection and predator-prey interactions.

## Related paper
More informations on the method can be found in:

Lecigne, B., Eitel, J. U., & Rachlow, J. L. (2020). viewshed3d: An r package for quantifying 3D visibility using terrestrial lidar data. Methods in Ecology and Evolution, 11(6), 733-738.
https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13385

## Visibility and cumulated viewsheds

Visibility within a single viewshed is calculated using the visibility() function. This function is designed to sample the point cloud in every direction of the 3D space from
a single user-defined location and to record the distance to the nearest point in each direction. Each direction is thus considered as a
sightline - of a user defined angle - that is assumed to end when an object is encountered.
The viewsheds() function computes the overlap between viewsheds calculated from different locations and returns a voxel cloud
quantifying for each voxel (i.e. each portion of the 3D scene) the number of times it was visible from any location.

## Ground reconstruction

In the point clouds, some portions of the ground is frequently not sampled by the sensor (especially in the case of a TLS).
That would result in infinite sightlines that continue
below the ground surface. To correct for this effect, the reconstruct_ground() function can be used to reconstruct the ground
before using the visibility() function. The reconstruct_ground() function computes the optimal resolution to
reconstruct the ground based on user-defined parameters for visibility calculation.

## 3D scene reshaping

Because 3D scenes might cover a large area, but the visibility analyses might be computed for smaller areas, the sample_scene() function can be used
at the beginning of the data preparation process to segment a scene, with the appropriate properties in terms of size and shape for visibility
calculation. This might be usefull to reduce computation time during the ground reconstruction process. 

## Noise and point cloud density filters

The downsample_scene() function can be use to reduce the point cloud density and the denoise_scene() function provides three different methods to filter isolated points from 3D point clouds.

## Workflow to compute visibility and cumulative viewsheds.

<p align="center">
<img src="https://github.com/Blecigne/viewshed3d/blob/master/Workflow_figure.png" width="500">
</p>

The preprocessing steps are indicated on the left, and the
functions used at each step are indicated in the centre. External functions from the lidR package are indicated by ‘lidR::’. Illustrations show (from top to
bottom) (i) the reshaped input terrestrial laser scanning scene after duplicated points were removed, (ii) the scene
with noise shown in red, (iii) the ground reconstructed with the reconstruct_ground() function (note the portion of the
ground reconstructed with a finer point resolution around the scene centre, on left, or the grid_terrain() function, on right)
and the reconstructed scene (middle) and (iv) the output of the visibility() and viewsheds() functions. On the left, 3D view
of visible points from a single location (red sphere) shown in black and non-visible points in grey, and the associated visibility
as a function of the distance to location. On the right, the cumulative viewshed computed from 10 locations is shown with
a colour gradient: darker colour indicates less visible voxels.
