# MPAS-A spectralnug v1.1

2026.06.04 Yiyuan Cheng
 
**Base model:** MPAS-A v8.2.2  


This repository provides a modified MPAS-A v8.2.2 patch for **analysis nudging (FDDA)** and **diffusion-based scale-selective nudging**. The nudging implementation was migrated from the MPAS-A v4.0 FDDA framework of Bullock Jr. et al. (2018), while the spatial filtering follows the diffusion-based filtering idea of Grooms et al. (2021) and the GCM-Filters implementation of Loose et al. (2022).

All original MPAS-A license terms are retained. This repository only contains the files needed to patch a clean MPAS-A v8.2.2 source tree.


## 1. Installation as a minimal patch

this repository is organized as a **minimal patch** rather than a complete modified MPAS-A source tree. After preparing a clean MPAS-A v8.2.2 source code, copy the contents of this repository into the MPAS-A directory and overwrite the corresponding files.

```bash
# Example only: run this from the patch repository
cp -r ./* /path/to/MPAS-Model/
```

Then recompile MPAS-A as usual.


## 2. Main features

**FDDA analysis nudging**

- Migrated from the MPAS-A v4.0 FDDA implementation of Bullock Jr. et al. (2018) to MPAS-A v8.2.2.
- Adds analysis nudging tendencies in `physics_tend`.
- Supports nudging for wind, temperature, and moisture fields.
- Requires preprocessed FDDA target data. The target field at the current time and the next forcing time is needed because the nudging target is linearly interpolated in time.

For standard grid/analysis nudging, use:

```fortran
config_fdda_scheme = 'analysis'
```

The previous `scaled` option should no longer be selected through `config_fdda_scheme`; in v1.1 it is handled together with the regional mask logic.

**Diffusion-based spatial filter**

- Converts the filter kernel into repeated diffusion iterations.
- Uses MPI halo exchanges, allowing filter scales much larger than the local grid spacing at relatively low communication cost.
- Supports `taper` and `gaussian` diffusion-filter shapes.
- Coefficients are computed following the GCM-Filters-style diffusion approximation.
- OpenMP synchronization is required around block/halo exchanges during the filtering iterations.

To enable scale-selective nudging with the diffusion filter:

```fortran
config_nudging_with_spatial_filter = true
config_gaussian_filter             = false  ! deprecated direct distance-convolution filter
config_diffusion_filter            = true
```

Example:

```fortran
config_filter_shape = 'taper'       ! or 'gaussian'
config_filter_scale = 1000000.0     ! m
config_n_steps      = 43
```


## 3. Diffusion type in v1.1

Anisotropic diffusion

Recommended for variable-resolution meshes. The diffusion coefficient is scaled by `areaCell`, so the effective filtering scale varies with mesh resolution. This helps maintain stability when applying diffusion filtering on variable-resolution grids.

```fortran
config_fdda_diffusion_anisotropic = true
config_areacellmax = 9126471680.0   ! m^2, maximum areaCell for the 92-25 km mesh
```

Isotropic diffusion

```fortran
config_fdda_diffusion_anisotropic = false
```

The isotropic option is available, but it has not yet been tested in long-term integrations.


## 4. Nudging variables and vertical layers

Basic nudging switches follow the MPAS-A FDDA-style namelist structure, for example:

```fortran
config_fdda_t         = true
config_fdda_q         = true
config_fdda_uv        = true
config_fdda_uv_in_pbl = false
config_fdda_t_coef    = 3.0e-4
config_fdda_q_coef    = 3.0e-4
config_fdda_uv_coef   = 3.0e-4
config_fdda_int       = 21600.
```

Since v1.1, minimum and maximum vertical layers can be specified for each nudged variable. For example:

```fortran
config_fdda_uv_min_layer = 5
config_fdda_uv_max_layer = 56
```

Use the corresponding `t`, `q`, and `uv` layer options as needed. Detailed namelist definitions are in:

```text
src/core_atmosphere/Registry.xml
```


## 5. Regional-mask nudging

v1.1 adds regional-mask nudging. Outside the selected region, the nudging tendency is set to zero. The previous resolution-dependent `scaled` logic is also handled in this mask framework.

Resolution-scaled range

```fortran
config_fdda_scale_min = 0.0
config_fdda_scale_max = 100000.0
```

Circular mask

```fortran
config_fdda_mask_lon0   = 110.0
config_fdda_mask_lat0   = 30.0
config_fdda_mask_radius = 1500000.0
```

Elliptical mask

```fortran
config_fdda_mask_lon0      = 110.0
config_fdda_mask_lat0      = 30.0
config_fdda_mask_radius_x  = 1800000.0
config_fdda_mask_radius_y  = 900000.0
config_fdda_mask_rotation  = 30.0   ! degree, counterclockwise from local east
```

Latitude-longitude box

```fortran
config_fdda_mask_lon_min = 90.0
config_fdda_mask_lon_max = 130.0   ! if smaller than lon_min, the box crosses the dateline
config_fdda_mask_lat_min = 15.0
config_fdda_mask_lat_max = 45.0
```

## 6. Example namelist block

for VR92-25km mesh with 56 vertical layers

```fortran
&nudging
config_fdda_scheme      = 'analysis'
config_fdda_t           = true
config_fdda_q           = true
config_fdda_uv          = true
config_fdda_uv_in_pbl   = false
config_fdda_t_coef      = 3.0e-4
config_fdda_q_coef      = 3.0e-4
config_fdda_uv_coef     = 3.0e-4
config_fdda_int         = 21600. ! The nudging files are updated every 6 hours.

config_fdda_uv_min_layer = 5
config_fdda_uv_max_layer = 56

config_nudging_with_spatial_filter = true
config_gaussian_filter             = false
config_diffusion_filter            = true
config_filter_shape                = 'taper'
config_filter_scale                = 1000000.0
config_n_steps                     = 43

config_fdda_diffusion_anisotropic  = true
config_areacellmax                 = 9126471680.0
/
```

The `streams.atmosphere` file should include:

```
<immutable_stream name="fdda"
                  type="input"
                  io_type="pnetcdf,cdf5"
                  filename_template="/path/to/your/nudging/target/FDDA_VR92-25km.$Y-$M-$D_$h.$m.$s.nc"
                  filename_interval="input_interval"
                  input_interval="6:00:00" />
```

**For a 6-hour nudging interval, the target file read at 00:00:00 should contain the 06:00:00 target field. The nudging tendency is then linearly interpolated between 00:00:00 and 06:00:00. This is not an ideal design.**

## 7. Notes

- `config_gaussian_filter` refers to an older direct physical-distance convolution filter. It is retained only for testing and should normally be kept `false`.
- The diffusion filter is the recommended spatial-filter option.
- Some scripts for preparing FDDA target fields may still be workflow-dependent.
- Please check `src/core_atmosphere/Registry.xml` for the complete and latest namelist definitions and variable input requirements.


## 8. Attribution

- **Base model:** MPAS-A v8.2.2, Copyright (c) LANS and UCAR, BSD 3-Clause style license.
- **Nudging reference:** Bullock Jr. et al. (2018), MPAS-A v4.0 analysis nudging.
- **Diffusion filter references:** Grooms et al. (2021) and Loose et al. (2022).
- **This work:** Migration and extensions by Yiyuan Cheng.

Please cite the relevant original works when using this code:

- Bullock Jr., O. R., Foroutan, H., Gilliam, R. C., and Herwehe, J. A. (2018): *Adding four-dimensional data assimilation by analysis nudging to the Model for Prediction Across Scales – Atmosphere (version 4.0)*, Geosci. Model Dev., 11, 2897–2922. https://doi.org/10.5194/gmd-11-2897-2018
- Grooms, I. et al. (2021): *Diffusion-Based Smoothers for Spatial Filtering of Gridded Geophysical Data*, JAMES, 13, e2021MS002552. https://doi.org/10.1029/2021MS002552
- Loose, N. et al. (2022): *GCM-Filters: A Python Package for Diffusion-based Spatial Filtering of Gridded Data*, JOSS, 7(70), 3947. https://doi.org/10.21105/joss.03947
- MPAS-A model: https://mpas-dev.github.io
