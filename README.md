# MPAS-A spectralnug_v1.0

This repository contains a modified version of the **Model for Prediction Across Scales – Atmosphere (MPAS-A, v8.2.2)** that integrates:  
- **analysis nudging (FDDA)** migrated from the implementation of *Bullock Jr. et al. (2018)* for MPAS-A v4.0, and  
- **a diffusion-based spatial filter** inspired by *Grooms et al. (2021)* and implemented consistently with the **[GCM-Filters](https://doi.org/10.21105/joss.03947)** Python package.  

All original MPAS-A license terms are retained (BSD-style license).  
All modifications are clearly marked, and `.original` files are provided for comparison with the unmodified MPAS source.

---

## 1. Nudging scheme (FDDA)

- **Implemented by porting the Bullock Jr. (2018) v4.0 code to MPAS-A v8.2.2.**  
- Introduces non-selective analysis nudging and resolution-scaled nudging tendencies in `physics_tend`.  
- Requires **preprocessed FDDA data**:  
  - Input at current time **plus 6 hours ahead**, since linear interpolation is applied in time.  

**Development record**  
- *2025.06.22* – Initial migration and implementation completed.  
- Every modified file is accompanied by a `.original` backup.  

---

## 2. Diffusion-based spatial filter

- **Approach:**  
  - Converts the GCM-filter convolution kernel (Grooms et al. 2021) into a **multi-step diffusion iteration**.  
  - Performs **MPI halo exchanges** to allow filter scales **much larger than grid spacing** at low communication cost.  
  - Implements **tapered Gaussian kernel approximation**.  
  - Coefficients computed by `diffusion_filter_coefficients.F` using interpolation consistent with the Python GCM-Filters package (Loose et al. 2022).  
  - Supports **anisotropic diffusion** via area weighting.  

**Implementation notes**  
- Requires **OpenMP** to synchronize block exchanges after each diffusion step.  
- An alternate `gaussian_filter` implementation exists for testing, but has high MPI cost.  

**Manual namelist parameters:**  
- `config_fdda_scheme` - including None、analysis，scaled (for setting approximation strength based on grid resolution, see Bullock Jr. et al. 2018)
- `config_fdda_uv` - With the variables involved in the approximation turned on, the wind speed needs to be compiled to the Cell center
- `config_fdda_int` — characteristic timescale (time interval) for nudging target field
- `n_steps` — number of diffusion iterations (computed from a Python demo mimicking gcm-filter)
- `config_filter_shape` - we have two filter shape now: taper and gaussian (see Grooms et al. 2021)
- `config_dx_min` — minimum grid spacing  
- `config_areacellmax` — maximum cell area (manual input to carry out an anisotropic diffusion of norm in areaCell
)  
- `config_gaussian_filter` — performance test only, gauss filtering based on physical distance is performed at each computational block
**keep false in production**  

Example configuration:  
```fortran
&nudging
config_fdda_scheme    = 'analysis'
config_fdda_t         = true
config_fdda_q         = true
config_fdda_uv        = true
config_fdda_uv_in_pbl = false
config_fdda_t_coef    = 3.0e-4
config_fdda_q_coef    = 3.0e-4
config_fdda_uv_coef   = 3.0e-4
config_fdda_uv_min_layer = 5
config_fdda_int       = 21600.
config_nudging_with_spatial_filter = true ! turn on spectral nudging
config_gaussian_filter = false ! test only
config_diffusion_filter = true
config_filter_shape = 'taper'
config_filter_scale = 1000000.0
config_dx_min = 92000.
config_n_steps = 130
config_areacellmax = 9126471680.0 ! maximum areacell for 92-25km mesh
```


---

## 3.License and attribution

- **Base model:** MPAS-A v8.2.2, Copyright (c) 2013 LANS and UCAR, BSD 3-Clause style license.  
- **Nudging reference:** Bullock Jr. et al. (2018), available at [Zenodo](https://zenodo.org/records/1101204) — DOI: [10.5281/zenodo.1101203](https://doi.org/10.5281/zenodo.1101203).  
- **Diffusion filter references:**  
  - Grooms et al. (2021), *Diffusion-Based Smoothers for Spatial Filtering of Gridded Geophysical Data*, [https://doi.org/10.1029/2021MS002552](https://doi.org/10.1029/2021MS002552)  
  - Loose et al. (2022), *GCM-Filters: A Python Package for Diffusion-based Spatial Filtering of Gridded Data*, [https://doi.org/10.21105/joss.03947](https://doi.org/10.21105/joss.03947)  
- **This work:** Migration and extensions by *Yiyuan Cheng (2025)*.  

If you use this modified code, please cite **all** of the following:  
- Bullock Jr., O. R., Foroutan, H., Gilliam, R. C., and Herwehe, J. A. (2018):  
  *Adding four-dimensional data assimilation by analysis nudging to the Model for Prediction Across Scales – Atmosphere (version 4.0)*,  
  Geosci. Model Dev., 11, 2897–2922, [https://doi.org/10.5194/gmd-11-2897-2018](https://doi.org/10.5194/gmd-11-2897-2018)  
- Grooms, I. et al. (2021):  
  *Diffusion-Based Smoothers for Spatial Filtering of Gridded Geophysical Data*,  
  JAMES, 13, e2021MS002552, [https://doi.org/10.1029/2021MS002552](https://doi.org/10.1029/2021MS002552)  
- Loose, N. et al. (2022):  
  *GCM-Filters: A Python Package for Diffusion-based Spatial Filtering of Gridded Data*,  
  JOSS, 7(70), 3947, [https://doi.org/10.21105/joss.03947](https://doi.org/10.21105/joss.03947)  
- MPAS-A model (https://mpas-dev.github.io)  

---

## 4.Notes / Development Status:

yycheng 2025.08
- Due to current skill level, the code may be incomplete and lacks full standardization. Related experiments are ongoing.
- Target fields (fddanew) for nudging still need to be extracted via scripts such as atmosphere_init. Many of these scripts are currently manual and are being organized.
- Some namelist variables are documented in src/core_atmosphere/Registry.xml.
Please compare with the original Registry.xml for differences.