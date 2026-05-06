# Using Ocean Bottom Pressure Measurements of Internal Waves in the South China Sea to Track the Annual Cycle of Stratification
_Amanda Syamsul, Emily E. Brodsky, Ethan F. Williams, Heather R. Crume, Kristen A. Davis, Shiou-Ya Wang, Shu-Kun Hsu_

## Research summary
- Measurements of internal solitary wave period reveal variations in local stratification and its seasonal cycle.
- Ocean-bottom pressure under large-amplitude internal solitary waves is dominated by non-hydrostatic effects. 
- Geostationary satellite imagery provides decadal-scale measurements of internal wave speed and arrival direction.


## Contents
**1. OBS_himawari_analysis.m**
- Analyzes nonlinear internal wave propagation near the Dongsha OBS using satellite detections, bathymetry, and tidal models.
- Quantifies relationships between period, amplitude, and velocity, including seasonal and depth-dependent trends.
- Produces figures showing temporal variability and propagation characteristics across the OBS region.

**2. OBS_analysis.m**
- Processes ocean-bottom seismic and pressure data from Dongsha Atoll to identify and characterize internal wave signals.
- Applies filtering, artifact removal, and time-windowing to extract amplitude, period, and phase of each detected wave.
- Merges detections with seismic channels to validate internal wave signatures and exports a curated event catalog.
   
**3. himawari_analysis.m**
- Processes 2015–2024 satellite detections of nonlinear internal waves near Dongsha Atoll to derive propagation speed and back azimuth.
- Combines annual datasets with bathymetry to quantify depth-dependent and seasonal trends in wave speed and direction.
- Generates multi-year visualizations, including velocity–depth slopes, monthly variability, and polar plots of propagation direction.

**4. calculate_velocity.m**
- Computes propagation speed and back-azimuth of nonlinear internal waves from sequential satellite SVG detections.
- Tracks wavefront displacement through interpolation and perpendicular distance measurement between consecutive images.
- Outputs per-frame estimates of speed, uncertainty, and direction, with options for spatial constraints and diagnostic visualization.
