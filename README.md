# Decadal Analysis of Nonlinear Internal Waves in the South China Sea Using Seafloor Instrumentation and Satellite

## Research summary
- Ocean-bottom pressure measurements reveal variations in local sea state
- Seafloor instrumentation can measure pressure changes related to nonlinear internal waves and be used to observe wave period and amplitude changes over time, which cannot be obtained from satellite imagery
- Geostationary satellite imagery can monitor decadal-scale variations, or lack thereof, in nonlinear internal wave propagation speeds and arrival directions
- Integrating observation methods reveals nonhydrostatic processes influencing seafloor pressure signals


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
