FUNCTION get_lissan_grid_parameters, pitch, orientation, phase, moire_width_mm=moire_width_mm
  
  
  ;+
  ; :project:
  ;   LISSAN imaging system simple simulation
  ;
  ; :description:
  ;   This function takes the average pitch and orientation (and phase) for LISSAN subcollimators
  ;   and return the exact pitch and orientation of the rear and front grids needed to have a moire pattern
  ;   sampled on a LISSAN array of detectors
  ;
  ; :inputs:
  ;   pitch: mean pitch used for the visibility in meters
  ;   orientation: mean orientation used for the visibility in degree
  ;   phase: phase of the grid (1 or -1)
  ;   
  ; :output:
  ;   The function returns a 4 element array with the values of: pitch of the front grid in mm, pitch of the rear grid in mm, orientation of the front grid in degrees, orientation of the rear grid in degrees
  ;
  ; :keywords:
  ;   moire_width_mm: default is 100. Width of the moire pattern on detector plane, in mm
  ;   
  ; :history:
  ;   2022/06/22: Sophie Musset (ESTEC): initial release, modified version of get_stix_grid_parameters
  ;-
 
  DEFAULT, moire_width_mm, 100.

  ; convert pitch in mm
  pitch = pitch*1d3

  du = 1./moire_width_mm/2.

  deg2rad = !pi/180
  rad2deg = 180/!pi

  ; calcul les composantes u,v de la visibilite
  module = 1./pitch ; module de la visibilite
  u = module*cos(orientation*deg2rad)
  v = module*sin(orientation*deg2rad)

  ; calcul des composantes u pour la front and rear grids
  uf = u - phase*du
  ur = u + phase*du

  modulef = sqrt( uf^2 + v^2)
  moduler = sqrt( ur^2 + v^2)

  pitchf = 1/modulef
  pitchr = 1/moduler

  orientationf = atan(v/uf)*rad2deg
  orientationr = atan(v/ur)*rad2deg
  IF orientationf LT 0 THEN orientationf=180.+orientationf
  IF orientationr LT 0 THEN orientationr=180.+orientationr

  RETURN, [pitchf, pitchr, orientationf, orientationr] ; pitches are in mm here!!!

END