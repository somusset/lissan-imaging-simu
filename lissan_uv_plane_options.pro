PRO lissan_uv_plane_options, pitches=pitches, angles=angles, u, v
  
  deg2rad = !pi/180
  rad2deg = 180/!pi
  rad2arcsec = 3600d*180d/!pi
  
  ;-------------------------------------
  ; LISSAN constants
  ;-------------------------------------
  
  ; distance between the grids
  D = 3 ; meters
  ncol = 15 ; number of collimators used for imaging
  

  ; ------------------------------------
  ; u,v coverage
  ;-------------------------------------

;  ; option 4 different pitches
;  resolutions = 1d*[8,15,30,60] ; arcsec
;  start_angles = 1d*[0,45,90,135] ; in degrees
;  angle_step = 25d ; deg
  
;  pitch_values = 2*resolutions/3600*deg2rad*D ; pitch = hole+slit in meters here
;  pitches = [replicate(pitch_values[0],4),replicate(pitch_values[1],4),replicate(pitch_values[2],4),replicate(pitch_values[3],4)]
;  angles = [start_angles, start_angles+angle_step, start_angles+2*angle_step, start_angles+3*angle_step]


;  ; option 8 different pitches FAILED STILL WITH 4 PITCHES
;  resolutions = 1d*[8,12,17,23,30,40,55,70] ; arcsec
;  start_angles = 1d*indgen(8)*44 ; in degrees
;  angle_step = 30d ; deg ; option2
;
;  pitch_values = 2*resolutions/3600*deg2rad*D ; pitch = hole+slit in meters here
;  pitches = [replicate(pitch_values[0],4),replicate(pitch_values[1],4),replicate(pitch_values[2],4),replicate(pitch_values[3],4)]
;  angles = [start_angles, start_angles+angle_step, start_angles+2*angle_step, start_angles+3*angle_step]

;  ; option 8 different pitches FAILED STILL WITH 4 PITCHES
;  resolutions = 1d*[8,12,17,23,30,40,55,70] ; arcsec
;  start_angles = 1d*indgen(8)*40 ; in degrees ; option 3 was with 25d
;  angle_step = 90d ; deg ; option2
;
;  pitch_values = 2*resolutions/3600*deg2rad*D ; pitch = hole+slit in meters here
;  pitches = [pitch_values,pitch_values]
;  angles = [start_angles, start_angles+angle_step]

  resolutions = 1d*[8,16,32,64,128] ; arcsec
  start_angles = 1d*indgen(3)*60+15. ; in degrees ; option 3 was with 25d
  angle_step = 12d ; deg ; option2

  pitch_values = 2*resolutions/3600*deg2rad*D ; pitch = hole+slit in meters here
  pitches = [replicate(pitch_values[0],3),replicate(pitch_values[1],3), replicate(pitch_values[2],3), replicate(pitch_values[3],3),replicate(pitch_values[4],3)]
  angles = [start_angles, start_angles+angle_step, start_angles+2*angle_step,  start_angles+3*angle_step,  start_angles+4*angle_step]

 
  ; calcul les composantes u,v de la visibilite
  
  module = 2*D/pitches ; module de la visibilite
  u = module*cos(angles*deg2rad)/rad2arcsec ; u coordinate in arcsec-1
  v = module*sin(angles*deg2rad)/rad2arcsec ; v coordinate in arcsec-1
  
  s=scatterplot([u,-u],[v,-v], aspect_ratio=1, title='LISSAN possible (u,v) plane coverage', axis_style=2, dimensions=[600,600], sym_thick=2, sym_filled=1, sym_size=1.2, xtitle='arcsec-1', ytitle='arcsec-1')
 ; FOR k=0, n_elements(module)-1 DO e=ellipse(0,0,major=module[k]/rad2arcsec,/overplot,/data, fill_transparency=100) 
  
  ; phase calculatation to come, for the moment assume all 1
  phases = replicate(1,ncol)
  stop
  ; -----------------------------------
  ; calculate parameters for each grids
  ; -----------------------------------
  
;  pitch_front = fltarr(16)
;  pitch_rear = fltarr(16)
;  orientation_front = fltarr(16)
;  orientation_rear = fltarr(16)
;  
;  FOR k=0,15 DO BEGIN
;    temp = get_lissan_grid_parameters(pitches[k],angles[k],phases[k])
;    pitch_front[k] = temp[0]
;    pitch_rear[k] = temp[1]
;    orientation_front[k] = temp[2]
;    orientation_rear[k] = temp[3]
;  ENDFOR
  
  ; initialize a variable containing the grid parameters
  grid_param = DBLARR(ncol,4) ; pitch of the front grid, pitch of the rear grid, orientation of the front grid, orientation of the rear grid
  ; for each collimator, calculate the grid parameters by calling calcul_stix_grid_param with the pitch, orientation and phase for that collimator
  FOR k=0,ncol-1 DO BEGIN
    res = get_lissan_grid_parameters(pitches[k], angles[k], phases[k])
    grid_param[k,*] = res
  ENDFOR

  ;----------------------------------------------------------
  ; IF NEEDED calculate the pixel values for each  collimator
  ;----------------------------------------------------------
   angles_sec = [0.,0.] ; define incidence angles in arcsec
   photon_flux = 1d6
   angles = angles_sec / 3600d ; in degree

;  ; make two calculatations just to produce figures to show examples
;  k=14
;  angles_sec = [100.,0.] ; define incidence angles in arcsec
;  angles = angles_sec / 3600d ; in degree
;  r = lissan_sim_single_grid(angles, reform(grid_param[k,*]), loud=1, motif=motif, n_counts=photon_flux, plot_interm=0)
;  stop
;  angles_sec = [0.,0.] ; define incidence angles in arcsec
;  angles = angles_sec / 3600d ; in degree
;  r = lissan_sim_single_grid(angles, reform(grid_param[k,*]), loud=1, motif=motif, n_counts=photon_flux, plot_interm=0)
;  stop
  
  
;  IF calculate_moire EQ 1 THEN BEGIN
    pixels = dblarr(ncol,4)
    nameangles = strtrim(round(angles[0]),2)+'_'+strtrim(round(angles[1]),2)
    row=0; here we choose the first row of pixels
    FOR k=0, ncol-1 DO BEGIN
      r = lissan_sim_single_grid(angles, reform(grid_param[k,*]), loud=1, motif=motif, n_counts=photon_flux, plot_interm=0, /noplot)
      pixels[k,*] = reform(r[*,row]) ; here we choose the first row of pixels
     ; print, 'flux ', total(r), ' counts/s'
     ; stop
    ENDFOR
    ;save it because it takes time
    save, pixels, filename='save_pixels_'+nameangles+'.sav'
    print, 'save_pixels_'+nameangles+'.sav', ' has been saved'
;  ENDIF ELSE BEGIN
;    restore, moire_file
;  ENDELSE
  
  ;-------------------------------------
  ; LISSAN constants
  ;-------------------------------------
  l = 25 ; largeur pixel in mm
  h = 25 ; longueur pixel in mm
  eff_area = l*h*100 ; effective area in cm2
  
  ;-----------------------------------------------------
  ; calculate count flux
  ;-----------------------------------------------------
  flux = dblarr(ncol) ; will be in counts/s/cm2
  FOR k=0,ncol-1 DO flux[k] = total(reform(pixels[k,*])) /eff_area
  print, flux, 'counts/s/cm2'

  ;-----------------------------------------------------
  ; calculate the corresponding visibilities
  ;-----------------------------------------------------
  
  visibilities = lissan_sim_calculate_visibilities(pixels)

  ;-----------------------------------------------------
  ; create dirty map and plot it
  ;-----------------------------------------------------

  weights = 1/pitches ; default weights
  DEFAULT, npix, 80
  DEFAULT, pixsize, 0.5
  DEFAULT, xc, 0
  DEFAULT, yc, 0 
  facto=2*!pi
  ;-------------------------------------------------------------
  ; calculate the x and y values of the image
  ;-------------------------------------------------------------
  xvalues = indgen(npix)*pixsize - (npix*pixsize/2 - xc)
  yvalues = indgen(npix)*pixsize - (npix*pixsize/2 - yc)
  ;-------------------------------------------------------------
  ; calculate the amplitude and phase of each visibitily
  ;-------------------------------------------------------------
  real = real_part(visibilities)
  imgn = imaginary(visibilities)
  amplitude = sqrt( real^2 + imgn^2 )
  phase = atan(visibilities, /phase)
  
  ; initialisation
  dmap = dblarr(npix,npix)
  ; calculation of the value in each pixel
  FOR i=0, npix-1 DO BEGIN
    FOR j=0, npix-1 DO BEGIN
      dmap[i,j] = total(weights*amplitude*cos(facto*(u*xvalues[i]+v*yvalues[j])+phase))/total(weights)
    ENDFOR
  ENDFOR
  
  xmap = indgen(npix)*pixsize + xc - npix/2*pixsize
  ymap = indgen(npix)*pixsize + yc - npix/2*pixsize
  i=image(dmap, rgb_table=5, margin=[0.05,0.05,0.05,0.05])
  i=image(dmap, xmap,ymap, rgb_table=5, margin=[0.2,0.2,0.2,0.2], AXIS_STYLE=1, xtickdir=1, ytickdir=1, xtitle='arcsec', ytitle='arcsec')
  
  ;-------------------------------------------------------------
  ; do intermediate maps (for info plot)
  ;-------------------------------------------------------------

  dmap1 = dblarr(npix,npix)
  FOR i=0, npix-1 DO BEGIN
    FOR j=0, npix-1 DO BEGIN
      dmap1[i,j] = total(weights[0]*amplitude[0]*cos(facto*(u[0]*xvalues[i]+v[0]*yvalues[j])+phase[0]))/total(weights[0])
    ENDFOR
  ENDFOR

  dmap2 = dblarr(npix,npix)
  FOR i=0, npix-1 DO BEGIN
    FOR j=0, npix-1 DO BEGIN
      dmap2[i,j] = total(weights[0:1]*amplitude[0:1]*cos(facto*(u[0:1]*xvalues[i]+v[0:1]*yvalues[j])+phase[0:1]))/total(weights[0:1])
    ENDFOR
  ENDFOR

  dmap4 = dblarr(npix,npix)
  FOR i=0, npix-1 DO BEGIN
    FOR j=0, npix-1 DO BEGIN
      dmap4[i,j] = total(weights[0:3]*amplitude[0:3]*cos(facto*(u[0:3]*xvalues[i]+v[0:3]*yvalues[j])+phase[0:3]))/total(weights[0:3])
    ENDFOR
  ENDFOR

  dmap8 = dblarr(npix,npix)
  FOR i=0, npix-1 DO BEGIN
    FOR j=0, npix-1 DO BEGIN
      dmap8[i,j] = total(weights[0:7]*amplitude[0:7]*cos(facto*(u[0:7]*xvalues[i]+v[0:7]*yvalues[j])+phase[0:7]))/total(weights[0:7])
    ENDFOR
  ENDFOR


  stop
END