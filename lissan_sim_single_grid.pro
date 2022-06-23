FUNCTION lissan_sim_transmission_patterns, collimator_param, xfront, yfront, incidence_angles, pass_grid_1=pass_grid_1, pass_grid_2=pass_grid_2

  ;+
  ; :project:
  ;   LISSAN imaging system simple simulation
  ;
  ; :description:
  ;   This function returns the transmission value through a pair of LISSAN grids, for a point-source at infinite distance at a given incidence angle,
  ;   for a position xfront and yfront on the front grid, as well as the position of arrival of the photon on the detector plane
  ;
  ; :inputs:
  ;   collimator_param: [p1, omega1, p2, omega2, D1, D2, h_fgrid, w_fgrid] were lengths are in mm and angles in rad
  ;   xfront, yfront, position in the transmission patterm
  ;   incidence_angles: two-element array containing the incidence angles in two planes, in degrees
  ;
  ; :output:
  ;   The function returns the value (x,y) of the position of the photon on the detector
  ;
  ; :keywords:
  ;   pass_grid_1: out, optional: 1 if the photon passes the front grid, 0 oterhwise
  ;   pass_grid_2: out, optional: 1 if the photon passes the rear grid, 0 oterhwise
  ;
  ; :called by:
  ;   
  ;
  ; :history:
  ;   2022/06/22 Sophie Musset (ESTEC) - initial release, written from stix_sim_transmission_patterns
  ;
  ; :supporting documentation:
  ;   To follow the calculation you will need to refer to a document called "STIX imaging simulation in IDL"
  ;-

  ;--------------------------------------------------------
  ; read grid parameters
  ;--------------------------------------------------------

  p1 = collimator_param[0]
  omega1 = collimator_param[1]
  p2 = collimator_param[2]
  omega2 = collimator_param[3]
  D1 = collimator_param[4]
  D2 = collimator_param[5]
  h_grid = collimator_param[6]
  w_grid = collimator_param[7]

  ;--------------------------------------------------------
  ; read incidence angles
  ;--------------------------------------------------------

  deg2rad = !pi/180
  phi_x = incidence_angles[0]*deg2rad ; in rad
  phi_y = incidence_angles[1]*deg2rad ; in rad

  ;--------------------------------------------------------
  ; calculate position on front grid
  ;--------------------------------------------------------

  ; coordinate transformation
  rfront = sqrt(xfront^2+yfront^2)
  alphafront = atan(xfront/yfront)
  ; hfront*p1 is the height of the photon arrival location in the grid frame, i.e. in a frame where the slits are horizontal
  ; and with one of the slit aligned with the position xfront=0 and yfront=0
  hfront = rfront*sin(alphafront-omega1)/p1

  ;--------------------------------------------------------
  ; determine if we pass the front grid
  ;--------------------------------------------------------

  res1 = hfront mod 2
  IF res1 LT 0 THEN res1 = res1+2 ; take care of negative numbers
  ;IF res1 LT 1 THEN pass_grid_1=1 ELSE pass_grid_1 = 0
  IF res1 LT 0.5 OR res1 GT 1.5 THEN pass_grid_1 = 0 ELSE pass_grid_1 = 1

  ;--------------------------------------------------------
  ; calculate position on rear grid
  ;--------------------------------------------------------

  ; position of photon given it's angle of incidence and the distance between the two grids
  xrear = xfront + D1*tan(phi_x)
  yrear = yfront + D1*tan(phi_y)

  ;--------------------------------------------------------
  ; determine if we pass the rear grid
  ;--------------------------------------------------------

  ; assuming that the rear grid as the same dimension than the front grid (which is NOT the case but is taking care of later)
  ; then if the photons goes outside of the boundaries of the grid, then it does not pass.
  IF xrear GT w_grid/2. OR xrear LT -w_grid/2. or yrear GT h_grid/2. or yrear LT -h_grid/2. THEN BEGIN
    pass_grid_2 = 0
    ;      print, 'there will be a shadow effect'
  ENDIF ELSE BEGIN
    ; do the same kind of calculation seen for front grid
    arear = sqrt(xrear^2+yrear^2)
    alpharear = atan(xrear/yrear)
    hrear = arear*sin(alpharear-omega2)/p2
    res2 = hrear mod 2
    IF res2 lt 0 then res2 = res2+2 ; take care of negative numbers
    IF res2 LT 0.5 OR res2 GT 1.5 theN pass_grid_2 = 0 else pass_grid_2 = 1
  ENDELSE

  ;--------------------------------------------------------
  ; calculate position on detector
  ;--------------------------------------------------------

  xdet = xfront + (D1+D2)*tan(phi_x)
  ydet = yfront + (D1+D2)*tan(phi_y)

  RETURN, [xdet,ydet]
END



FUNCTION lissan_sim_single_grid, incidence_angles, grid_parameters, npoints=npoints, n_counts=n_counts, loud=loud, plot_interm=plot_interm, $
  noplot=noplot, motif=motif, m_pixels=m_pixels, positions=positions, theorique=theorique

  ;+
  ; :project:
  ;   LISSAN imaging system simple simulation
  ;   
  ; :description:
  ;   This function simulate output of a subcollimator of LISSAN (a pair or grids and a detector) to a point source located at infinite distance, arriving at a certain incidence angle.
  ;   It will produce the moiré pattern in an 2D array of size npoints
  ;   A list of (x,y) position on the detector for n_counts photons (random position)
  ;   The value of the intensiy seen in each of the 8 pixels of the detector (we ignored the 4 small pixels so far)
  ;
  ; :inputs:
  ;   incidence_angles: two-element array containing the incidence angles in two planes, in degrees
  ;   grid_parameters: four-element array containing the pitches p1 and p2, and the orientation omega1 and omega2, of the front and rear grid respectively. pitches are in mm, angles are in degrees
  ;
  ; :output:
  ;   The function returns a 2*4 element array with the values in the 8 pixels of the detector
  ;
  ; :keywords:
  ;   npoints:   double,                in, default=5d2, number of pixels in the pattern array
  ;   n_counts:  double,                in, default=2d4, number of photons to be considered. This is a number of photons that would arrive on the detector surface if there was no grid
  ;                                                                                          The final number of photons that arrive on the detector is lower, probably by a factor 4.
  ;   loud:      int,                   in, default=0, set to 1 to print intermediate messages
  ;   plot_interm, int,                 in, default=0, set to 1 to plot intermediate plots
  ;   noplot, int                       in, default=0, set to 1 to not plot anything
  ;   motif:     npoints*npoints array, out, optional, image of the moire pattern created
  ;   m_pixels                      out, optional, image of the pixels
  ;   positions: 2*n array,             out, optional, positions in (x,y) of each of the n photons arriving on the detector
  ;   theorique, int                    in, default=0, set to 1 to return theorical values of pixels instead of simulated ones from random counts
  ;
  ; :call:
  ;   lissan_sim_transmission_patterns
  ;
  ; :example:
  ;   lissan_sim_single_grid, [0.001, 0.0], [0.320259, 0.330680, 150.52, 149.462], /loud
  ;
  ; :history:
  ;   2022/06/22: Sophie Musset (ESTEC): initial release, modified version of stix_sim_single_grid
  ;
  ; :supporting documentation:
  ;   To follow the caculation you will need to refer to a document called "STIX imaging simulation in IDL"
  ;-

  ;----------------------------------------------------------------------------------------
  ; define default values
  ;----------------------------------------------------------------------------------------

  DEFAULT, npoints, 5d2
  DEFAULT, n_counts, 2d4
  DEFAULT, loud, 1
  DEFAULT, plot_interm, 1
  DEFAULT, noplot, 0
  DEFAULT, theorique, 0

  ;----------------------------------------------------------------------------------------
  ; define constant values
  ;----------------------------------------------------------------------------------------

  deg2rad = !pi/180

  ;----------------------------------------------------------------------------------------
  ; define machanical values in STIX
  ;----------------------------------------------------------------------------------------

  p1 = grid_parameters[0]/2. ; width of the slit of front grid ;**********************
  p2 = grid_parameters[1]/2. ; width of the slit of rear grid  ;**********************
  omega1 = grid_parameters[2]*deg2rad ; orientation of the front grid in rad
  omega2 = grid_parameters[3]*deg2rad ; orientation of the rear grid in rad

  D1 = 3d3 + 30.   ; mm distance between grids + width of grid
  D2 = 20. + 30.    ; mm distance between rear grid and detector + width of grid
  h_det = 100.       ; mm detector height
  w_det = 100.       ; mm detector width
  h_fgrid = 150.     ; mm front grid height
  w_fgrid = 150.     ; mm front grid width
  h_rgrid = 110.     ; mm rear grid height
  w_rgrid = 110.     ; mm rear grid width

  collimator_param = [p1, omega1, p2, omega2, D1, D2, h_fgrid, w_fgrid]

  ;----------------------------------------------------------------------------------------
  ; Create two array representing the front and rear grids, of the size of the front grid
  ;----------------------------------------------------------------------------------------

  xfront = indgen(npoints+1, /double)*w_fgrid/(npoints) - w_fgrid/2.
  yfront = indgen(npoints+1, /double)*h_fgrid/(npoints) - h_fgrid/2.

  frontgrid_transmission = fltarr(n_elements(xfront), n_elements(yfront))
  reargrid_transmission = fltarr(n_elements(xfront), n_elements(yfront))


  print, 'inside proc, incidence angles = ', incidence_angles
  ;----------------------------------------------------------------------------------------
  ; Calculate the moire pattern (tranmission patterns)
  ;----------------------------------------------------------------------------------------

  FOR k = 0, npoints DO BEGIN                                                         ; go through all the x positions
    IF k mod 100 EQ 0 and loud EQ 1 THEN print, 'enter k =', k, '/' ,npoints          ; print intermediate message if needed
    FOR i = 0, npoints DO BEGIN                                                       ; go through the y positions

      pos_ondet = lissan_sim_transmission_patterns(collimator_param, xfront[k], yfront[i], incidence_angles, pass_grid_1=pass_grid_1, pass_grid_2=pass_grid_2)
      frontgrid_transmission[k,i] = pass_grid_1
      reargrid_transmission[k,i] = pass_grid_2

    ENDFOR
  ENDFOR

  moire_pattern = frontgrid_transmission*reargrid_transmission

  x_on_det = where(xfront LT w_det/2. and xfront GT -1*w_det/2., nx )
  y_on_det = where(yfront LT h_det/2. and yfront GT -1*h_det/2., ny )
  xdet = xfront[x_on_det]
  ydet = yfront[y_on_det]
  moire_on_detector = moire_pattern[ x_on_det[0]:x_on_det[nx-1], y_on_det[0]:y_on_det[ny-1] ] ; detector has the size of the detector

  ;----------------------------------------------------------------------------------------
  ; Plot moire pattern
  ;----------------------------------------------------------------------------------------

  IF plot_interm EQ 1 AND noplot NE 1 THEN BEGIN
    imoire = image(moire_pattern, title='Transmission (moire pattern)')
  ENDIF

  IF plot_interm EQ 1 AND noplot NE 1 THEN  idet = image(moire_on_detector, title='Moire on det.')

  ;----------------------------------------------------------------------------------------
  ; Calculate resulting values in detector pixels
  ;----------------------------------------------------------------------------------------

  pixel_val = fltarr(4,4)
  detector_pixelise = moire_on_detector
  FOR k=0,3 DO BEGIN ; k is on the x axis {horizontal}
    FOR i =0,3 DO BEGIN ; i is on the y axis {vertical}
      sel_x = where(xdet GT (i*w_det/4. - w_det/2.) AND xdet LT ((i+1)*w_det/4. - w_det/2.))
      sel_y = where(ydet GT (k*h_det/4. - h_det/2.) AND ydet LT ((k+1)*h_det/4. - h_det/2.))
      pix = moire_on_detector[ min(sel_x):max(sel_x), min(sel_y):max(sel_y)]
      pixel_val[i,k] = total(pix)
      detector_pixelise[i*nx/4.:(i+1)*nx/4.-1, k*ny/4.:(k+1)*ny/4.-1] = pixel_val[i,k]
    ENDFOR
  ENDFOR

  ;----------------------------------------------------------------------------------------
  ; Plot resulting detector pixel values
  ;----------------------------------------------------------------------------------------

  IF plot_interm EQ 1 AND noplot NE 1 THEN idet = image(detector_pixelise, axis_style=2, xmajor=0, ymajor=0, xminor=0, yminor=0, title='theoretical values')

  ;----------------------------------------------------------------------------------------
  ; Generate random counts on detector using the moire pattern
  ;----------------------------------------------------------------------------------------

  listx = list()
  listy = list()

  ; distribute photons position on detector: random location over the detector are
  random = randomu(1,n_counts,2,/uniform)
  randx = random[*,0]*n_elements(xdet)
  randy = random[*,1]*n_elements(ydet)

  FOR k=0, n_counts-1 DO BEGIN
    IF moire_on_detector[randx[k],randy[k]] EQ 1 THEN BEGIN
      listx.add, randx[k]
      listy.add, randy[k]
    ENDIF
  ENDFOR

  arrayx = listx.ToArray()
  arrayy = listy.ToArray()

  nc_final = n_elements(listx)
  open_area = nc_final/n_counts
  IF loud EQ 1 THEN print, 'open area for random counts is ', open_area, ' with ', nc_final, ' counts'

  ;----------------------------------------------------------------------------------------
  ; Plot repartition of randon counts on detector
  ;----------------------------------------------------------------------------------------

  dx=[0,w_det/2]
  dy=[0,h_det/2]
  xx=[dx[0],dx[1],dx[1],dx[0],dx[0]]
  yy=[dy[0],dy[0],dy[1],dy[1],dy[0]]

  IF plot_interm EQ 1 AND noplot NE 1 THEN BEGIN
    s = scatterplot(xdet[arrayx], ydet[arrayy] ,$
      symbol='.', sym_size=5, sym_filled=1, aspect_ratio=1, xr=minmax(xdet), yr=minmax(ydet), axis_style=2, xmajor=0, ymajor=0, xminor=0, yminor=0, title='Random counts')
    p1 = polyline(xx, yy, target=s, /data)
    p2 = polyline(xx-w_det/4, yy, target=s, /data)
    p3 = polyline(xx, -yy, target=s, /data)
    p4 = polyline(xx-w_det/4, -yy, target=s, /data)
    p5 = polyline(xx-w_det/2, yy, target=s, /data)
    p6 = polyline(xx-w_det/2, -yy, target=s, /data)
    p7 = polyline([-1,1]*dx[1], [dy[1],dy[1]]/2, target=s, /data)
    p8 = polyline([-1,1]*dx[1], [dy[1],dy[1]]/(-2), target=s, /data)
  ENDIF

  ;----------------------------------------------------------------------------------------
  ; Calculate resulting values in detector pixels
  ;----------------------------------------------------------------------------------------

  r_pixel_val = fltarr(4,4)
  r_detector_pixelise = moire_on_detector
  FOR k=0,3 DO BEGIN ; k is on the x axis {horizontal}
    FOR i =0,3 DO BEGIN ; i is on the y axis {vertical}
      r_pix = where(xdet[arrayx] GT i*w_det/4. - w_det/2. AND xdet[arrayx] LT (i+1)*w_det/4. - w_det/2. AND ydet[arrayy] GT k*h_det/4. - h_det/2. AND ydet[arrayy] LT (k+1)*h_det/4. - h_det/2.)
      r_pixel_val[i,k] = N_elements(r_pix)
      IF loud EQ 1 THEN print, 'pixel ',i,N_elements(r_pix),' counts/s'
      r_detector_pixelise[i*nx/4.:(i+1)*nx/4.-1, k*ny/4.:(k+1)*ny/4.-1] = r_pixel_val[i,k]
    ENDFOR
  ENDFOR
  IF loud EQ 1 THEN print, total(r_pixel_val), ' counts in pixels'

  ;----------------------------------------------------------------------------------------
  ; Plot resulting detector pixel values
  ;----------------------------------------------------------------------------------------

  IF plot_interm EQ 1 AND noplot NE 1 THEN idet = image(r_detector_pixelise, axis_style=2, xmajor=0, ymajor=0, xminor=0, yminor=0, title= 'random counts')

  ;-----------------------------------------------------------------------------------------
  ; everything in one plot
  ;-----------------------------------------------------------------------------------------

  IF noplot NE 1 THEN BEGIN
    imoire = image(moire_pattern, title='Moire pattern', layout=[5,1,1], dimensions = [2100,500])
    idet = image(moire_on_detector, title='Moire on det.', layout=[5,1,2], /current)
    idet = image(detector_pixelise, axis_style=2, xmajor=0, ymajor=0, xminor=0, yminor=0, title='Theoretical values', layout=[5,1,3], /current)
    s= scatterplot(xdet[arrayx], ydet[arrayy] , symbol='.', sym_size=5, sym_filled=1, aspect_ratio=1, $
      axis_style=2, xmajor=0, ymajor=0, xminor=0, yminor=0, xr=[-1*w_det/2., w_det/2.], yr=[-1*h_det/2., h_det/2.], $
      title='Random counts', layout=[5,1,4], /current)
    p1 = polyline(xx, yy, target=s, /data)
    p2 = polyline(xx-w_det/4, yy, target=s, /data)
    p3 = polyline(xx, -yy, target=s, /data)
    p4 = polyline(xx-w_det/4, -yy, target=s, /data)
    p5 = polyline(xx-w_det/2, yy, target=s, /data)
    p6 = polyline(xx-w_det/2, -yy, target=s, /data)
    i=image(r_detector_pixelise, axis_style=2, xmajor=0, ymajor=0, xminor=0, yminor=0, title= 'Value with random', layout=[5,1,5], /current)
  ENDIF

  ;-----------------------------------------------------------------------------------------
  ; Return
  ;-----------------------------------------------------------------------------------------

  motif = moire_on_detector
  positions = [[xdet[arrayx]], [ydet[arrayy]]]
  m_pixels = r_detector_pixelise
  IF theorique EQ 1 THEN RETURN, pixel_val ELSE RETURN, r_pixel_val

END