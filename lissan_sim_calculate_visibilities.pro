FUNCTION lissan_sim_calculate_visibilities, pixel_val, phase=phase

  ;+
  ; :project:
  ;   LISSAN imaging system simple simulation
  ;
  ; :description:
  ;   This function takes the values in four STIX pixels in 30 STIX collimators and return an array of 30 complex visibilities
  ;
  ; :inputs:
  ;   pixel_val: 30*4 array of pixel values
  ;
  ; :keyword:
  ;   phase: default is 0, set to 1 to take into account the phase of the collimator
  ;
  ;   
  ; :history:
  ;   2022/06/22, Sophie Musset (ESTEC), initial release, adapted from stix_sim_calculate_visibilities
  ;-

  ;-----------------------------------------------------------------------------------------------
  ; initialisation
  ;-----------------------------------------------------------------------------------------------

  l = 2.5 ; largeur pixel in cm
  h = 2.5 ; longueur pixel in cm
  m1 = 4d/(!pi)^3*l*h*sin(!pi/4.)
  ncol = 15

  vis_real_proxy = dblarr(16)
  vis_imgn_proxy = dblarr(16)

  ;-----------------------------------------------------------------------------------------------
  ; this calculates the values of C-A and D-B, where A, B, C, D are the values of the four pixels.
  ;-----------------------------------------------------------------------------------------------

  FOR k=0, ncol-1 DO BEGIN
    vis_real_proxy[k] = pixel_val[k,2]-pixel_val[k,0] ; C-A
    vis_imgn_proxy[k] = pixel_val[k,3]-pixel_val[k,1] ; D-B
  ENDFOR

  vis_complex = complex(vis_real_proxy, vis_imgn_proxy)/(4*M1)*complex(cos(!pi/4.), sin(!pi/4.)) ; visibility in counts/s/cm2

  RETURN, vis_complex
END
