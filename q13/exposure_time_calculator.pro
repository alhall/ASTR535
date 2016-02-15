PRO exposure_time_calculator, diameter, wavelengths, magnitudes, background, skyarea, pixelscale, time 
;------------------------------------------------------------------------------------------------------------------
; Agnar Hall, ASTR 535
;
; this will eventually calculate exposure times
;
; INPUTS: diameter - telescope diameter in meters
;         wavelengths - array of wavelengths included in the desired bandpass
;         magnitudes - array corresponding to the magnitude of the target at each wavelength
;         background - array corresponding to the magnitude per square arcsec of the background at each wavelength
;         skyarea - area of the aperture in square arcsec
;         pixelscale - plate/pixel scale of the CCD being used in square arcsec per pixel
;         time - exposure time in seconds
; 
; OUTPUT: The signal to noise for the given exposure time in the specified bandpass 
;------------------------------------------------------------------------------------------------------------------
  FORWARD_FUNCTION telescope_tp
  FORWARD_FUNCTION filter_tp
  FORWARD_FUNCTION instrument_tp
  FORWARD_FUNCTION detector_efficiency
  FORWARD_FUNCTION telescope_area
  FORWARD_FUNCTION flux_from_mag
  FORWARD_FUNCTION atmo_transmission
  FORWARD_FUNCTION system_transmission
  FORWARD_FUNCTION count_equation
  FORWARD_FUNCTION readout_noise

  S_over_t = count_equation(wavelengths, magnitudes, diameter)    ; the signal per time S/t
  B_over_t = count_equation(wavelengths, background, diameter)    ; the background per square arcsec
  BA_over_t = B_over_t * skyarea                                  ; the background B/t
  rn= readout_noise(pixelscale, skyarea)                          ; readout noise
  ; now get signal to noise
  StoN = S_over_t*time / sqrt( (S_over_t*time) + (BA_over_t*time) + rn )
  
END 






 

FUNCTION count_equation, wavelengths, magnitudes, diameter
;------------------------------------------------------------------------------
; given all the terms in the count equation except S/t, 
;   calculates S/t
; 
; takes as input an array of wavelengths in Angstroms, an equally-large 
;   array of magnitudes corresponding to each wavelength, 
;   and a telescope diameter in meters
;------------------------------------------------------------------------------
  FORWARD_FUNCTION telescope_tp
  FORWARD_FUNCTION filter_tp
  FORWARD_FUNCTION instrument_tp
  FORWARD_FUNCTION detector_efficiency
  FORWARD_FUNCTION telescope_area
  FORWARD_FUNCTION flux_from_mag
  FORWARD_FUNCTION atmo_transmission
  FORWARD_FUNCTION system_transmission

  hc = 6.6260755d-27*2.99792458d10  ; erg cm

  n = n_elements(wavelengths) - 1  ; index of last element of wavelength array
  dlambda = (wavelengths[n]-wavelengths[0])/float(n)   ; interval between wavelengths in array
  Area = telescope_area(diameter)
  q = system_transmission(wavelengths)
  a = atmo_transmission(wavelengths)
  F = flux_from_mag(wavelengths, magnitudes)
  
  ; CHANGE THIS TO BE A SUM FOR THE INTEGRAL
  ; also, note a 1d-8: this accounts for converting Angstroms to cm
  ;
  ; count equation: 
  integrand = q*a*F*wavelengths
  S_over_t = Area/hc*total(integrand)*dlambda*1d-8
  return, S_over_t
END





 
FUNCTION telescope_area, diameter
;------------------------------------------------------------------------------
; accepts a telescope diameter in meters and 
;  returns the area of that telescope's primary  
;  mirror in square centimeters 
;------------------------------------------------------------------------------
  ; first, convert meters to centimeters
  d_cm = float(diameter)*100.
  ; get radius, use radius to get area
  r = d_cm/2.
  A = !PI*r^2.
  ; output the area
  return, A
END





 
FUNCTION flux_from_mag, wavelengths, magnitudes  
;------------------------------------------------------------------------------
; accepts two arrays of equal length:
;  one is an array of wavelengths,
;  the other is the magnitude of an object at each of those wavelengths
;
; uses these wavelengths to get Vega's flux at the same wavelengths, 
;  then uses that Vega flux and the input magnitudes to 
;  determine the flux of the object at all specified wavelengths
;  and returns this flux (an array) 
;------------------------------------------------------------------------------
  ; Vega's F as function of wavelength in optical is roughly:
  F_Vega = (3.6d-9)*(wavelengths/5500.)^2.   ; erg/cm^2/s/A

  ; now get the object's flux
  flux = F_Vega * 10.^(-0.4*magnitudes)
  return, flux 
END





 
FUNCTION readout_noise, pixelscale, skyarea
;------------------------------------------------------------------------------
; accepts a plate/pixel scale and an aperture area and 
;  returns the readout noise in counts
;------------------------------------------------------------------------------
  Npixels = float(skyarea) / float(pixelscale)   ; total number of pixels 
  rn = 5.*Npixels                  ; just a number for now, 5 counts per pixel
  return, rn
END




 
FUNCTION atmo_transmission, wavelengths
;------------------------------------------------------------------------------
; accepts an array of wavelengths and returns the
;  atmospheric transmission for each wavelength as an array
;------------------------------------------------------------------------------
  atmo_trans = 0.8+dblarr(n_elements(wavelengths))    ; just to get #1 of q10
  return, atmo_trans
END






FUNCTION system_transmission, wavelengths
;------------------------------------------------------------------------------
; accepts an array of wavelengths and returns a
;  system transmission for each wavelength as an array by 
;  multiplying the detector efficiency and the telescope, 
;  filter, and instruments throughputs for those wavelengths
; 
; uses four other functions: see the "forward function" block below
;  (these four other functions give the throughputs and detector 
;   efficiency mentioned above) 
;------------------------------------------------------------------------------
  FORWARD_FUNCTION telescope_tp
  FORWARD_FUNCTION filter_tp
  FORWARD_FUNCTION instrument_tp
  FORWARD_FUNCTION detector_efficiency

  ttp = telescope_tp(wavelengths)
  ftp = filter_tp(wavelengths)
  itp = instrument_tp(wavelengths)
  deff = detector_efficiency(wavelengths) 

  sys_trans = ttp*ftp*itp*deff    ; taking advantage of IDL's array arithmetic
  return, sys_trans
END






FUNCTION telescope_tp, wavelengths
;------------------------------------------------------------------------------
; accepts an array of wavelengths and returns a 
;  telescope throughput for each wavelength as an array
;------------------------------------------------------------------------------
  tel_tp = 0.5 + dblarr(n_elements(wavelengths))     ; just to get #1 of q10
  return, tel_tp
END



FUNCTION filter_tp, wavelengths
;------------------------------------------------------------------------------
; accepts an array of wavelengths and returns a 
;  filter throughput for each wavelength as an array
;------------------------------------------------------------------------------
  filt_tp = 1.+ dblarr(n_elements(wavelengths))
  return, filt_tp
END



FUNCTION instrument_tp, wavelengths
;------------------------------------------------------------------------------
; accepts an array of wavelengths and returns an 
;  instrument throughput for each wavelength as an array
;------------------------------------------------------------------------------
  ins_tp = 1.+ dblarr(n_elements(wavelengths))
  return, ins_tp
END



FUNCTION detector_efficiency, wavelengths
;------------------------------------------------------------------------------
; accepts an array of wavelengths and returns a 
;  detector efficiency for each wavelength as an array
;------------------------------------------------------------------------------
  det_eff = 1.+ dblarr(n_elements(wavelengths)) 
  return, det_eff
END