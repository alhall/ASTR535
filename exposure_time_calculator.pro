PRO exposure_time_calculator
;------------------------------------------------------------------------------------------------------------------
; this will eventually calculate exposure times, but for now, it'll just
;   give an S/t from the count equation
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

  diameter = 1.  ; 1m telescope
  wavelengths = [5000., 6000.]
  magnitudes = [22., 22.] 

  S_over_t = count_equation(wavelengths, magnitudes, diameter) 

  print, S_over_t, ' photons per second' 

  ; for these numbers (given in q10), the program gets: 
  ; 4.9629263131563395  photons per second
  ;   (for both wavelengths in the array)
  ; Success! 
  
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

  Area = telescope_area(diameter)
  q = system_transmission(wavelengths)
  a = atmo_transmission(wavelengths)
  F = flux_from_mag(wavelengths, magnitudes)
  n = n_elements(wavelengths) - 1  ; index of last element of wavelength array
  
  ; every part of the count equation has been integrated over except
  ;   lambda, so they don't have to be included, and the
  ;   lambda integral just gives lambda squared
  ; also, note a 1d-8: this accounts for converting Angstroms to cm
  ;
  ; count equation: 
  S_over_t = (Area*F*q*a/hc/2.)*( (wavelengths[n])^2. - (wavelengths[0])^2. )*1d-8
  return, S_over_t
END





 
FUNCTION telescope_area, diameter
;------------------------------------------------------------------------------
; accepts a telescope diameter in meters and 
;  returns the area of that telescope's primary  
;  mirror in square centimeters 
;------------------------------------------------------------------------------
  ; first, account for user possibly entering an integer
  d = float(diameter)
  ; now convert meters to centimeters
  d_cm = d*100.
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
  ; (eventually, will have Vega's F as function of wavelength,
  ;  but for now we'll assume it's constant at all wavelengths)
  F_Vega = (3.6d-9)+dblarr(n_elements(wavelengths))

  ; now get the object's flux
  flux = F_Vega * 10.^(-0.4*magnitudes)
  return, flux 
END





 
FUNCTION atmo_transmission, wavelengths
;------------------------------------------------------------------------------
; accepts and array of wavelengths and returns the
;  atmospheric transmission for each wavelength as an array
;------------------------------------------------------------------------------
  atmo_trans = 0.8+dblarr(n_elements(wavelengths))    ; just to get #1 of q10
  return, atmo_trans
END






FUNCTION system_transmission, wavelengths
;------------------------------------------------------------------------------
; accepts and array of wavelengths and returns a
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