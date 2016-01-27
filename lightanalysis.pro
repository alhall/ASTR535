;-------------------------------------------------------------------------------------------------------
;
; lightanalysis.pro
;
; Uses the lightconvert function to determine the frequency and energy of a 5500 Angstrom photon and
;  uses the fluxconvert function to determine the conversion factor between Flambda and Fnu at 3000A, 5500A, and 8000A.
;  Also tests the magconvert function on 5500A to convert between STMAG and ABNU magnitude systems.
;
; INPUTS: none (the functions used by this code can be used independently with user input, though)
;
; OUTPUTS: returns the frequency and energy of a 5500A photon
;          returns the conversion factor(s) between Flambda and Fnu at 3000A, 5500A, and 8000A.
;          also returns the difference between m(ABNU) and m(STMAG) for a 5500A photon. 
;
;-------------------------------------------------------------------------------------------------------
PRO lightanalysis

  FORWARD_FUNCTION lightconvert
  FORWARD_FUNCTION fluxconvert
  FORWARD_FUNCTION magconvert

  ; use lightconvert to determine the frequency and energy of a 5500A photon 
  lightinfo = lightconvert(wavelength=5500.)
  print, "For a 5500-Angstom photon, "
  print, "lambda = 5500 Angstroms"
  print, "nu = ", lightinfo.frequency, " Hz"
  print, "E = ", lightinfo.energy, " ergs"
  print, " " 

 ;Use fluxconvert to determine the conversion factor between Flambda and Fnu at 3000A, 5500A, and 8000A.
 conversion3000A = fluxconvert(3000.)
 conversion5500A = fluxconvert(5500.)
 conversion8000A = fluxconvert(8000.)
 print, "At 3000 Angstroms, Flambda [erg/cm2/s/A]= ", conversion3000A.nu_to_lambda, " * Fnu [erg/cm2/s/Hz]."
 print, "  (Therefore, at 3000 Angstroms, Fnu [erg/cm2/s/Hz] = ", conversion3000A.lambda_to_nu, " * Flambda [erg/cm2/s/A].)"
 print, "At 5500 Angstroms, Flambda [erg/cm2/s/A]= ", conversion5500A.nu_to_lambda, " * Fnu [erg/cm2/s/Hz]."
 print, "  (Therefore, at 5500 Angstroms, Fnu [erg/cm2/s/Hz] = ", conversion5500A.lambda_to_nu, " * Flambda [erg/cm2/s/A].)"
 print, "At 8000 Angstroms, Flambda [erg/cm2/s/A]= ", conversion8000A.nu_to_lambda, " * Fnu [erg/cm2/s/Hz]."
 print, "  (Therefore, at 8000 Angstroms, Fnu [erg/cm2/s/Hz] = ", conversion8000A.lambda_to_nu, " * Flambda [erg/cm2/s/A].)"
 print, " "

 ; test out the magconvert function using 5500 Angstroms
 magconversion = magconvert(5500.)
 print, "At 5500 Angstroms, m(STMAG) = m(ABNU) + ", magconversion
 print, "Therefore, at 5500 Angstroms, m(ABNU) = m(STMAG) - ", magconversion

END 





;-------------------------------------------------------------------------------------------------------
;
; lightconvert function
;
; Converts between wavelength, frequency, and energy. Used in the lightanalysis procedure.  
;
; Returns the wavelength, frequency, and energy of a photon by converting one of those properties 
;   (pre-specified by the user) into the other two. 
; Which two conversions are performed is determined by one of three possible keyword choices.
;   Exactly one of these keywords must be specified.
;
; INPUTS: wavelength: the wavelength of the photon in Angstroms
;         frequency: the frequency of the photon in Hertz
;         energy: the energy of the photon in ergs
;   Only one parameter out of these three needs to be specified.
;
; OUTPUTS: creates a structure containing the wavelength in Angstroms, frequency in Hertz, AND
;            energy in ergs of the photon. 
;
;-------------------------------------------------------------------------------------------------------
FUNCTION lightconvert, wavelength=wavelength, frequency=frequency, energy=energy
  c = double(2.99792458*10.^10.)
  h = double(6.62607554*10.^(-27.))

  ; first, make sure the user only specified one of the three parameters

  IF (( n_elements(wavelength) NE 0 ) AND ( n_elements(frequency) NE 0 )) THEN BEGIN         ; if both wavelength and frequency are specified
    print, "Error: too many constraints. Specify only one of the following: wavelength, frequency, energy."

  ENDIF ELSE IF (( n_elements(wavelength) NE 0 ) AND ( n_elements(energy) NE 0 )) THEN BEGIN ; if both wavelength and energy are specified
    print, "Error: too many constraints. Specify only one of the following: wavelength, frequency, energy."

  ENDIF ELSE IF (( n_elements(frequency) NE 0 ) AND ( n_elements(energy) NE 0 )) THEN BEGIN  ; if both frequency and energy are specified
    print, "Error: too many constraints. Specify only one of the following: wavelength, frequency, energy."

  ; now, on to the actual conversions 

  ENDIF ELSE IF ( n_elements(wavelength) NE 0 ) THEN BEGIN    ; if wavelength is specified
    lambda = double(wavelength)
    nu = (1d8)*c/lambda                                       ; use c = lambda*nu and convert Angstroms to cm for calculations in cgs
    E = h*nu                                                  ; use E = h*nu 
    output = {lightproperties, wavelength:lambda, frequency:nu, energy:E}  ; make a structure containing all three properties
    return, output

  ENDIF ELSE IF ( n_elements(frequency) NE 0 ) THEN BEGIN     ; if frequency is specified
    nu = double(frequency)
    lambda = (1d8)*c/nu                                       ; use c = lambda*nu and convert from cgs to Angstroms
    E = h*nu                                                  ; use E = h*nu 
    output = {lightproperties, wavelength:lambda, frequency:nu, energy:E}  ; make a structure containing all three properties
    return, output

  ENDIF ELSE IF ( n_elements(energy) NE 0 ) THEN BEGIN        ; if energy is specified
    E = double(energy)
    nu = E/h                                                  ; use E = h*nu 
    lambda = (1d8)*c/nu                                       ; use c = lambda*nu and convert from cgs to Angstroms
    output = {lightproperties, wavelength:lambda, frequency:nu, energy:E}  ; make a structure containing all three properties
    return, output

  ENDIF      ; wavelength, frequency, or energy has been specified

END          ; lightconvert function





;-------------------------------------------------------------------------------------------------------
;
; fluxconvert function
;
; Calculates the conversion factors between Flambda (in ergs/cm2/s/A) and Fnu (in ergs/cm2/s/Hz).  
;   Used in the lightanalysis procedure.  
;
; Given a wavelength in Angstroms, returns the conversion factors between fluxes Flambda and Fnu
;   (in both directions) in a structure. 
;
; INPUTS: wavelength: the wavelength of interest in Angstroms
;
; OUTPUTS: creates a structure containing the conversion factor from Flambda to Fnu and the
;            conversion factor from Fnu to Flambda.
;-------------------------------------------------------------------------------------------------------
FUNCTION fluxconvert, wavelength 
  c = double(2.99792458*10.^10.)
  lambda = double(wavelength) 
  ; Fnu = Flambda * c / (nu^2) = Flambda * (lambda^2) / c
  ;  so Flambda = Fnu * c / (lambda^2)
  ; NOTE: The strange units of Flambda, erg/cm2/s/A, mean that in the conversion, one of the lambdas in lambda^2 must
  ;   be in A and the other must be in cm (since c is in cm/s). This is accounted for by including a factor of 10^8
  ;   along with lambda^2 : 1 cm = 10^8 A, so 1 A = 10^-8 cm, so lambda(cm) = lambda(A) / 10^8
  lambda_to_nu = c / (lambda^2. / 1d8)
  nu_to_lambda = (lambda^2. / 1d8) / c 
  output = {FluxConversionFactors, lambda_to_nu:lambda_to_nu, nu_to_lambda:nu_to_lambda}
  return, output
END      ; fluxconvert function





;-------------------------------------------------------------------------------------------------------
;
; magconvert function
;
; Calculates the conversion factors between STMAG and ABNU magnitude systems.
;  (Tested in the lightanalysis procedure.) 
;
; Given a wavelength in Angstroms, calculates the additional factor to add/subtract to 
;  m(ABNU)/m(STMAG) to get m(STMAG)/m(ABNU). 
;
; INPUTS: wavelength: the wavelength of interest in Angstroms 
;
; OUTPUTS: calculates the factor to add to m(ABNU) to get m(STMAG) (which is the same as the 
;            factor to subtract from m(STMAG) to get m(ABNU) of course)
;-------------------------------------------------------------------------------------------------------
FUNCTION magconvert, wavelength 
  c = double(2.99792458*10.^10.)
  lambda = double(wavelength)
  conversionfactor = 2.5*alog10(lambda^2. / 1d8 / c) + 27.5   ; note the factor of 1d8 from the fluxconvert funtion is here too
  ; for STMAG to ABNU, the conversion factor is positive; for ABNU to STMAG, the conversion factor is negative
  return, conversionfactor
END 
