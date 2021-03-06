FUNCTION combined_mag, m1, m2
;-------------------------------------------------------------------------------------------------------
;  Short IDL function that takes as input the apparent magnitude of each of two stars, 
;    and computes and returns the apparent magnitude of the two stars combined.
;
;  INPUTS: m1, the magnitude of the first star
;          m2, the magnitude of the second star
; 
;  OUTPUT: returns the magnitude of the two stars combined if they were seen
;             as a single object. 
; 
;  This is a stand-alone function and can be used without compiling other programs. 
;
;  (This function has been tested against by-hand calculations and shown to work correctly.)
;-------------------------------------------------------------------------------------------------------
mag1 = float(m1)         ; ensures that the 1st input magnitude is a float for calculations
mag2 = float(m2)         ; ensures that the 2nd input magnitude is a float for calculations

flux1 = 10.^(-0.4*mag1)  ; find the flux from the 1st star (F0 can be taken to be 1; it cancels out) 
flux2 = 10.^(-0.4*mag2)  ; find the flux from the 2nd star (F0 can be taken to be 1; it cancels out) 

fluxtotal = flux1 + flux2    ; calculate the total flux from both stars (flux is additive)

mag_combo = ( alog10(fluxtotal) ) / (-0.4)  ; convert the total flux to a magnitude (F0 cancels here)
return, mag_combo                           ; return the magnitude of the two stars combined

END   ; combined_mag function