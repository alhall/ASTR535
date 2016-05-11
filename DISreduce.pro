PRO DISreduce, datared, datablue, Hered, Heblue, Nered, Neblue, Arred, Arblue, flatlistred, flatlistblue, biaslistred, biaslistblue 
; Agnar Hall
; ASTR535 
;
; INPUTS: datared - the filename of a red DIS spectrum (fits file) 
;         datablue - filename of a blue DIS spectrum (fits file) 
;         Hered & Nered & Arred - HeNeAr spectra containing lines of known wavelength in red 
;         Heblue & Neblue & Arblue - HeNeAr spectra containing lines of known wavelength in blue 
;         flatlistred - a text file containing names of all the fits files of flats (quartz lamps) to use 
;                            for reduction for red, separated by commas 
;         flatlistblue - a text file containing names of all the fits files of flats (quartz lamps) to use for  
;                            reduction for blue, separated by commas  
;         biaslistred - a text file containing names of all the fits files of biases to use for reduction 
;                            for red, separated by commas 
;         biaslistblue - a text file containing names of all the fits files of biases to use for reduction 
;                            for blue, separated by commas 
;
; OUTPUT: reduces DIS spectra (both red and blue) and calibrates wavelengths for each using linear interpolation 
;             of the reference/calibration (HeNeAr) spectra
;
 FORWARD_FUNCTION CCDread
 FORWARD_FUNCTION overscansub
 FORWARD_FUNCTION biassub
 FORWARD_FUNCTION normalize
 FORWARD_FUNCTION flatreduce
 FORWARD_FUNCTION imcombine
 ; read in text files containing list of filenames of flats and biases in red
 rflats = read_csv(flatlistred)
 rflats = flats[0]             ; this converts the 1-field structure into a 1D array of filenames
 rbiases = read_csv(biaslistred)
 rbiases = biases[0]           ; this converts the 1-field structure into a 1D array of filenames
 ; prepare red flats
 rflatlist = flatreduce(rflats, rbiases)
 rmegaflat = imcombine(rflatlist) 
 ; read in text files containing list of filenames of flats and biases in blue
 bflats = read_csv(flatlistblue)
 bflats = flats[0]             ; this converts the 1-field structure into a 1D array of filenames
 bbiases = read_csv(biaslistblue)
 bbiases = biases[0]           ; this converts the 1-field structure into a 1D array of filenames
 ; prepare blue flats
 bflatlist = flatreduce(bflats, bbiases)
 bmegaflat = imcombine(bflatlist) 
 ; apply bias and flat analysis to data spectra
 rdataspec = CCDread(datared)                ; red 
 rdataspec = overscansub(rdataspec)
 rdataspec = biassub(rdataspec, rbiases) 
 rdataspec = rdataspec/rmegaflat
 bdataspec = CCDread(datablue)               ; blue 
 bdataspec = overscansub(bdataspec)
 bdataspec = biassub(bdataspec, bbiases) 
 bdataspec = bdataspec/bmegaflat
 ; now do the same for the reference spectra 
    ; red 
 rHespec = CCDread(Hered)        
 rNespec = CCDread(Nered)                
 rArspec = CCDread(Arred)               
 rrefspec = rHespec + rNespec + rArspec  ; combine HeNeArs into one reference spectrum
 rrefspec = overscansub(rrefspec)
 rrefspec = biassub(rrefspec, rbiases) 
 rrefspec = rrefspec/rmegaflat
    ; blue
 bHespec = CCDread(Heblue)              
 bNespec = CCDread(Neblue)                
 bArspec = CCDread(Arblue)    
 brefspec = bHespec + bNespec + bArspec  ; combine HeNeArs into one reference spectrum           
 brefspec = overscansub(brefspec)
 brefspec = biassub(brefspec, bbiases) 
 brefspec = brefspec/bmegaflat
 ; wavelength calibration 
 imagex = 2098.    ; number of pixels in x direction in fits spectra
 imagey = 1078.    ; number of pixels in y direction in fits spectra 
 pixel = indgen(imagex)   ; number of each pixel in the x direction
 pixelf = findgen(imagex) ; same as above, but floats instead of integers 
 rrefmiddle = rrefspec[*,imagey/2.]   ; this gives a strip in the middle of the spectrum (red) 
 brefmiddle = brefspec[*,imagey/2.]   ; this gives a strip in the middle of the spectrum (blue) 
 thresh = 4.5   ; anything above this (times the median level) counts as a spectral line
 ; now, any line is defined by 1) being high enough and 2) being the peak of a line 
 rpeaks = WHERE( (rrefmiddle GT thresh*median(rrefmiddle)) AND (rrefmiddle[pixel] GT   $
     rrefmiddle[pixel-1]) AND (rrefmiddle[pixel] GT rrefmiddle[pixel+1]) )
 bpeaks = WHERE( (brefmiddle GT thresh*median(brefmiddle)) AND (brefmiddle[pixel] GT   $
     brefmiddle[pixel-1]) AND (brefmiddle[pixel] GT brefmiddle[pixel+1]) )
 ; lines I want to use for calibration: 
    ; red 
 rlinesn = [3, 6, 12, 16, 24, 28]
 rlinesreal = [8424.65, 8115.31, 7635.11, 7245.17, 6402.25, 5875.62]
    ; blue 
 blinesn = [0, 1, 2]
 blinesreal = [3888.65, 4471.48, 5015.68]
 ; select only those lines out of the whole image
 rline_indices = ( pixelf[rpeaks] )[rlinesn]
 bline_indices = ( pixelf[bpeaks] )[blinesn] 
 ; fit wavelengths to pixels with a linear interpolation (y = mx + b) 
 rinterpolation = linfit(rline_indices, rlinesreal)
 rwavelengths = rinterpolation[1]*pixelf + rinterpolation[0]
 binterpolation = linfit(bline_indices, blinesreal) 
 bwavelengths = binterpolation[1]*pixelf + binterpolation[0]
 ; now plot 
 !P.MULTI = [0, 2, 1, 0, 0]  ; 2x1 arrangement of plots 
 window, retain=2, xsize=1500, ysize=1000
    ; red
 plot, rwavelengths, rrefmiddle, title='Wavelength Calibration: Red', xtitle = 'Wavelength [Angstroms]', ytitle='Relative Intensity', $
     background=255., color=0.
 FOR i=0, n_elements(rpeaks)-1 DO BEGIN 
   xposition = rpeaks[i]
   xyouts, rwavelengths[xposition], rrefmiddle[xposition], i, /data, color=0.,alignment=0.5
 ENDFOR 
    ; blue 
 plot, bwavelengths, brefmiddle, title='Wavelength Calibration: Blue', xtitle = 'Wavelength [Angstroms]', ytitle='Relative Intensity', $
     background=255., color=0.
 FOR i=0, n_elements(bpeaks)-1 DO BEGIN 
   xposition = bpeaks[i]
   xyouts, bwavelengths[xposition], brefmiddle[xposition], i, /data, color=0., alignment=0.5
 ENDFOR 
  
END    ; DISreduce 




FUNCTION CCDread, filename
; reads in fits files from the SPICAM data in holtz apo directory
  file = '/home/holtz/raw/apo/mar16/UT160327/' + filename
  return, readfits(file)
END 




FUNCTION overscansub, image 
; subtracts the overscan of an image from the image 
  print, "overscan is", mean(image[2055:2095,50:1000])  ; this is a safety check
  return, image - mean(image[2055:2095,50:1000])
END




FUNCTION biassub, image, biases
; combines biases and subtracts them from an image
  FORWARD_FUNCTION imcombine
  FORWARD_FUNCTION CCDread
  FOR i=0, n_elements(biases)-1 DO BEGIN
    bias = CCDread(biases[i])
    IF (i EQ 0) THEN biasims=bias ELSE biasims=[[[biasims]],[[bias]]]
  ENDFOR
  combinedbias = imcombine(biasims)
  return, image-combinedbias
END




FUNCTION normalize, image
; normalizes an image (quartz spectrum) by dividing it by max value in image
  return, image / max(image)
END 




FUNCTION flatreduce, filelist, biases
  FORWARD_FUNCTION CCDread
  FORWARD_FUNCTION overscansub
  FORWARD_FUNCTION normalize 
  FORWARD_FUNCTION biassub
  FOR i=0, n_elements(filelist)-1 DO BEGIN   
   reducedfile = overscansub(CCDread(filelist[i]))
   reducedfile = biassub(reducedfile, biases)
   reducedfile = normalize(reducedfile)
   IF (i EQ 0) THEN cube=reducedfile ELSE cube=[[[cube]],[[reducedfile]]]
  ENDFOR 
  return, cube
END 




FUNCTION imcombine, cube 
  med = median(cube, dimension=3)
  return, med
END 