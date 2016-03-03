PRO CCDreduce, data 
; ASTR535 labwork 
 FORWARD_FUNCTION CCDread
 FORWARD_FUNCTION biassub
 FORWARD_FUNCTION normalize
 FORWARD_FUNCTION flatreduce
 FORWARD_FUNCTION flatcombine
 dataimage = CCDread(data)
 dataimage = biassub(dataimage)
 dataimage = normalize(dataimage) 
 help, dataimage    ; this is a safety check
 flatlist = flatreduce(['flat_r.0013.fits','flat_r.0014.fits','flat_r.0015.fits'])
 megaflat = flatcombine(flatlist) 
END


FUNCTION CCDread, filename
; reads in fits files from the SPICAM data in holtz apo directory
  file = '/home/holtz/raw/apo/dec06/UT061215/' + filename
  return, readfits(file)
END 


FUNCTION biassub, image 
; subtracts the overscan of an image from the image 
  print, "overscan is", mean(image[1037:1072,10:1065])  ; this is a safety check
  return, image - mean(image[1037:1072,10:1065])
END


FUNCTION normalize, image 
; normalizes an image by dividing it by a good mean value from the image
  print, "normalization is", mean(image[400:600,400:600])  ; this is a safety check
  return, image / mean(image[400:600,400:600])
END 


FUNCTION flatreduce, filelist
  FORWARD_FUNCTION CCDread
  FORWARD_FUNCTION biassub
  FORWARD_FUNCTION normalize 
  FOR i=0, n_elements(filelist)-1 DO BEGIN   
   reducedfile = normalize(biassub(CCDread(filelist[i])))
   IF (i EQ 0) THEN cube=reducedfile ELSE cube=[[[cube]],[[reducedfile]]]
  ENDFOR 
  return, cube
END 


FUNCTION flatcombine, cube 
  med = median(cube, dimension=3)
  return, med
END 
  
  