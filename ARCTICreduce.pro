PRO ARCTICreduce, data, flatlist, biaslist 
; ASTR535 labwork 
; first remove overscan
; remove bias frames (avg them then subtract)
; divide by flat field (get flat to have counts of order unity first)
;
; INPUTS: data - the filename of a data image
;         flatlist - a text file containing names of all the flats to use for reduction
;         biaslist - a text file containing names of all the biases to use for reduction
 FORWARD_FUNCTION CCDread
 FORWARD_FUNCTION overscansub
 FORWARD_FUNCTION biassub
 FORWARD_FUNCTION normalize
 FORWARD_FUNCTION flatreduce
 FORWARD_FUNCTION imcombine
 flats = read_csv(flatlist)
 flats = flats[0]
 biases = read_csv(biaslist)
 biases = biases[0]
 dataimage = CCDread(data)
 dataimage = overscansub(dataimage)
 dataimage = biassub(dataimage, biases) 
 help, dataimage    ; this is a safety check
 flatlist = flatreduce(flats, biases)
 megaflat = imcombine(flatlist) 
 dataimage = dataimage/megaflat
 ; then display the data image or something
END


FUNCTION CCDread, filename
; reads in fits files from the SPICAM data in holtz apo directory
  file = '/home/holtz/raw/apo/mar16/UT160327/' + filename
  return, readfits(file)
END 


FUNCTION overscansub, image 
; subtracts the overscan of an image from the image 
  print, "overscan is", mean(image[2057:2093,14:2044])  ; this is a safety check
  return, image - mean(image[2057:2093,14:2044])
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
; normalizes an image (flat) by dividing it by a good mean value from the image
  print, "normalization is", mean(image[400:600,400:600])  ; this is a safety check
  return, image / mean(image[400:600,400:600])
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
  
  
