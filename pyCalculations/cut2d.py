# This file contains a function for retrieving a 2d cut from a 2d plane ( Retrieves the cell ids in that area )

#NOT IMPLEMENTED YET


def cut2d( vlsvReader, xmin, xmax, ymin, ymax, zmin, zmax ):
   ''' Retrieves cell ids from the given 2d cut
       :param vlsvReader         Some VlsvFile with a file open
       :param xmin               The minimum x coordinate of the 2d cut
       :param xmax               The maximum x coordinate of the 2d cut
       :param ymin               The minimum y coordinate of the 2d cut
       :param ymax               The maximum y coordinate of the 2d cut
       :param zmin               The minimum z coordinate of the 2d cut
       :param zmax               The maximum z coordinate of the 2d cut
       NOTE: ONE OF THE PAIRS *MUST* BE ZERO, SO FOR EX. XMIN=0 AND XMAX=0 or YMIN=0 AND YMAX=0 OR ZMIN=0 AND ZMAX=0
   '''
   # Check to make sure that one of the variables is zero
   isOk = False
   # First min/max
   f_min = 0
   f_max = 0
   # Second min/max
   s_min = 0
   s_max = 0
   if (xmax-xmin) == 0:
      isOk = True
   if (ymax-ymin) == 0:
      isOk = True
   if (zmax-zmin) == 0:
      isOk = True
   if isOk == False:
      print "BAD INPUT, ONE OF THE BOUNDS MUST BE LENGTH=0"
      return []
   
   # Retrieve the cell ids:
   #Iterate through the BBOX coordinates:
   minCoordinates = [BBOX[0], BBOX[2], BBOX[4]]
   maxCoordinates = [BBOX[1], BBOX[3], BBOX[5]]
   coordinates = np.array(minCoordinates)



