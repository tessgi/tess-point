#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  6 09:14:04 2019

@author: cjburke
"""

import numpy as np
from tess_stars2px import tess_stars2px_function_entry

if __name__ == '__main__':
    # This is an example of using tess_stars2px functionality 
    # from a program rather than the typical command line interface
    # the function example only takes ra and dec [deg] 
    #  Your program should provide these
    #  The other ways of getting coordinates 'TIC MAST query' 'By name' 'file list'
    #  are not supported in the wrapper function.  Just ra and decs
    # Location for pi mensae
    ra = 84.291188 
    dec = -80.469119
    ticid = 261136679 # code doesn't actually use ticid, so this can be
        # any integer you like.  It is included just for convenience of 
        # keeping track of target in the output because
        # for a given target it is not known ahead of time how many
        # output entries a target will have.  Thus the user
        # should make sure to match the output ticid (outID)
        #  with the ticid they are using
    outID, outEclipLong, outEclipLat, outSec, outCam, outCcd, \
            outColPix, outRowPix, scinfo = tess_stars2px_function_entry(
                    ticid, ra, dec)
    for i in range(len(outID)):
        print('{0:d} {1:d} {2:d} {3:d} {4:f} {5:f}'.format(outID[i], outSec[i], \
          outCam[i], outCcd[i], outColPix[i], outRowPix[i]))
        
    # For efficiency purposes if you save scinfo between calls
    #  you will save time in setting up the the telescope fields
    outID, outEclipLong, outEclipLat, outSec, outCam, outCcd, \
            outColPix, outRowPix, scinfo = tess_stars2px_function_entry(
                    ticid, ra, dec, scInfo=scinfo)
    print('Faster to re-use scinfo in repeated calls')
    for i in range(len(outID)):
        print('{0:d} {1:d} {2:d} {3:d} {4:f} {5:f}'.format(outID[i], outSec[i], \
          outCam[i], outCcd[i], outColPix[i], outRowPix[i]))
    
    print('Also accepts multiple targets')
    ra = np.array([219.90085,10.897379], dtype=np.float)
    dec = np.array([-60.835619,-17.986606], dtype=np.float)
    ticid = np.array([0,1], dtype=np.int)
    outID, outEclipLong, outEclipLat, outSec, outCam, outCcd, \
            outColPix, outRowPix, scinfo = tess_stars2px_function_entry(
                    ticid, ra, dec, scInfo=scinfo)
    for i in range(len(outID)):
        print('{0:d} {1:d} {2:d} {3:d} {4:f} {5:f}'.format(outID[i], outSec[i], \
          outCam[i], outCcd[i], outColPix[i], outRowPix[i]))
