
"""
tess_stars2px.py - High precision TESS pointing tool.
Convert target coordinates given in Right Ascension and Declination to
TESS detector pixel coordinates for the prime mission TESS observing 
sectors (Year 1 & 2) and Extendend mission Years 3-4.
Can also query MAST to obtain detector
pixel coordinates for a star by TIC ID only (must be online for this option).

USAGE to display command line arguments:
python tess_stars2px.py -h

AUTHORS: Original programming in C and focal plane geometry solutions 
    by Alan Levine (MIT)
 This python translation by Christopher J. Burke (MIT)
 Testing and focal plane geometry refinements by Michael Fausnaugh & 
         Roland Vanderspek (MIT)
 Testing by Thomas Barclay (NASA Goddard) &
         Jessica Roberts (Univ. of Colorado)
 Sesame queries by Brett Morris (UW)
 Proxy Support added by Dishendra Mishra 

VERSION: 0.6.2

WHAT'S NEW:
    -Sector 46 field update 2021 Oct. 1
    -Too close to edge Warning flag now output in column 11. If a target is within 6 pixels
        of the edge of the science region (edgeWarn==1), then the target is unlikely
        to be assigned a 2minute or 20s aperture. The science pixels range
        in column from 45-2092 and row from 1-2048
    -Year 4 Sectors 40-55 Now Available
    -Now has option to perform an approximate aberration correction.
        Uses astropy GCRS geocentric aberrated frame.  The Earth has velocity of 30km/s
        whereas TESS has velocity relative to Earth of <4km/s.  Thus, correcting
        for Earth aberration is most of the way there. Earth aberrates ~20 arcsec ~ 1 TESS pixel
        Aberration is only done in forward ra,dec-> pix direction. No aberratuib correction is done for inverse pix->ra,dec direction
    -Inverse transform (input Sector, Camera, CCD, pixel Column, pixel Row --> RA and Dec) is now 'analytic' rather than through brute force minimization.  The inverse transform is much faster and much more reliable.
    -MUCH FASTER NOW - skipped rough estimate step which was much much slower
         than just doing the matrix math for position.
    

NOTES:
    -Pointing table is for TESS Year 1 - 4 (Sectors 1-55) 
    -Pointing table is unofficial, and the pointings may change.
    -See https://tess.mit.edu/observations/ for latest TESS pointing table
    -Pointing prediction algorithm is similar to internally at MIT for
        target management.  However, hard coded focal plane geometry is not
        up to date and may contain inaccurate results.
    -Testing shows pointing with this tool should be accurate to better than
        a pixel, but without including aberration effects, ones algorithm
        adopted for centroiding highly assymmetric point-spread function
        at edge of
        camera, and by-eye source location, a 2 pixel accuracy estimate is
        warranted. There is an approximate aberration option now available
    -The output pixel coordinates assume the ds9 convention with
        1,1 being the middle of the lower left corner.
    -No corrections for velocity aberration are calculated.
       Potentially more accurate
        results can be obtained if the target RA and Declination coordinates
        have aberration effects applied.  
    -For proposals to the TESS science office or directors discretionary time,
      please consult the TESS prediction webtool available at
      https://heasarc.gsfc.nasa.gov/cgi-bin/tess/webtess/wtv.py
      for official identification of 'observable' targets.  However,
      if your proposal depends on a single or few targets, then this tool is 
      helpful to further refine the likelihood of the target being available
      on the detectors.
     -The calibrated FFI fits file release at MAST and calibrated by
        NASA Ames SPOC will have WCS information available to
        supplant this code.  The WCS generation is independent of the 
        focal plane geometry model employed in this code and will give
        different results at the pixel level.  However, the WCS information
        is not available until the FFI files are released, whereas
        this code can predict positions in advance of data release.
     -Hard coded focal plane geometry parameters from rfpg5_c1kb.txt

    
TODOS:
    -Time dependent Focal plane geometry

DEPENDENCIES:
    python 3+
    astropy
    numpy

SPECIAL THANKS TO:
    Includes code from the python MAST query examples 
    https://mast.stsci.edu/api/v0/pyex.html

IMPLEMENTATION DETAILS:
    In summary, the code begins with a space craft bore site pointing in RA, 
    Dec, and roll angle.  A series of Euler angle translation matrices
    are calculated based upon the space craft bore site.  Next the target
    coordinates in RA and Dec are translated to the space craft bore site
    frame.  Next, the target coordinates are translated to each of the four
    TESS camera frames.  Once target coordinates are translated to the
    camera frame the radial position of the target relative to the camera
    center is checked to see if it is potentially in the camera field of view.
    If so, the focal plane position is calculated using a radial polynomial
    model with a constant term and terms the even powers (2nd ,4th , and 8th).
    Rotations are applied to convert the on sky positions to the detector
    readout directions.

MIT License    
Copyright (c) 2018 Christopher J Burke

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import numpy as np
import os
import argparse
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
import sys
import datetime
import json
try: # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
    from urllib.parse import urlparse
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve
    from urlparse import urlparse
try: # Python 3.x
    import http.client as httplib 
except ImportError:  # Python 2.x
    import httplib
import scipy.optimize as opt
import base64


class Levine_FPG():
    """Al Levine Focal Plane Geometry Methods
        Translated from starspx6.c
        INPUT:
            sc_ra_dec_roll = numpy array of the SpaceCraft boresite (sc Z-axis)
            ra, dec, and roll [deg]
            The roll angle is in RA, Dec space clockwise relative to the celestial
            pole.  roll angle = 0 [deg] implies space craft X-axis points N celestial (increasing dec)
            roll angle = 90 [deg] implies sc X-axis points towards increasing/decreasing (?) RA
        *** In practice there is a separate fpg file for each of the four cameras ***
        rmat1[3,3] = is the rotation matrix from ra&dec to spacecraft boresite coords
        rmat4[NCAM,3,3] - is the rotation matrix from ra&dec to NCAM coords
    """
    parm_dict_list = [{}, {}, {}, {}]
    NCAM = 4 # Number of Cameras
    NCCD = 4 # Number of CCDs per Camera
    
    def __init__(self, sc_ra_dec_roll=None, fpg_file_list=None):
        self.eulcam = np.zeros((self.NCAM,3), dtype=np.double)
        self.optcon = np.zeros((self.NCAM,6), dtype=np.double)
        self.ccdxy0 = np.zeros((self.NCAM, self.NCCD, 2), dtype=np.double)
        self.pixsz = np.zeros((self.NCAM, self.NCCD, 2), dtype=np.double)
        self.ccdang = np.zeros((self.NCAM, self.NCCD), dtype=np.double)
        self.ccdtilt = np.zeros((self.NCAM, self.NCCD, 2), dtype=np.double)
        self.asymang = np.zeros((self.NCAM,), dtype=np.double)
        self.asymfac = np.zeros((self.NCAM,), dtype=np.double)
        self.rmat1 = np.zeros((3,3), dtype=np.double)
        self.rmat4 = np.zeros((self.NCAM,3,3), dtype=np.double)
        self.havePointing = False
        # Read in the fpg parameter files        
        self.read_all_levine_fpg_files(fpg_file_list)
        # Generate rotation matrices if ra dec and roll values given
        if not sc_ra_dec_roll is None:
            # go from sky to spacecraft
            self.sky_to_sc_mat(sc_ra_dec_roll)
            # Go from spacecraft to each camera's coords
            for icam in range(self.NCAM):
                cureul = self.eulcam[icam,:]
                rmat2 = self.sc_to_cam_mat(cureul)
                self.rmat4[icam] = np.matmul(rmat2, self.rmat1)
            self.havePointing = True
        
    def read_all_levine_fpg_files(self, fpg_file_list=None):                
        default_fpg_file_list = ['fpg_pars.txt-', \
                                 'fpg_pars.txt-', \
                                 'fpg_pars.txt-', \
                                 'fpg_pars.txt-']
        # For each camera read in the separate fpg parameter file
        for icam in range(self.NCAM):
            if fpg_file_list == None:
                fpg_file = default_fpg_file_list[icam]
            else:
                fpg_file = fpg_file_list[icam]
            self.read_levine_fpg_file(icam, fpg_file)
        # We now have parameters for all 4 cameras in the parm_dict_list
        # parse the dictionary values into the working numpy arrays
        for icam in range(self.NCAM):
            pd = self.parm_dict_list[icam]
            self.eulcam[icam][0] = pd['ang1_cam1']
            self.eulcam[icam][1] = pd['ang2_cam1']
            self.eulcam[icam][2] = pd['ang3_cam1']
            self.optcon[icam][0] = pd['fl_cam1']
            self.optcon[icam][1] = pd['opt_coef1_cam1']	
            self.optcon[icam][2] = pd['opt_coef2_cam1']	
            self.optcon[icam][3] = pd['opt_coef3_cam1']	
            self.optcon[icam][4] = pd['opt_coef4_cam1']	
            self.optcon[icam][5] = pd['opt_coef5_cam1']
            self.asymang[icam] = pd['asymang_cam1']
            self.asymfac[icam] = pd['asymfac_cam1']
            for iccd in range(self.NCCD):
                self.ccdxy0[icam][iccd][0] = pd['x0_ccd{0:1d}_cam1'.format(iccd+1)]
                self.ccdxy0[icam][iccd][1] = pd['y0_ccd{0:1d}_cam1'.format(iccd+1)]
                self.pixsz[icam][iccd][0] = pd['pix_x_ccd{0:1d}_cam1'.format(iccd+1)]
                self.pixsz[icam][iccd][1] = pd['pix_y_ccd{0:1d}_cam1'.format(iccd+1)]
                self.ccdang[icam][iccd] = pd['ang_ccd{0:1d}_cam1'.format(iccd+1)]
                self.ccdtilt[icam][iccd][0] = pd['tilt_x_ccd{0:1d}_cam1'.format(iccd+1)]
                self.ccdtilt[icam][iccd][1] = pd['tilt_y_ccd{0:1d}_cam1'.format(iccd+1)]
        
            
    def read_levine_fpg_file(self, icam, fpg_file):
        gotParm = False
        parm_dict = {}
        if os.path.isfile(fpg_file):        
            try:    
                fpin = open(fpg_file, 'r')
                # Read in parameters 
                dtypeseq = ['U20','i4','f16']
                dataBlock = np.genfromtxt(fpin, dtype=dtypeseq)
                parm_keys = dataBlock['f0']
                parm_fitted_flags = dataBlock['f1']
                parm_values = dataBlock['f2']
                # Now build dictionary of the parameters
                for i in range(len(parm_keys)):
                    parm_dict[parm_keys[i]] = parm_values[i]
                self.parm_dict_list[icam] = parm_dict
                gotParm = True
                print('Successful Focal Plane Geometry Read From {0}'.format(fpg_file))
            except:
                print('Could not open {0}!  Using Hard-coded Focal Plane Geometry from Levine_FPG read_levine_fpg_file()'.format(fpg_file))
        # If anything goes wrong with reading in parameters revert to hard coded version
        # or file was never given and default_fpg_file does not exist
        if not gotParm:
            #print('Using Hard-coded Focal Plane Geometry from Levine_FPG read_levine_fpg_file')
            # *** For now this hard code is just a filler need to actually fill in values for all cameras separately
            # to prepare parameters for dictionary
            # awk -v q="'" -v d=":" '{print q $1 q d $3 ",\"}' rfpg5_c1kb.txt
            if icam == 0:
                parm_dict = {'ang1_cam1' : 0.101588, \
                             'ang2_cam1' : -36.022035, \
                             'ang3_cam1' :  90.048315, \
                             'fl_cam1' :   145.948116, \
                             'opt_coef1_cam1' :   1.00000140, \
                             'opt_coef2_cam1' :  0.24779006, \
                             'opt_coef3_cam1' :  -0.22681254, \
                             'opt_coef4_cam1' :  10.78243356, \
                             'opt_coef5_cam1' :  -34.97817276, \
                             'asymang_cam1' :  0.00000000, \
                             'asymfac_cam1' :  1.00000000, \
                             'x0_ccd1_cam1' :  31.573417, \
                             'y0_ccd1_cam1' :  31.551637, \
                             'pix_x_ccd1_cam1' :  0.015000, \
                             'pix_y_ccd1_cam1' :  0.015000, \
                             'ang_ccd1_cam1' :  179.980833, \
                             'tilt_x_ccd1_cam1' :  0.000000, \
                             'tilt_y_ccd1_cam1' :  0.000000, \
                             'x0_ccd2_cam1' :  -0.906060, \
                             'y0_ccd2_cam1' :  31.536148, \
                             'pix_x_ccd2_cam1' :  0.015000, \
                             'pix_y_ccd2_cam1' :  0.015000, \
                             'ang_ccd2_cam1' :  180.000000, \
                             'tilt_x_ccd2_cam1' :  0.000000, \
                             'tilt_y_ccd2_cam1' :  0.000000, \
                             'x0_ccd3_cam1' :  -31.652818, \
                             'y0_ccd3_cam1' :  -31.438350, \
                             'pix_x_ccd3_cam1' :  0.015000, \
                             'pix_y_ccd3_cam1' :  0.015000, \
                             'ang_ccd3_cam1' :  -0.024851, \
                             'tilt_x_ccd3_cam1' :  0.000000, \
                             'tilt_y_ccd3_cam1' :  0.000000, \
                             'x0_ccd4_cam1' :  0.833161, \
                             'y0_ccd4_cam1' :  -31.458180, \
                             'pix_x_ccd4_cam1' :  0.015000, \
                             'pix_y_ccd4_cam1' :  0.015000, \
                             'ang_ccd4_cam1' :  0.001488, \
                             'tilt_x_ccd4_cam1' :  0.000000, \
                             'tilt_y_ccd4_cam1' :  0.000000}
                        

            if icam == 1:
                parm_dict =  {'ang1_cam1':-0.179412,\
                              'ang2_cam1':-12.017260,\
                              'ang3_cam1':90.046500,\
                              'fl_cam1':145.989933,\
                              'opt_coef1_cam1':1.00000140,\
                              'opt_coef2_cam1':0.24069345,\
                              'opt_coef3_cam1':0.15391120,\
                              'opt_coef4_cam1':4.05433503,\
                              'opt_coef5_cam1':3.43136895,\
                              'asymang_cam1':0.00000000,\
                              'asymfac_cam1':1.00000000,\
                              'x0_ccd1_cam1':31.653635,\
                              'y0_ccd1_cam1':31.470291,\
                              'pix_x_ccd1_cam1':0.015000,\
                              'pix_y_ccd1_cam1':0.015000,\
                              'ang_ccd1_cam1':180.010890,\
                              'tilt_x_ccd1_cam1':0.000000,\
                              'tilt_y_ccd1_cam1':0.000000,\
                              'x0_ccd2_cam1':-0.827405,\
                              'y0_ccd2_cam1':31.491388,\
                              'pix_x_ccd2_cam1':0.015000,\
                              'pix_y_ccd2_cam1':0.015000,\
                              'ang_ccd2_cam1':180.000000,\
                              'tilt_x_ccd2_cam1':0.000000,\
                              'tilt_y_ccd2_cam1':0.000000,\
                              'x0_ccd3_cam1':-31.543794,\
                              'y0_ccd3_cam1':-31.550699,\
                              'pix_x_ccd3_cam1':0.015000,\
                              'pix_y_ccd3_cam1':0.015000,\
                              'ang_ccd3_cam1':-0.006624,\
                              'tilt_x_ccd3_cam1':0.000000,\
                              'tilt_y_ccd3_cam1':0.000000,\
                              'x0_ccd4_cam1':0.922834,\
                              'y0_ccd4_cam1':-31.557268,\
                              'pix_x_ccd4_cam1':0.015000,\
                              'pix_y_ccd4_cam1':0.015000,\
                              'ang_ccd4_cam1':-0.015464,\
                              'tilt_x_ccd4_cam1':0.000000,\
                              'tilt_y_ccd4_cam1':0.000000}
                
            if icam == 2:
                parm_dict = {'ang1_cam1':0.066596,\
                             'ang2_cam1':12.007750,\
                             'ang3_cam1':-89.889085,\
                             'fl_cam1':146.006602,\
                             'opt_coef1_cam1':1.00000140,\
                             'opt_coef2_cam1':0.23452229,\
                             'opt_coef3_cam1':0.33552009,\
                             'opt_coef4_cam1':1.92009863,\
                             'opt_coef5_cam1':12.48880182,\
                             'asymang_cam1':0.00000000,\
                             'asymfac_cam1':1.00000000,\
                             'x0_ccd1_cam1':31.615486,\
                             'y0_ccd1_cam1':31.413644,\
                             'pix_x_ccd1_cam1':0.015000,\
                             'pix_y_ccd1_cam1':0.015000,\
                             'ang_ccd1_cam1':179.993948,\
                             'tilt_x_ccd1_cam1':0.000000,\
                             'tilt_y_ccd1_cam1':0.000000,\
                             'x0_ccd2_cam1':-0.832993,\
                             'y0_ccd2_cam1':31.426621,\
                             'pix_x_ccd2_cam1':0.015000,\
                             'pix_y_ccd2_cam1':0.015000,\
                             'ang_ccd2_cam1':180.000000,\
                             'tilt_x_ccd2_cam1':0.000000,\
                             'tilt_y_ccd2_cam1':0.000000,\
                             'x0_ccd3_cam1':-31.548296,\
                             'y0_ccd3_cam1':-31.606976,\
                             'pix_x_ccd3_cam1':0.015000,\
                             'pix_y_ccd3_cam1':0.015000,\
                             'ang_ccd3_cam1':0.000298,\
                             'tilt_x_ccd3_cam1':0.000000,\
                             'tilt_y_ccd3_cam1':0.000000,\
                             'x0_ccd4_cam1':0.896018,\
                             'y0_ccd4_cam1':-31.569542,\
                             'pix_x_ccd4_cam1':0.015000,\
                             'pix_y_ccd4_cam1':0.015000,\
                             'ang_ccd4_cam1':-0.006464,\
                             'tilt_x_ccd4_cam1':0.000000,\
                             'tilt_y_ccd4_cam1':0.000000}
                
            if icam == 3:
                parm_dict = {'ang1_cam1':0.030756,\
                             'ang2_cam1':35.978116,\
                             'ang3_cam1':-89.976802,\
                             'fl_cam1':146.039793,\
                             'opt_coef1_cam1':1.00000140,\
                             'opt_coef2_cam1':0.23920416,\
                             'opt_coef3_cam1':0.13349450,\
                             'opt_coef4_cam1':4.77768896,\
                             'opt_coef5_cam1':-1.75114744,\
                             'asymang_cam1':0.00000000,\
                             'asymfac_cam1':1.00000000,\
                             'x0_ccd1_cam1':31.575820,\
                             'y0_ccd1_cam1':31.316510,\
                             'pix_x_ccd1_cam1':0.015000,\
                             'pix_y_ccd1_cam1':0.015000,\
                             'ang_ccd1_cam1':179.968217,\
                             'tilt_x_ccd1_cam1':0.000000,\
                             'tilt_y_ccd1_cam1':0.000000,\
                             'x0_ccd2_cam1':-0.890877,\
                             'y0_ccd2_cam1':31.363511,\
                             'pix_x_ccd2_cam1':0.015000,\
                             'pix_y_ccd2_cam1':0.015000,\
                             'ang_ccd2_cam1':180.000000,\
                             'tilt_x_ccd2_cam1':0.000000,\
                             'tilt_y_ccd2_cam1':0.000000,\
                             'x0_ccd3_cam1':-31.630470,\
                             'y0_ccd3_cam1':-31.716942,\
                             'pix_x_ccd3_cam1':0.015000,\
                             'pix_y_ccd3_cam1':0.015000,\
                             'ang_ccd3_cam1':-0.024359,\
                             'tilt_x_ccd3_cam1':0.000000,\
                             'tilt_y_ccd3_cam1':0.000000,\
                             'x0_ccd4_cam1':0.824159,\
                             'y0_ccd4_cam1':-31.728751,\
                             'pix_x_ccd4_cam1':0.015000,\
                             'pix_y_ccd4_cam1':0.015000,\
                             'ang_ccd4_cam1':-0.024280,\
                             'tilt_x_ccd4_cam1':0.000000,\
                             'tilt_y_ccd4_cam1':0.000000}
                
            self.parm_dict_list[icam] = parm_dict

    def sky_to_sc_mat(self, sc_ra_dec_roll):
        """Calculate the rotation matrix that will convert a vector in ra&dec
            into the spacecraft boresite frame
        """
        deg2rad = np.pi / 180.0
        # Define the 3 euler angles of rotation
        xeul = np.zeros((3,), dtype=np.double)
        xeul[0] = deg2rad * sc_ra_dec_roll[0]
        xeul[1] = np.pi/2.0 - deg2rad*sc_ra_dec_roll[1]
        xeul[2] = deg2rad * sc_ra_dec_roll[2] + np.pi
        # Generate the rotation matrix from the 3 euler angles
        self.rmat1 = self.eulerm323(xeul)

    def sc_to_cam_mat(self, eul):
        """Calculate the rotation matrix that will convert a vector in spacecraft
            into the a camera's coords
        """
        deg2rad = np.pi / 180.0
        # Generate the rotation matrix from the 3 euler angles
        xeul = deg2rad * eul
        return self.eulerm323(xeul)
        
    def eulerm323(self, eul):
        mat1 = self.rotm1(2, eul[0])
        mat2 = self.rotm1(1, eul[1])
        mata = np.matmul(mat2, mat1)
        mat1 = self.rotm1(2, eul[2])
        rmat = np.matmul(mat1, mata)
        return rmat
        
    def rotm1(self, ax, ang):
        mat = np.zeros((3,3), dtype=np.double)
        n1 = ax
        n2 = np.mod((n1+1), 3)
        n3 = np.mod((n2+1), 3)
        sinang = np.sin(ang)
        cosang = np.cos(ang) 
        mat[n1][n1] = 1.0
        mat[n2][n2] = cosang
        mat[n3][n3] = cosang
        mat[n2][n3] = sinang
        mat[n3][n2] = -sinang
        return mat

    def sphereToCart(self, ras, decs):
        """ Convert 3d spherical coordinates to cartesian
        """
        deg2rad = np.pi / 180.0
        rarads = deg2rad * ras
        decrads = deg2rad * decs
        sinras = np.sin(rarads)
        cosras = np.cos(rarads)
        sindecs = np.sin(decrads)
        cosdecs = np.cos(decrads)
        vec0s = cosras * cosdecs
        vec1s = sinras * cosdecs
        vec2s = sindecs
        return vec0s, vec1s, vec2s

    def cartToSphere(self, vec):
        ra = 0.0
        dec = 0.0
        norm = np.sqrt(np.sum(vec*vec))
        if (norm > 0.0):
            dec = np.arcsin(vec[2] / norm)
            if (not vec[0] == 0.0) or (not vec[1] == 0.0):
                ra = np.arctan2(vec[1], vec[0])
                ra = np.mod(ra, 2.0*np.pi)
        return ra, dec

    def star_in_fov(self, lng, lat):
        deg2rad = np.pi / 180.0
        inView = False
        if lat > 70.0:
            vec0, vec1, vec2 = self.sphereToCart(lng, lat)
            vec = np.array([vec0, vec1, vec2], dtype=np.double)
            norm = np.sqrt(np.sum(vec*vec))
            if norm > 0.0:
                vec = vec / norm
                xlen = np.abs(np.arctan(vec[0]/vec[2]))
                ylen = np.abs(np.arctan(vec[1]/vec[2]))
                if (xlen <= (12.5 * deg2rad)) and (ylen <= (12.5 * deg2rad)):
                    inView = True
        return inView

    def optics_fp(self, icam, lng_deg, lat_deg):
        deg2rad = np.pi / 180.0
        thetar = np.pi / 2.0 - (lat_deg * deg2rad)
        tanth = np.tan(thetar)
        cphi = np.cos(deg2rad*lng_deg)
        sphi = np.sin(deg2rad*lng_deg)
        rfp0 = self.optcon[icam][0]*tanth
        noptcon = len(self.optcon[icam])
        ii = np.arange(1, noptcon)
        rfp = np.sum(self.optcon[icam][1:] * np.power(tanth, 2.0*(ii-1)))
        xytmp = np.zeros((2,), dtype=np.double)
        xytmp[0] = -cphi*rfp0*rfp
        xytmp[1] = -sphi*rfp0*rfp
        return self.make_az_asym(icam, xytmp)

    def fp_optics(self, icam, xyfp):
        deg2rad = np.pi/ 180.0
        xy = self.rev_az_asym(icam, xyfp)
        rfp_times_rfp0 = np.sqrt(xy[0]*xy[0] + xy[1]*xy[1])
        phirad = np.arctan2(-xy[1], -xy[0])
        phideg = phirad / deg2rad
        if (phideg < 0.0):
            phideg += 360.0
        thetarad = self.tanth_of_r(icam, rfp_times_rfp0)
        thetadeg = thetarad / deg2rad
        lng_deg = phideg
        lat_deg = 90.0 - thetadeg
        return lng_deg, lat_deg
    
    def r_of_tanth(self, icam, z):
        tanth = np.tan(z)
        rfp0 = self.optcon[icam][0]*tanth
        noptcon = len(self.optcon[icam])
        ii = np.arange(1, noptcon)
        rfp = np.sum(self.optcon[icam][1:] * np.power(tanth, 2.0*(ii-1)))
        return rfp0*rfp
        
    def tanth_of_r(self, icam, rprp0):

        if np.abs(rprp0) > 1.0e-10:
            c0 = self.optcon[icam][0]
            zi = np.arctan(np.sqrt(rprp0) / c0)
            def minFunc(z, icam, rp):
                rtmp = self.r_of_tanth(icam, z)
                return (rtmp - rprp0)*(rtmp - rprp0)
            
            optResult = opt.minimize(minFunc, [zi], \
                                     args=(icam, rprp0), method='Nelder-Mead', \
                                     tol=1.0e-10, \
                                     options={'maxiter':500})
            #print(optResult)
            return optResult.x[0]
        else:
            return 0.0

    def make_az_asym(self, icam, xy):
        xyp = self.xyrotate(self.asymang[icam], xy)
        xypa = np.zeros_like(xyp)
        xypa[0] = self.asymfac[icam] * xyp[0]
        xypa[1] = xyp[1]
        xyout = self.xyrotate(-self.asymang[icam], xypa)
        return xyout

    def rev_az_asym(self, icam, xyin):
        xyp = self.xyrotate(self.asymang[icam], xyin)
        xypa = np.zeros_like(xyp)
        xypa[0] = xyp[0] / self.asymfac[icam]
        xypa[1] = xyp[1]
        xyout = self.xyrotate(-self.asymang[icam], xypa)
        return xyout
        
    def xyrotate(self, angle_deg, xin):
        deg2rad = np.pi / 180.0
        ca = np.cos(deg2rad * angle_deg)
        sa = np.sin(deg2rad * angle_deg)
        xyout = np.zeros_like(xin)
        xyout[0] = ca*xin[0] + sa*xin[1]
        xyout[1] = -sa*xin[0] + ca*xin[1]
        return xyout

    def mm_to_pix(self, icam, xy):
        """Convert focal plane to pixel location also need to add in the
            auxillary pixels added into FFIs
        """
        CCDWD_T=2048
        CCDHT_T=2058
        ROWA=44
        ROWB=44
        COLDK_T=20
        xya = np.copy(xy)
        xyb = np.zeros_like(xya)
        ccdpx = np.zeros_like(xya)
        fitpx = np.zeros_like(xya)
        if xya[0] >= 0.0:
            if xya[1] >= 0.0:
                iccd = 0
                xyb[0] = xya[0] - self.ccdxy0[icam][iccd][0]
                xyb[1] = xya[1] - self.ccdxy0[icam][iccd][1]
                xyccd = self.xyrotate(self.ccdang[icam][iccd], xyb)
                ccdpx[0] = (xyccd[0] / self.pixsz[icam][iccd][0]) - 0.5
                ccdpx[1] = (xyccd[1] / self.pixsz[icam][iccd][1]) - 0.5
                fitpx[0] = (CCDWD_T - ccdpx[0]) + CCDWD_T + 2*ROWA + ROWB - 1.0
                fitpx[1] = (CCDHT_T - ccdpx[1]) + CCDHT_T + 2*COLDK_T - 1.0
            else:
                iccd = 3
                xyb[0] = xya[0] - self.ccdxy0[icam][iccd][0]
                xyb[1] = xya[1] - self.ccdxy0[icam][iccd][1]
                xyccd = self.xyrotate(self.ccdang[icam][iccd], xyb)
                ccdpx[0] = (xyccd[0] / self.pixsz[icam][iccd][0]) - 0.5
                ccdpx[1] = (xyccd[1] / self.pixsz[icam][iccd][1]) - 0.5
                fitpx[0] = ccdpx[0] + CCDWD_T + 2*ROWA + ROWB
                fitpx[1] = ccdpx[1]
        else:
            if xya[1] >= 0.0:
                iccd = 1
                xyb[0] = xya[0] - self.ccdxy0[icam][iccd][0]
                xyb[1] = xya[1] - self.ccdxy0[icam][iccd][1]
                xyccd = self.xyrotate(self.ccdang[icam][iccd], xyb)
                ccdpx[0] = (xyccd[0] / self.pixsz[icam][iccd][0]) - 0.5
                ccdpx[1] = (xyccd[1] / self.pixsz[icam][iccd][1]) - 0.5
                fitpx[0] = (CCDWD_T - ccdpx[0]) + ROWA - 1.0
                fitpx[1] = (CCDHT_T - ccdpx[1]) + CCDHT_T + 2*COLDK_T - 1.0
            else:
                iccd = 2
                xyb[0] = xya[0] - self.ccdxy0[icam][iccd][0]
                xyb[1] = xya[1] - self.ccdxy0[icam][iccd][1]
                xyccd = self.xyrotate(self.ccdang[icam][iccd], xyb)
                ccdpx[0] = (xyccd[0] / self.pixsz[icam][iccd][0]) - 0.5
                ccdpx[1] = (xyccd[1] / self.pixsz[icam][iccd][1]) - 0.5
                fitpx[0] = ccdpx[0] + ROWA
                fitpx[1] = ccdpx[1]
                
        return iccd, ccdpx, fitpx

    def mm_to_pix_single_ccd(self, icam, xy, iccd):
        """Convert focal plane to pixel location also need to add in the
            auxillary pixels added into FFIs
        """
        CCDWD_T=2048
        CCDHT_T=2058
        ROWA=44
        ROWB=44
        COLDK_T=20
        xya = np.copy(xy)
        xyb = np.zeros_like(xya)
        ccdpx = np.zeros_like(xya)
        fitpx = np.zeros_like(xya)
        if iccd == 0:
            xyb[0] = xya[0] - self.ccdxy0[icam][iccd][0]
            xyb[1] = xya[1] - self.ccdxy0[icam][iccd][1]
            xyccd = self.xyrotate(self.ccdang[icam][iccd], xyb)
            ccdpx[0] = (xyccd[0] / self.pixsz[icam][iccd][0]) - 0.5
            ccdpx[1] = (xyccd[1] / self.pixsz[icam][iccd][1]) - 0.5
            fitpx[0] = (CCDWD_T - ccdpx[0]) + CCDWD_T + 2*ROWA + ROWB - 1.0
            fitpx[1] = (CCDHT_T - ccdpx[1]) + CCDHT_T + 2*COLDK_T - 1.0
        if iccd == 3:
            xyb[0] = xya[0] - self.ccdxy0[icam][iccd][0]
            xyb[1] = xya[1] - self.ccdxy0[icam][iccd][1]
            xyccd = self.xyrotate(self.ccdang[icam][iccd], xyb)
            ccdpx[0] = (xyccd[0] / self.pixsz[icam][iccd][0]) - 0.5
            ccdpx[1] = (xyccd[1] / self.pixsz[icam][iccd][1]) - 0.5
            fitpx[0] = ccdpx[0] + CCDWD_T + 2*ROWA + ROWB
            fitpx[1] = ccdpx[1]
        if iccd == 1:
            xyb[0] = xya[0] - self.ccdxy0[icam][iccd][0]
            xyb[1] = xya[1] - self.ccdxy0[icam][iccd][1]
            xyccd = self.xyrotate(self.ccdang[icam][iccd], xyb)
            ccdpx[0] = (xyccd[0] / self.pixsz[icam][iccd][0]) - 0.5
            ccdpx[1] = (xyccd[1] / self.pixsz[icam][iccd][1]) - 0.5
            fitpx[0] = (CCDWD_T - ccdpx[0]) + ROWA - 1.0
            fitpx[1] = (CCDHT_T - ccdpx[1]) + CCDHT_T + 2*COLDK_T - 1.0
        if iccd == 2:
            xyb[0] = xya[0] - self.ccdxy0[icam][iccd][0]
            xyb[1] = xya[1] - self.ccdxy0[icam][iccd][1]
            xyccd = self.xyrotate(self.ccdang[icam][iccd], xyb)
            ccdpx[0] = (xyccd[0] / self.pixsz[icam][iccd][0]) - 0.5
            ccdpx[1] = (xyccd[1] / self.pixsz[icam][iccd][1]) - 0.5
            fitpx[0] = ccdpx[0] + ROWA
            fitpx[1] = ccdpx[1]
                
        return ccdpx, fitpx

    def pix_to_mm_single_ccd(self, icam, ccdpx, iccd):
        """convert pixel to mm focal plane position
        """
        xyccd = np.zeros_like(ccdpx)
        xyccd[0] = (ccdpx[0] + 0.5) * self.pixsz[icam][iccd][0]
        xyccd[1] = (ccdpx[1] + 0.5) * self.pixsz[icam][iccd][1]
        xyb = self.xyrotate(-self.ccdang[icam][iccd], xyccd)
        xya = np.zeros_like(xyb)
        xya[0] = xyb[0] + self.ccdxy0[icam][iccd][0]
        xya[1] = xyb[1] + self.ccdxy0[icam][iccd][1]
                
        return xya
        
    def radec2pix(self, ras, decs):
        """ After the rotation matrices are defined to the actual
            ra and dec to pixel coords mapping
        """
        nStar = len(ras)
        inCamera = np.array([], dtype=np.int)
        ccdNum = np.array([], dtype=np.int)
        fitsxpos = np.array([], dtype=np.double)
        fitsypos = np.array([], dtype=np.double)
        ccdxpos = np.array([], dtype=np.double)
        ccdypos = np.array([], dtype=np.double)
        
        deg2rad = np.pi / 180.0
        if self.havePointing == True:
            # Convert ra and dec spherical coords to cartesian
            vec0s, vec1s, vec2s = self.sphereToCart(ras, decs)
            for i in range(nStar):
                curVec = np.array([vec0s[i], vec1s[i], vec2s[i]], dtype=np.double)
                # Find the new vector in all cameras
                for j in range(self.NCAM):
                    # Do the rotation from ra dec coords to camera coords
                    camVec = np.matmul(self.rmat4[j], curVec)
                    # Get the longitude and latitude of camera coords position
                    lng, lat = self.cartToSphere(camVec)
                    lng = lng / deg2rad
                    lat = lat / deg2rad
                    if self.star_in_fov(lng, lat):
                        # Get the xy focal plane position in mm
                        xyfp = self.optics_fp(j, lng, lat)
                        # Convert mm to pixels
                        iccd, ccdpx, fitpx = self.mm_to_pix(j, xyfp)
                        inCamera = np.append(inCamera, j+1) # Als code is base 0 convert to base 1
                        ccdNum = np.append(ccdNum, iccd+1) # ""
                        fitsxpos = np.append(fitsxpos, fitpx[0])
                        fitsypos = np.append(fitsypos, fitpx[1])
                        ccdxpos = np.append(ccdxpos, ccdpx[0])
                        ccdypos = np.append(ccdypos, ccdpx[1])

        else:
            print('Spacecraft Pointing Not specified!')
        
        return inCamera, ccdNum, fitsxpos, fitsypos, ccdxpos, ccdypos

    def radec2pix_nocheck_single(self, ras, decs, cam, iccd):
        """
            ra and dec to pixel coords mapping
            With no checks and assuming a single target and detector
            Supports minimizing for reverse mode
        """
        deg2rad = np.pi / 180.0
        # Convert ra and dec spherical coords to cartesian
        vec0s, vec1s, vec2s = self.sphereToCart(ras, decs)
        curVec = np.array([vec0s, vec1s, vec2s], dtype=np.double)
        j = cam
        # Do the rotation from ra dec coords to camera coords
        camVec = np.matmul(self.rmat4[j], curVec)
        # Get the longitude and latitude of camera coords position
        lng, lat = self.cartToSphere(camVec)
        lng = lng / deg2rad
        lat = lat / deg2rad
        # Get the xy focal plane position in mm
        xyfp = self.optics_fp(j, lng, lat)
        # Convert mm to pixels
        ccdpx, fitpx = self.mm_to_pix_single_ccd(j, xyfp, iccd)
        ccdNum = iccd+1
        fitsxpos = fitpx[0]
        fitsypos = fitpx[1]
        ccdxpos = ccdpx[0]
        ccdypos = ccdpx[1]
        return ccdNum, fitsxpos, fitsypos, ccdxpos, ccdypos, lat
        
    def pix2radec_nocheck_single(self, cam, iccd, ccdpx):
        """
            Reverse the transform going from pixel coords to Ra & Dec
        """
        deg2rad = np.pi / 180.0
        # Convert pixels to mm
        xyfp = self.pix_to_mm_single_ccd(cam, ccdpx, iccd)
        lng_deg, lat_deg = self.fp_optics(cam, xyfp)
        vcam0, vcam1, vcam2 = self.sphereToCart(lng_deg, lat_deg)
        vcam = np.array([vcam0, vcam1, vcam2], dtype=np.double)
        curVec = np.matmul(np.transpose(self.rmat4[cam]), vcam)
        ra, dec = self.cartToSphere(curVec)
        
        return ra/deg2rad, dec/deg2rad
        
class TESS_Spacecraft_Pointing_Data:
    #Hard coded spacecraft pointings by Sector
    # When adding sectors the arg2 needs to end +1 from sector
    #  due to the np.arange function ending at arg2-1
    sectors = np.arange(1,56, dtype=np.int)

    # Arrays are broken up into the following sectors:
    # Line 1: Sectors 1-5
    # Line 2: Secotrs 6-9
    # Line 3: Sectors 10-13
    # Line 4: Sectors 14-17
    # Line 5: Sectors 18-22
    # Line 6: Sectors 23-26
    # Line 7: Sectors 27-30
    # Line 8: Sectors 31-34
    # Line 9: Sectors 35-38
    # Line 10: Sectors 39
    # Line 11: Sectors 40-43
    # Line 12: S 44-47
    # Line 13: S 48-51
    # Line 14: S 52-55
    ### NOTE IF you add Sectors be sure to update the allowed range
    ### for sectors in argparse arguments!!!
    ras = np.array([352.6844,16.5571,36.3138,55.0070,73.5382, \
                    92.0096,110.2559,128.1156,145.9071,\
                    165.0475,189.1247,229.5885,298.6671, \
                    276.7169,280.3985,282.4427,351.2381,\
                    16.1103,60.2026,129.3867,171.7951,197.1008,\
                    217.2879,261.4516,265.6098,270.1381,\
                    326.8525,357.2944,18.9190,38.3564,\
                    57.6357,77.1891,96.5996,115.2951,\
                    133.2035,150.9497,170.2540,195.7176,\
                    242.1981, \
                    273.0766, 277.6209, 13.0140, 49.5260, \
                    89.6066, 130.2960, 157.6997, 143.3807,\
                    179.4254, 202.6424, 221.8575, 239.4257, \
                    266.3618, 270.8126, 290.1210, 307.8655], dtype=np.float)
            
    decs = np.array([-64.8531,-54.0160,-44.2590,-36.6420,-31.9349, \
                     -30.5839,-32.6344,-37.7370,-45.3044,\
                     -54.8165,-65.5369,-75.1256,-76.3281,\
                     62.4756,64.0671,66.1422,57.8456, \
                     67.9575,76.2343,75.2520,65.1924,53.7434, \
                     43.8074,63.1181,61.9383,61.5637,\
                     -72.4265,-63.0056,-52.8296,-43.3178,\
                     -35.7835,-31.3957,-30.7848,-33.7790,\
                     -39.6871,-47.7512,-57.3725,-67.8307,\
                     -76.3969, \
                     61.7450, 62.7640, 6.3337, 18.9737,\
                     24.1343, 19.0181, 10.0922, 73.1125, \
                     62.1038, 50.9532, 41.7577, 35.2333, \
                     61.8190, 61.5761, 32.6073, 37.6464], dtype=np.float) 

    rolls = np.array([-137.8468,-139.5665,-146.9616,-157.1698,-168.9483, \
                      178.6367,166.4476,155.3091,145.9163,\
                      139.1724,138.0761,153.9773,-161.0622,\
                      32.2329,55.4277,79.4699,41.9686,\
                      40.5453,19.6463,334.5689,317.9495,319.6992,\
                      327.4246,317.2624,339.5293,0.6038,\
                      214.5061,222.5216,219.7970,212.0441,\
                      201.2334,188.6263,175.5369,163.1916,\
                      152.4006,143.7306,138.1685,139.3519,\
                      161.5986,\
                      14.1539, 37.2224, 292.8009, 284.9617,\
                      270.1557, 255.0927, 248.4063, 327.1020,\
                      317.4166, 321.3516, 329.7340, 339.8650,\
                      343.1429, 3.6838, 13.4565, 24.5369], dtype=np.float) 

    midtimes = np.array([ 2458339.652778, 2458368.593750, 2458396.659722, 2458424.548611, 2458451.548611,\
                         2458478.104167, 2458504.697917, 2458530.256944, 2458556.722222,\
                         2458582.760417, 2458610.774306, 2458640.031250, 2458668.618056,\
                         2458697.336806, 2458724.934028, 2458751.649306, 2458777.722222,\
                         2458803.440972, 2458828.958333, 2458856.388889, 2458884.916667, 2458913.565972,\
                         2458941.829861, 2458969.263889, 2458996.909722, 2459023.107639,\
                         2459049.145833, 2459075.166667, 2459102.319444, 2459130.201389,\
                         2459158.854167, 2459186.940972, 2459215.427083, 2459241.979167,\
                         2459268.579861, 2459295.301177, 2459322.577780, 2459349.854382,\
                         2459377.130985,\
                         2459404.407588, 2459431.684191, 2459458.960794, 2459486.237397,\
                         2459513.514000, 2459540.790602, 2459568.067205, 2459595.343808,\
                         2459622.620411, 2459649.897014, 2459677.173617, 2459704.450219,\
                         2459731.726822, 2459759.003425, 2459786.280028, 2459813.556631], dtype=np.float)

    camSeps = np.array([36.0, 12.0, 12.0, 36.0], dtype=np.float)
    

    def __init__(self, trySector=None, fpgParmFileList=None):
        # Convert S/C boresite pointings to ecliptic coords for each camera
        # If trySector is set only keep the single requested sector
        if not trySector is None:
            idx = np.where(self.sectors == trySector)[0]
            self.sectors = self.sectors[idx]
            self.ras = self.ras[idx]
            self.decs = self.decs[idx]
            self.rolls = self.rolls[idx]
            self.midtimes = self.midtimes[idx]
        nPoints = len(self.sectors)
        self.camRa = np.zeros((4, nPoints), dtype=np.float)
        self.camDec = np.zeros((4, nPoints), dtype=np.float)
        # Convert S/C boresite ra and dec to camera ra and dec
        for iPnt in range(nPoints):
            curra = self.ras[iPnt]
            curdec = self.decs[iPnt]
            curroll = self.rolls[iPnt]
            camposangs = np.array([180.0-curroll, 180.0-curroll, \
                                   360.0-curroll, 360.0-curroll])
            camposangs = np.mod(camposangs, 360.0)
            for iCam in range(4):
                # Need to correct s/c roll to posang
                pang = camposangs[iCam]
                camra, camdec = get_radec_from_posangsep(curra, curdec, \
                                            pang, self.camSeps[iCam])
                self.camRa[iCam,iPnt] = camra 
                self.camDec[iCam,iPnt] = camdec 
                # Just for testing camera coords 
                # compare to published values
#                print('{:d} {:d} {:f} {:f}'.format(self.sectors[iPnt],iCam+1,\
#                         self.camRa[iCam,iPnt], self.camDec[iCam,iPnt]))
        # For every pointing make a Levine pointing class object
        self.fpgObjs = []
        fpg_file_list=None
        if not fpgParmFileList is None:
            fpg_file_list=fpgParmFileList
        for iPnt in range(nPoints):
            sc_ra_dec_roll =  np.array([self.ras[iPnt], self.decs[iPnt], self.rolls[iPnt]])
            self.fpgObjs.append(Levine_FPG(sc_ra_dec_roll, fpg_file_list=fpg_file_list))

def get_radec_from_posangsep(ra, dec, pa, sep):
    deg2rad = np.pi/180.0
    rad2deg = 180.0/np.pi
    twopi = 2.0*np.pi
    pidtwo = np.pi/2.0
    rar = ra*deg2rad
    decr = dec*deg2rad
    par = pa*deg2rad
    sepr = sep*deg2rad
    c = pidtwo - decr
    bigB = par
    a = sepr
    b = np.arccos((np.cos(c)*np.cos(a) + np.sin(c)*np.sin(a)*np.cos(bigB)))
    newdec = pidtwo - b
    delalp = np.arccos(np.min([(np.cos(sepr)-np.sin(decr)*np.sin(newdec))/(np.cos(decr)*np.cos(newdec)),1.0]))
    if pa > 180.0:
        newra = rar - delalp
    else:
        newra = rar + delalp
#    print(pa, newra*rad2deg, rar*rad2deg, delalp*rad2deg)
    newra = np.mod(newra, twopi)
    return newra*rad2deg, newdec*rad2deg
                

class target_info:
    def __init__(self):
        self.ticid = int(0)
        self.ra = 0.0
        self.dec = 0.0
        self.eclipLong = 0.0
        self.eclipLat = 0.0
        self.sectors = np.array([], dtype=np.int)
        self.onSiliconFlag = np.array([], dtype=np.int)
        self.possibleOnSiliconFlag = np.array([], dtype=np.int)
        self.cameras = np.array([], dtype=np.int)
        self.xpxs = np.array([], dtype=np.float)
        self.ypxs = np.array([], dtype=np.float)

def make_target_objects(tic, ra, dec):
    starList = []
    for i, curTic in enumerate(tic):
        curRa = ra[i]
        curDec = dec[i]
        # instantiate target object
        curTarg = target_info()
        curTarg.ra = curRa
        curTarg.dec = curDec
        curTarg.ticid = tic[i]
        # Convert ra and dec coords to ecliptic
        planCoords = SkyCoord(curRa, curDec, unit='deg')
        planEclipCoords = planCoords.transform_to(frame='barycentrictrueecliptic')
        curTarg.eclipLat = planEclipCoords.lat.deg
        curTarg.eclipLong = planEclipCoords.lon.deg
        starList.append(curTarg)
    return starList

def doRoughPosition(targinfo, scinfo):
    # Return the combinations of sector and detectors that can possibly observe
    #  target
    # go through each position in the spacecraft info class
    tRa = targinfo.ra
    tDec = targinfo.dec
    FOVDeg = 17.68
    nPoints = len(scinfo.sectors)
    targCoords = SkyCoord(tRa, tDec, unit='deg', frame='icrs')
    for iPnt in range(nPoints):
        for iCam in range(4):
            camRa = scinfo.camRa[iCam,iPnt]
            camDec = scinfo.camDec[iCam,iPnt]
            camCenter = SkyCoord(camRa, np.max([-89.99999,camDec]), unit='deg', frame='icrs')
            posAngs = camCenter.position_angle(targCoords)
            seps = camCenter.separation(targCoords)
            # Check for potentially on silicon
            if seps.deg < FOVDeg:
                # Append potential pointing camera combo to targets list
                targinfo.sectors = np.append(targinfo.sectors, scinfo.sectors[iPnt])
                targinfo.onSiliconFlag = np.append(targinfo.onSiliconFlag, 0)
                targinfo.possibleOnSiliconFlag = np.append(targinfo.possibleOnSiliconFlag, 1)
                targinfo.cameras = np.append(targinfo.cameras, iCam+1)
                targinfo.xpxs = np.append(targinfo.xpxs, 0.0)
                targinfo.ypxs = np.append(targinfo.ypxs, 0.0)
    return targinfo

## [Mast Query]
def mastQuery(request,proxy_uri=None):

    host='mast.stsci.edu'
    # Grab Python Version 
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent":"python-requests/"+version}

    # Encoding the request as a json string
    requestString = json.dumps(request)
    requestString = urlencode(requestString)
    
    # opening the https connection
    if None == proxy_uri:
        conn = httplib.HTTPSConnection(host)
    else:
        port = 443
        url = urlparse(proxy_uri)
        conn = httplib.HTTPSConnection(url.hostname,url.port)

        if url.username and url.password:
            auth = '%s:%s' % (url.username, url.password)
            headers['Proxy-Authorization'] = 'Basic ' + str(base64.b64encode(auth.encode())).replace("b'", "").replace("'", "")
        conn.set_tunnel(host, port, headers)

    # Making the query
    conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)

    # Getting the response
    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode('utf-8')

    # Close the https connection
    conn.close()

    return head,content
## [Mast Query]

def fileOutputHeader(fp, fpgParmFileList=None):
    # output a header to the file
    fp.write('# stars2px.py - Convert Target RA and Dec to TESS spacecraft pixel coordinates\n')
    fp.write('#   Original starspx.c by Alan Levine (MIT Kavli Institute)\n')
    fp.write('#   Python translation by Christopher Burke (MIT Kavli Institute)\n')
    fp.write('# Output columns Pipe Delimited; 16 header lines\n')
    fp.write('# File Creation: {:}\n'.format(datetime.datetime.now()))
    if fpgParmFileList is None:
        fp.write('# FPG Model Default\n')
    else:
        fp.write('# FPG Model {:s} {:s} {:s} {:s}\n'.format(fpgParmFileList[0], \
                 fpgParmFileList[1],fpgParmFileList[2],fpgParmFileList[3]))
    fp.write('# 1 [int] Input TIC ID\n')
    fp.write('# 2 [degree] Input or MAST TIC query target RA\n')
    fp.write('# 3 [degree] Input or MAST TIC query target Dec\n')
    fp.write('# 4 [degree] Ecliptic Longitude\n')
    fp.write('# 5 [degree] Ecliptic Latitude\n')
    fp.write('# 6 [int] - Observing Sector number for target\n')
    fp.write('# 7 [int] - Camera number for target\n')
    fp.write('# 8 [int] - Detector number for target\n')
    fp.write('# 9 [float] - Column pixel location for target\n')
    fp.write('# 10 [float] - Row pixel location for target\n')
             

def tess_stars2px_function_entry(starIDs, starRas, starDecs, trySector=None, scInfo=None, \
                              fpgParmFileList=None, combinedFits=False,\
                              noCollateral=False, aberrate=False):
    if scInfo == None:
        # Instantiate Spacecraft position info
        scinfo = TESS_Spacecraft_Pointing_Data(trySector=trySector, fpgParmFileList=fpgParmFileList)
    else:
        scinfo = scInfo
    # Now make list of the star objects
    starList = make_target_objects(np.atleast_1d(starIDs), \
                                   np.atleast_1d(starRas), np.atleast_1d(starDecs))
    # Make rough determination as to which pointing camera combos are worth
    # Checking in detail and then do detailed checking
    findAny=False
    outID = np.array([-1], dtype=np.int64)
    outEclipLong = np.array([-1.0], dtype=np.float)
    outEclipLat = np.array([-1.0], dtype=np.float)
    outSec = np.array([-1], dtype=np.int)
    outCam = np.array([-1], dtype=np.int)
    outCcd = np.array([-1], dtype=np.int)
    outColPix = np.array([-1.0], dtype=np.float)
    outRowPix = np.array([-1.0], dtype=np.float)
    for i, curTarg in enumerate(starList):
        starRas = np.array([curTarg.ra])
        starDecs =  np.array([curTarg.dec])
        for curSec in scinfo.sectors:
            idxSec = np.where(scinfo.sectors == curSec)[0][0]
            # Apply an approximate aberration correction
            if aberrate:
                useTime = Time(scinfo.midtimes[idxSec], format='jd')
                # Make coordinate object in ICRS coordinates
                ccat = SkyCoord(ra=starRas * u.deg,
                                 dec=starDecs * u.deg,
                                 obstime=useTime, frame='icrs')
                # convert to Geocentric aberrated coordinates
                # This is only an approximation to TESS
                #  because TESS orbits Earth and has 
                #  velocity <=4km/s relative to Earth whereas Earth is 30km/s
                cgcrs = ccat.transform_to('gcrs')
                starRas = np.array(cgcrs.ra.degree)
                starDecs = np.array(cgcrs.dec.degree)

            starInCam, starCcdNum, starFitsXs, starFitsYs, starCcdXs, starCcdYs = scinfo.fpgObjs[idxSec].radec2pix(\
                       starRas, starDecs)
            for jj, cam in enumerate(starInCam):
                # SPOC calibrated FFIs have 44 collateral pixels in x and are 1 based  
                xUse = starCcdXs[jj] + 45.0
                yUse = starCcdYs[jj] + 1.0
                xMin = 44.0
                ymaxCoord = 2049
                xmaxCoord = 2093
                if combinedFits:
                    xUse = starFitsXs[jj]
                    yUse = starFitsYs[jj]
                    xmaxCoord = 4097
                    ymaxCoord = 4097
                    xMin = 0.0
                if noCollateral:
                    xUse = starCcdXs[jj]
                    yUse = starCcdYs[jj]
                    xMin = 0.0
                if xUse>xMin and yUse>0 and xUse<xmaxCoord and yUse<ymaxCoord:
                    if findAny==False:
                        outID[0] = int(curTarg.ticid)
                        outEclipLong[0] = curTarg.eclipLong
                        outEclipLat[0] = curTarg.eclipLat
                        outSec[0] = curSec
                        outCam[0] = starInCam[jj]
                        outCcd[0] = starCcdNum[jj]
                        outColPix[0] = xUse
                        outRowPix[0] = yUse
                        findAny=True
                    else:
                        outID = np.append(outID, int(curTarg.ticid))
                        outEclipLong = np.append(outEclipLong, curTarg.eclipLong)
                        outEclipLat = np.append(outEclipLat, curTarg.eclipLat)
                        outSec = np.append(outSec, curSec)
                        outCam = np.append(outCam, starInCam[jj])
                        outCcd = np.append(outCcd, starCcdNum[jj])
                        outColPix = np.append(outColPix, xUse)
                        outRowPix = np.append(outRowPix, yUse)
    return outID, outEclipLong, outEclipLat, outSec, outCam, outCcd, \
            outColPix, outRowPix, scinfo

def tess_stars2px_reverse_function_entry(trySector, camera, ccd, colWant, rowWant, scInfo=None, \
                              fpgParmFileList=None):
    if scInfo == None:
        # Instantiate Spacecraft position info
        scinfo = TESS_Spacecraft_Pointing_Data(trySector=trySector, fpgParmFileList=fpgParmFileList)
    else:
        scinfo = scInfo
        
    idxSec = np.where(scinfo.sectors == trySector)[0][0]
    starRa, starDec = scinfo.fpgObjs[idxSec].pix2radec_nocheck_single(\
                       camera-1, ccd-1, [colWant-45.0, rowWant-1.0])
    return starRa, starDec, scinfo
    
if __name__ == '__main__':
    # Parse the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--ticId", type=int, \
                        help="TIC Id [int] for MAST coordinate query.  MUST Be Online For this option to Work!")
    parser.add_argument("-c", "--coord", type=float, nargs=2, \
                        help="RA and Dec of target [deg]")
    parser.add_argument("-f", "--inputFile", type=argparse.FileType('r'), \
                        help="Filename for input Target TIC [int]; RA[deg]; Dec[dec]; in white space delimited text file Column 1, 2, and 3 respectively")
    parser.add_argument("-o", "--outputFile", type=argparse.FileType('w'), \
                        help="Optional filename for output.  Default is output to stdout ")
    parser.add_argument("-s", "--sector", type=int, choices=range(1,56),\
                        help="Search a single sector Number [int]")
    parser.add_argument("-x", "--combinedFits", action='store_true', \
                        help="Output detector pixel coordinates for the 'Big' multi-detector combined fits file format")
    parser.add_argument("-xin","--noCollateral", action='store_true', \
                        help="Output detector pixel coordinates for an internal format where there are not leading collateral pixels and zero based")
    parser.add_argument("-fpg", "--fpgParameterFiles", nargs=4,\
                        help="Instead of default focal plane geometry parameters, list the 4 filenames for the fpg files to use.  Expects files in Al's format in camera numerical order")
    parser.add_argument("-r", "--reverse", nargs=5,\
                        help="Do reverse search.  Return RA Dec for a given pixel position.  5 parameters sector cam ccd colpix rowpix")
    parser.add_argument("-n", "--name", nargs=1, type=str, \
                        help="Search for a target by resolving its name with SESAME")
    parser.add_argument("-p", "--proxy_uri", nargs=1, type=str, \
                        help="Use proxy e.g.\"http://<user>:<passwd>@<proxy_server>:<proxy_port>\" ")
    parser.add_argument("-a", "--aberrate", action='store_true',\
                        help="Apply approximate aberration correction to input coordinates. Uses astropy GCRS coordinate frame to approximate TESS aberration")
    args = parser.parse_args()

# DEBUG BLOCK for hard coding input parameters and testing
#    class test_arg:
#        def __init__(self):
#            #self.ticId = 281541555
#            self.ticId = None
#            self.coord = None
#            self.name = ['KIC 6443093']
#            self.coord = [330.6803807390524, 42.27777178]
#            self.inputFile = None
#            self.sector = None
#            self.fpgParameterFiles = None
#            self.outputFile = None
#            self.combinedFits = False
#            self.noCollateral = False
#            self.reverse = None
#            self.reverse = [2,1,2,2092.0,1.0]
#    args = test_arg()
    
    # At least one Mode -t -c -f -r -n must have been specified
    if (args.ticId is None) and (args.coord is None) and (args.inputFile is None) and (args.reverse is None) and (args.name is None):
        print('You must specify one and only one mode -t, -c, -f, -r, -n')
        print('`python stars2px.py -h\' for help')
        sys.exit(1)

    # Check for reverse mode
    if (args.reverse is None):        
        # Do single coords first
        if args.coord is not None and args.name is None:
            nTarg = 1
            starTics = np.array([0], dtype=np.int64)
            starRas = np.array([args.coord[0]], dtype=np.float)
            starDecs = np.array([args.coord[1]], dtype=np.float)
        elif args.coord is None and args.name is not None:
            nTarg = 1
            starTics = np.array([0], dtype=np.int64)

            # Name resolve  in try except  for detecting problem
            try:
                coordinate = SkyCoord.from_name(args.name[0])
                print("Coordinates for {0}: ({1}, {2})"
                  .format(args.name[0], coordinate.ra.degree,
                          coordinate.dec.degree))
                starRas = np.array([coordinate.ra.degree], dtype=np.float)
                starDecs = np.array([coordinate.dec.degree], dtype=np.float)
            except:
                print("Could not resolve: {0}".format(args.name[0]))
                sys.exit(1)
        else:
            if not (args.inputFile is None): # Check for input file list next
                # Read in star positions in input
                # Now go through stars
                starFile = args.inputFile
                dataBlock = np.genfromtxt(starFile, dtype=['i4','f8','f8'])
                starTics = np.atleast_1d(dataBlock['f0'])
                starRas = np.atleast_1d(dataBlock['f1'])
                starDecs = np.atleast_1d(dataBlock['f2'])
            else:
                # Must have requested MAST query with TIC ID
                # Make a list of TICs using strings
                starTics = np.array([args.ticId], dtype=np.int64)
                ticStringList = ['{0:d}'.format(x) for x in starTics]    
                # Setup mast query
                request = {'service':'Mast.Catalogs.Filtered.Tic', \
                   'params':{'columns':'*', 'filters':[{ \
                            'paramName':'ID', 'values':ticStringList}]}, \
                    'format':'json', 'removenullcolumns':True}
                if args.proxy_uri is None:
                    headers, outString = mastQuery(request)
                else:
                    headers, outString = mastQuery(request,args.proxy_uri[0])
                outObject = json.loads(outString)
                starRas = np.array([x['ra'] for x in outObject['data']])
                starDecs = np.array([x['dec'] for x in outObject['data']])
        
        trySector = None
        if not (args.sector is None):
            trySector = args.sector
        fpgParmFileList = None
        if not (args.fpgParameterFiles is None):
            fpgParmFileList = [x for x in args.fpgParameterFiles]
    
            
        # Instantiate Spacecraft position info
        scinfo = TESS_Spacecraft_Pointing_Data(trySector=trySector, fpgParmFileList=fpgParmFileList)
        # Open output file if requested
    #    if not (args.outputFile is None):
    #        fout = open(args.outputFile, 'w')
    
        # Add header to outputfile
        if not (args.outputFile is None):
            fileOutputHeader(args.outputFile, fpgParmFileList=fpgParmFileList)
        else:
            # add single line header to stdout
            print('# TIC     |   RA      |   Dec     | EclipticLong | EclipticLat | Sector | Camera | Ccd | ColPix | RowPix | EdgeWarn')
        # Now make list of the star objects
        starList = make_target_objects(starTics, starRas, starDecs)
        #print('Finished converting coords to ecliptic')
        # Make rough determination as to which pointing camera combos are worth
        # Checking in detail and then do detailed checking
        findAny=False
        for i, curTarg in enumerate(starList):
            starRas = np.array([curTarg.ra])
            starDecs =  np.array([curTarg.dec])
            for curSec in scinfo.sectors:
                idxSec = np.where(scinfo.sectors == curSec)[0][0]
                # Apply an approximate aberration correction
                if args.aberrate:
                    useTime = Time(scinfo.midtimes[idxSec], format='jd')
                    # Make coordinate object in ICRS coordinates
                    ccat = SkyCoord(ra=starRas * u.deg,
                                     dec=starDecs * u.deg,
                                     obstime=useTime, frame='icrs')
                    # convert to Geocentric aberrated coordinates
                    # This is only an approximation to TESS
                    #  because TESS orbits Earth and has 
                    #  velocity <=4km/s relative to Earth whereas Earth is 30km/s
                    cgcrs = ccat.transform_to('gcrs')
                    starRas = np.array(cgcrs.ra.degree)
                    starDecs = np.array(cgcrs.dec.degree)
                starInCam, starCcdNum, starFitsXs, starFitsYs, starCcdXs, starCcdYs = scinfo.fpgObjs[idxSec].radec2pix(\
                           starRas, starDecs)
                for jj, cam in enumerate(starInCam):
                    # SPOC calibrated FFIs have 44 collateral pixels in x and are 1 based  
                    xUse = starCcdXs[jj] + 45.0
                    yUse = starCcdYs[jj] + 1.0
                    xMin = 44.0
                    ymaxCoord = 2049
                    xmaxCoord = 2093
                    if args.combinedFits:
                        xUse = starFitsXs[jj]
                        yUse = starFitsYs[jj]
                        xmaxCoord = 4097
                        ymaxCoord = 4097
                        xMin = 0.0
                    if args.noCollateral:
                        xUse = starCcdXs[jj]
                        yUse = starCcdYs[jj]
                        xMin = 0.0
                    if xUse>xMin and yUse>0 and xUse<xmaxCoord and yUse<ymaxCoord:
                        findAny=True
                        edgeWarn = 0
                        edgePix = 6
                        if xUse<=(xMin+edgePix):
                            edgeWarn = 1
                        if xUse>=(xmaxCoord-edgePix):
                            edgeWarn = 1
                        if yUse<=edgePix:
                            edgeWarn = 1
                        if yUse>=(ymaxCoord-edgePix):
                            edgeWarn = 1
                        strout = '{:09d} | {:10.6f} | {:10.6f} | {:10.6f} | {:10.6f} | {:2d} | {:1d} | {:1d} | {:11.6f} | {:11.6f} | {:1d}'.format(\
                           curTarg.ticid, curTarg.ra, curTarg.dec, curTarg.eclipLong,\
                           curTarg.eclipLat, curSec, starInCam[jj], starCcdNum[jj], xUse, yUse, edgeWarn)
                        if not (args.outputFile is None):
                            args.outputFile.write('{:s}\n'.format(strout))
                        else:
                            print(strout)
    
        if not findAny:
            print('No Target/s were found to be on detectors')
    # Do reverse mode
    else:
        trySector = int(args.reverse[0])
        fpgParmFileList = None
        if not (args.fpgParameterFiles is None):
            fpgParmFileList = [x for x in args.fpgParameterFiles]
        # Instantiate Spacecraft position info
        scinfo = TESS_Spacecraft_Pointing_Data(trySector=trySector, fpgParmFileList=fpgParmFileList)
        # Camera center ra and dec
        iCam = int(args.reverse[1])-1
        iCcd = int(args.reverse[2])-1
        colWnt = float(args.reverse[3])
        rowWnt = float(args.reverse[4])

        # Lets try going direct route with inversion of codes
        starCcdXs = colWnt - 45.0
        starCcdYs = rowWnt - 1.0
        idxSec = np.where(scinfo.sectors == trySector)[0][0]
        ra_deg, dec_deg = scinfo.fpgObjs[idxSec].pix2radec_nocheck_single(iCam, iCcd, [starCcdXs, starCcdYs])
        print(ra_deg, dec_deg)
        
#  OLD Way with minimizer       
#        def minFunc(x, iCam, iCcd, colWnt, rowWnt):
#            starCcdNum, starFitsXs, starFitsYs, starCcdXs, starCcdYs, lat = scinfo.fpgObjs[0].radec2pix_nocheck_single(\
#                               x[0], x[1], iCam, iCcd)
#            xUse = starCcdXs + 45.0
#            yUse = starCcdYs + 1.0
#            # Penalize latitudes <70
#            latBad = 0.0
#            if lat < 70.0:
#                latBad = 70.0-lat
#            zUse = np.power(xUse-colWnt,2) + np.power(yUse-rowWnt,2) + np.power(latBad,4) + 1.0
#            return zUse
#        # Use an initial minimize with bounds to keep things from going haywire
#        # start with camera center
#        optResult = opt.minimize(minFunc, [camRa, camDec], \
#                                 args=(iCam, iCcd, colWnt, rowWnt), method='TNC', \
#                                 bounds=[[0.0,410.0],[-90.0, 90.0]], tol=1.0e-6, \
#                                 options={'maxiter':500})
#        newRa = optResult.x[0]
#        newDec = optResult.x[1]
#        #  Refine minimization to hopefully converge
#        optResult2 = opt.minimize(minFunc, [newRa, newDec], \
#                                 args=(iCam, iCcd, colWnt, rowWnt), method='Nelder-Mead', \
#                                 tol=1.0e-6, \
#                                 options={'maxiter':500})
#        newRa2 = optResult2.x[0]
#        newDec2 = optResult2.x[1]
#        minQuality = optResult2.fun
#        newRa2 = np.mod(newRa2, 360.0)
#        if np.abs(minQuality - 1.0) > 1.0e-4:
#
#            newRa = camRa + (newRa-camRa)/2.0
#            newDec = camDec + (newDec-camDec)/2.0
#            optResult3 = opt.minimize(minFunc, [newRa, newDec], \
#                                 args=(iCam, iCcd, colWnt, rowWnt), method='Nelder-Mead', \
#                                 tol=1.0e-6, \
#                                 options={'maxiter':500}) 
#            newRa3 = optResult3.x[0]
#            newDec3 = optResult3.x[1]
#            minQuality = optResult3.fun
#            newRa3 = np.mod(newRa3, 360.0)
#            if newDec3 < -90.0:
#                newDec3 = -90.0 + (-90.0 - newDec3)
#                newRa3 = newRa3 + 180.0
#                newRa3 = np.mod(newRa3, 360.0)
#            if newDec3 > 90.0:
#                newDec3 = 90.0 + (90.0 - newDec3)
#                newRa3 = newRa3 + 180.0
#                newRa3 = np.mod(newRa3, 360.0)
#                
#            print(newRa3, newDec3, minQuality)
#            
#        else:
#            if newDec2 < -90.0:
#                newDec2 = -90.0 + (-90.0 - newDec2)
#                newRa2 = newRa2 + 180.0
#                newRa2 = np.mod(newRa2, 360.0)
#            if newDec2 > 90.0:
#                newDec2 = 90.0 + (90.0 - newDec2)
#                newRa2 = newRa2 + 180.0
#                newRa2 = np.mod(newRa2, 360.0)
#
#            print(newRa2, newDec2, minQuality)
        
