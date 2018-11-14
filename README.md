# TESS-Point
High Precision TESS pointing tool.

Convert target coordinates given in Right Ascension and Declination to TESS detector pixel coordinates for the first 13 TESS observing sectors (Year 1) focused on the southern ecliptic plane.  Can also query MAST to obtain detector pixel coordinates for a star by TIC ID only (must be online for this option).  Provides the target ecliptic coordinates, Sector number, camera number, detector number, and pixel column and row.  If there is no output, then the target is not visible to TESS.

### Install
pip install tess-point


### Examples
- Display command line arguments and features

`python -m tess_stars2px -h`

- Return pixel coordinates for Pi Mensae ra and dec in degrees

`python -m tess_stars2px -c 84.291188 -80.469119`

- Return pixel coordinates for target by TIC ID (must be online)

`python -m tess_stars2px -t 261136679`

- Multi-target pixel coordinates results.  List the target TIC ID or other integer identifier [can be zero]; ra [deg]; dec [deg] in a whitespace delimited text file.  Process the target list.

`python -m tess_stars2px -f <target_list>`

Alternatively, the python module is a single file, tess_stars2px.py, so one can avoid pip install.  Just download tess_stars2px.py from github and use it in a local directory.  The above commands would then be python tess_starspx.py -t 261136679

### AUTHORS
Original programming in C and focal plane geometry solutions by Alan Levine (MIT).  This python translation by Christopher J. Burke (MIT).  Testing and focal plane geometry refinements by Michael Fausnaugh & Roland Vanderspek (MIT).  Testing by Thomas Barclay (NASA Goddard)

### VERSION: 0.1

### NOTES
- Pointing table is only for TESS Year 1 (Sectors 1-13) in Southern Ecliptic.

- Pointing prediction algorithm is same as employed internally at MIT for target management.  However, hard coded focal plane geometry is not up to date and may contain inaccurate results.

- Testing shows pointing with this tool should be accurate to better than a pixel, but without including aberration effects, ones algorithm adopted for centroiding highly assymmetric point-spread function at edge of camera, and by-eye source location, a 2 pixel accuracy estimate is warranted.

- The output pixel coordinates assume the ds9 convention with 1,1 being the middle of the lower left corner.

- Pointing table is unofficial, and the pointings may change.

- See https://tess.mit.edu/observations/ for latest TESS pointing table

- No corrections for velocity aberration are calculated. Potentially more accurate results can be obtained if the target RA and Declination coordinates have aberration effects applied.

- For proposals to the TESS science office or directors discretionary time, please consult the TESS prediction webtool available at https://heasarc.gsfc.nasa.gov/cgi-bin/tess/webtess/wtv.py for official identification of 'observable' targets.  However, if your proposal depends on a single or few targets, then this tool is helpful to further refine the likelihood of the target being available on the detectors.

- The calibrated FFI fits file release at MAST and calibrated by NASA Ames SPOC will have WCS information available to supplant this code.  The WCS generation is independent of the focal plane geometry model employed in this code and will give different results at the pixel level.  However, the WCS information is not available until the FFI files are released, whereas this code can predict positions in advance of data release.

- Hard coded focal plane geometry parameters from rfpg5_c1kb.txt

### TODOS:
1. Check python 2.7 compatability
1. Currently the python module expects to be run on command line provide a wrapper function such that the module can be more readily used with external python codes
1. Include approximate or detailed velocity aberration corrections
2. Provide estimated pointing table for TESS Year 2
3. Time dependent Focal plane geometry
4. Do the reverse transormation go from pixel to RA and Dec

### DEPENDENCIES:
- python 3+
- astropy
- numpy

### SPECIAL THANKS TO:
Includes code from the python MAST query examples 
https://mast.stsci.edu/api/v0/pyex.html

### IMPLEMENTATION DETAILS:
In summary, the code begins with a space craft bore site pointing in RA, Dec, and roll angle.  A series of Euler angle translation matrices are calculated based upon the space craft bore site.  Next the target coordinates in RA and Dec are translated to the space craft bore site frame.  Next, the target coordinates are translated to each of the four TESS camera frames.  Once target coordinates are translated to the  camera frame the radial position of the target relative to the camera center is checked to see if it is potentially in the camera field of view. If so, the focal plane position is calculated using a radial polynomial model with a constant term and terms the even powers (2nd ,4th , and 8th).  Rotations are applied to convert the on sky positions to the detector readout directions.
