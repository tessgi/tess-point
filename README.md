# TESS-Point
High Precision TESS pointing tool.

Convert target coordinates given in Right Ascension and Declination to TESS detector pixel coordinates for the TESS prime mission 26 observing sectors (Year 1 & 2) focused on the southern and northern ecliptic planes.  Can also query MAST to obtain detector pixel coordinates for a star by TIC ID only (must be online for this option).  Provides the target ecliptic coordinates, Sector number, camera number, detector number, and pixel column and row.  If there is no output, then the target is not visible to TESS.

### Install
`pip install tess-point`

### Examples
- Display command line arguments and features

`python -m tess_stars2px -h`

- Return pixel coordinates for Pi Mensae ra and dec in degrees

`python -m tess_stars2px -c 84.291188 -80.469119`

- Return pixel coordinates for target by TIC ID (must be online)

`python -m tess_stars2px -t 261136679`

- Return pixel coordinates by target name (name resolved via [SESAME](http://cds.u-strasbg.fr/cgi-bin/Sesame))

`python -m tess_stars2px -n "pi Mensae"`

- Multi-target pixel coordinates results.  List the target TIC ID or other integer identifier [can be zero]; ra [deg]; dec [deg] in a whitespace delimited text file.  Process the target list.

`python -m tess_stars2px -f <target_list>`

Alternatively, the python module is a single file, tess_stars2px.py, so one can avoid pip install.  Just download tess_stars2px.py from github and use it in a local directory.  The above commands would then be python tess_starspx.py -t 261136679

- tess_stars2px can be called from a python program.  See example_use_tess_strs2py_byfunction.py for this way to use tess_stars2px

### AUTHORS
Original programming in C and focal plane geometry solutions by Alan Levine (MIT).  This python translation by Christopher J. Burke (MIT).  Testing and focal plane geometry refinements by Michael Fausnaugh & Roland Vanderspek (MIT).  Testing by Thomas Barclay (NASA Goddard) and Jessica Roberts (Univ. of Colorado).  By target name resolving implemented by Brett Morris (UW).  Python help from Brigitta Sipocz and Martin Owens.

### VERSION: 0.3.3

### WHAT'S NEW:
- TESS Year 2 (Sectors 14-26) potential pointings added. Assumes unshift Sector 16 pointing.

### NOTES
- Pointing table is for TESS Year 1 and 2 (Sectors 1-26).

- Pointing prediction algorithm is same as employed internally at MIT for target management.  However, hard coded focal plane geometry is not up to date and may contain inaccurate results.

- Testing shows pointing with this tool should be accurate to better than a pixel, but without including aberration effects, ones algorithm adopted for centroiding highly assymmetric point-spread function at edge of camera, and by-eye source location, a 2 pixel accuracy estimate is warranted.

- The output pixel coordinates assume the ds9 convention with 1,1 being the middle of the lower left corner.

- Pointing table is unofficial, and the pointings may change.

- See https://tess.mit.edu/observations/ for latest TESS pointing table

- No corrections for velocity aberration are calculated. Potentially more accurate results can be obtained if the target RA and Declination coordinates have aberration effects applied.

- For proposals to the TESS science office or directors discretionary time, please consult the TESS prediction webtool available at https://heasarc.gsfc.nasa.gov/cgi-bin/tess/webtess/wtv.py for official identification of 'observable' targets.  However, if your proposal depends on a single or few targets, then this tool is helpful to further refine the likelihood of the target being available on the detectors.

- The calibrated FFI fits file release at MAST and calibrated by NASA Ames SPOC will have WCS information available to supplant this code.  The WCS generation is independent of the focal plane geometry model employed in this code and will give different results at the pixel level.  However, the WCS information is not available until the FFI files are released, whereas this code can predict positions in advance of data release.

- Hard coded focal plane geometry parameters from rfpg5_c1kb.txt

### OLD NOTES:
- Query by name using Sesame by Brett Morris

- Wrapper function implemented tess_stars2px_function_entry() with an example program, example_use_tess_stars2py_byfunction.py for using tess_stars2px in your own python program rather than on the command line.

- Pre filter step previously depended on the current mission profile of pointings aligned with ecliptic coordinates to work.  The pre filter step was rewritten in order to support mission planning not tied to ecliptic alignment.  End users should not see any change in results with this change.  However, local copies can be modified for arbitrary spacecraft ra,dec, roll and get same functionality.

- A reverse option is added to find the ra and dec for a given sector, camera, ccd, colpix, rowpix.  This is most useful for planning arbitrary pointing boundaries and internal use to identify targets on uncalibrated images that don't have WCS info available.  For precision work one should defer to WCS information on calibrated FFIs rather than this tool.  The reverse is a brute force 'hack' that uses a minimizer on the forward direction code to find ra and dec.  In principle it is possible to reverse the matrix transforms to get the ra and dec directly, but I chose this less efficient method for expediency.  The minimizer is not guaranteed to converge at correct answer.  The current method is a slow way to do this.



### TODOS:
1. Check python 2.7 compatability
2. Include approximate or detailed velocity aberration corrections
3. Time dependent Focal plane geometry

### DEPENDENCIES:
- python 3+
- astropy
- numpy

### SPECIAL THANKS TO:
Includes code from the python MAST query examples 
https://mast.stsci.edu/api/v0/pyex.html

### IMPLEMENTATION DETAILS:
In summary, the code begins with a space craft bore site pointing in RA, Dec, and roll angle.  A series of Euler angle translation matrices are calculated based upon the space craft bore site.  Next the target coordinates in RA and Dec are translated to the space craft bore site frame.  Next, the target coordinates are translated to each of the four TESS camera frames.  Once target coordinates are translated to the  camera frame the radial position of the target relative to the camera center is checked to see if it is potentially in the camera field of view. If so, the focal plane position is calculated using a radial polynomial model with a constant term and terms the even powers (2nd ,4th , and 8th).  Rotations are applied to convert the on sky positions to the detector readout directions.

### Notes to self
1. Make code changes
2. Update version number in README.md, code, and setup.py
3. git add, commit, push
4. upload to PyPI - python setup.py sdist upload -r pypi
5. Make release on github
