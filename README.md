# TESS-Point
High Precision TESS pointing tool.

Convert target coordinates given in Right Ascension and Declination to TESS detector pixel coordinates for the TESS prime mission 26 observing sectors (Year 1 & 2) and Year 3-4 up to Sectors 55.  Can also query MAST to obtain detector pixel coordinates for a star by TIC ID only (must be online for this option).  Provides the target ecliptic coordinates, Sector number, camera number, detector number, and pixel column and row.  If there is no output, then the target is not visible to TESS.

### Install or Upgrade
`pip install tess-point`

`pip install tess-point --upgrade`

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
Original programming in C and focal plane geometry solutions by Alan Levine (MIT).  This python translation by Christopher J. Burke (MIT).  Testing and focal plane geometry refinements by Michael Fausnaugh & Roland Vanderspek (MIT).  Testing by Thomas Barclay (NASA Goddard) and Jessica Roberts (Univ. of Colorado).  By target name resolving implemented by Brett Morris (UW).  Python help from Brigitta Sipocz and Martin Owens.  Bug reports by Adina Feinstein (Univ. Chicago). Proxy implementation by Dishendra Mishra.

### VERSION: 0.6.2

### WHAT'S NEW:
- Sector 46 field update 2021 October
- Too close to edge Warning flag now output in column. If a target is within 6 pixels of the edge of the science region (edgeWarn==1), then the target is unlikely to be assigned a 2minute or 20s aperture. The science pixels range in column from 45-2092 and row from 1-2048
- Year 4 pointings for Sectors 40-55 now available
- An approximate aberration correction is available with command line option.  Uses astropy GCRS Earth based frame which is close to TESS aberration

- Inverse transform (input Sector, Camera, CCD, pixel Column, pixel Row --> RA and Dec) is now 'analytic' rather than through brute force minimization.  The inverse transform is much faster and much more reliable.

### CITATION:
A citation for tess-point is available through the [Astrophysics Source Code Library](http://www.ascl.net/2003.001) entry. More complete BibTeX at bottom of page.

Burke, C. J., Levine, A., Fausnaugh, M., Vanderspek, R., Barclay, T., Libby-Roberts, J. E., Morris, B., Sipocz, B., Owens, M., Feinstein, A. D., Camacho, J., 2020, 0.4.1, Astrophysics Source Code Library, record ascl:2003:001

### NOTES
- Pointing table is for TESS Year 1 - 4(Sectors 1-55) .

- Testing shows pointing with this tool should be accurate to better than a pixel, but without including aberration effects, ones algorithm adopted for centroiding highly assymmetric point-spread function at edge of camera, and by-eye source location, a 2 pixel accuracy estimate is warranted. Use aberration option for better accuracy

- The output pixel coordinates assume the ds9 convention with 1,1 being the middle of the lower left corner.

- Pointing table is unofficial, and the pointings may change.

- See https://tess.mit.edu/observations/ for latest TESS pointing table

- No corrections for velocity aberration are calculated by default. Potentially more accurate results can be obtained if the target RA and Declination coordinates have aberration effects applied.  The aberrate option uses the astropy GCRS Earth based frame in order to approximate a TESS frame.  Earth has a velocity of 30km/s in solar system whereas TESS moves <4km/s relative to Earth, thus the GCRS correction should largely remove the 20arcsecond Earth induced aberration amplitude

- For proposals to the TESS science office or directors discretionary time, please consult the TESS prediction webtool available at https://heasarc.gsfc.nasa.gov/cgi-bin/tess/webtess/wtv.py for official identification of 'observable' targets.  However, if your proposal depends on a single or few targets, then this tool is helpful to further refine the likelihood of the target being available on the detectors.

- The calibrated FFI fits file release at MAST and calibrated by NASA Ames SPOC will have WCS information available to supplant this code.  The WCS generation is independent of the focal plane geometry model employed in this code and will give different results at the pixel level.  However, the WCS information is not available until the FFI files are released, whereas this code can predict positions in advance of data release.

- Hard coded focal plane geometry parameters from rfpg5_c1kb.txt

### OLD NOTES:
- Query by name using Sesame by Brett Morris

- Wrapper function implemented tess_stars2px_function_entry() with an example program, example_use_tess_stars2py_byfunction.py for using tess_stars2px in your own python program rather than on the command line.

- Pre filter step previously depended on the current mission profile of pointings aligned with ecliptic coordinates to work.  The pre filter step was rewritten in order to support mission planning not tied to ecliptic alignment.  End users should not see any change in results with this change.  However, local copies can be modified for arbitrary spacecraft ra,dec, roll and get same functionality.

- A reverse option is added to find the ra and dec for a given sector, camera, ccd, colpix, rowpix.  This is most useful for planning arbitrary pointing boundaries and internal use to identify targets on uncalibrated images that don't have WCS info available.  For precision work one should defer to WCS information on calibrated FFIs rather than this tool.


### TODOS:
1. Time dependent Focal plane geometry

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

### BibTeX

```
@MISC{2020ascl.soft03001B,
author = {{Burke}, C.~J. and {Levine}, A. and {Fausnaugh}, M. and {Vanderspek}, R. and {Barclay}, T. and {Libby-Roberts}, J.~E. and {Morris}, B. and {Sipocz}, B. and {Owens}, M. and {Feinstein}, A.~D. and {Camacho}, J.
},
title = "{TESS-Point: High precision TESS pointing tool}",
keywords = {Software },
howpublished = {Astrophysics Source Code Library},
year = 2020,
month = mar,
archivePrefix = "ascl",
eprint = {2003.001},
adsurl = {http://adsabs.harvard.edu/abs/2020ascl.soft03001B},
adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```