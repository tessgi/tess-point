from dataclasses import dataclass
from typing import Optional, List

@dataclass
class TESSpoint(object):    
    coord : Optional = None # SkyCoord?
    ra : Optional = None
    dec : Optional = None
    row : Optional = None
    column : Optional = None
    targetname : Optional = None
    filename : Optional[str] = None
    
    def __post_init__(self):
        if self.targetname is not None:
            self._read_name()
        elif self.filename is not None:
            self._read_csv()
        elif self.coord is not None:
            self._read_skycoord()
        self.validate()
    
    def _read_skycoord(self):
        if isinstance(self.coord, list):
            if isinstance(self.coord[0], SkyCoord):
                self.coord = SkyCoord(self.coord)
            else:
                raise ValueError("Must pass either a `astropy.coordinate.SkyCoord` object or a list of `astropy.coordinate.SkyCoord` objects to `coord`")
        elif not isinstance(self.coord, SkyCoord):
            raise ValueError("Must pass either a `astropy.coordinate.SkyCoord` object or a list of `astropy.coordinate.SkyCoord` objects to `coord`")
        if len(self.coord.shape) == 0:
            self.coord = SkyCoord([self.coord])
        self.ra, self.dec =  self.coord.ra.deg, self.coord.dec.deg
    
    def _read_name(self):
        if isinstance(self.targetname, str):
            c = SkyCoord.from_name(self.targetname)
            self.ra, self.dec = c.ra.deg, c.dec.deg
        elif isinstance(self.targetname, (list, np.ndarray)):
            self.ra, self.dec = [], []
            for name in self.targetname:
                c = SkyCoord.from_name(name)
                self.ra.append(c.ra.deg)
                self.dec.append(c.dec.deg)
    
    def _read_csv(self):
        df = pd.read_csv(self.filename)
        cols = np.asarray([c.lower().strip() for c in df.columns])
        if np.any(cols == 'col'):
            cols[cols == 'col'] = 'column'
        
        if not np.in1d(['ra', 'dec', 'row', 'column'], cols).any():
            raise ValueError('Must pass a dataframe with column headers of "ra", "dec", "column", or "row".')
        
        [setattr(self, attr, np.asarray(df[attr])) for attr in ['ra', 'dec', 'row', 'column'] if attr in cols]
        
        
    def validate(self):
        attrs = np.asarray(['ra', 'dec', 'row', 'column'])
        
        isnone = np.asarray([getattr(self, attr) is None for attr in attrs])
        # Passed in something
        if isnone.all():
            raise ValueError(f"Must pass either RA and Dec, Column and Row, a target name, or a filename.")

        if np.atleast_1d(np.where(~isnone)[0] == [0, 1]).all():
            self.coord_type = 'radec'
        elif np.atleast_1d(np.where(~isnone)[0] == [2, 3]).all():
            self.coord_type = 'pixels'
        else:
            raise ValueError("Must pass either RA and Dec, or Column and Row.")

        # Correct length
        valid_lengths = len(np.unique([len(np.atleast_1d(getattr(self, attr))) for attr in attrs[np.where(~isnone)]])) == 1
        if not valid_lengths:
            raise ValueError("Must pass arrays of the same length.")
        [setattr(self, attr, np.atleast_1d(getattr(self, attr))) for attr in attrs[np.where(~isnone)]]
        self.nvals = len(np.atleast_1d(getattr(self, attrs[np.where(~isnone)[0][0]])))
             
    def __len__(self):
        return self.nvals
    
    def __repr__(self):
        if self.coord_type == 'pixels':
            return f'TESSpoint ({self.nvals} Row/Column pairs)'
        elif self.coord_type == 'radec':
            return f'TESSpoint ({self.nvals} RA/Dec pairs)'
        else:
            return 'TESSpoint'
        
    def __getitem__(self, idx):
        if self.coord_type == 'radec':
            return TESSpoint(ra=self.ra[idx], dec=self.dec[idx])
        elif self.coord_type == 'pixels':
            return TESSpoint(row=self.row[idx], column=self.column[idx])
        else:
            raise ValueError('No `coord_type` set.')
        
            
    def to_RADec(self, sector:Optional[List] = None, camera:Optional[List] = None, ccd:Optional[List] = None):
        raise NotImplementedError
        
        #return df # with dates?
        
    def to_Pixel(self, sector:Optional[List] = None):
        raise NotImplementedError
        #return df # with dates?
        
    def ObservabilityMask(self, sectors:Optional[List] = None, cycle:Optional[List] = None):  
        raise NotImplementedError
        #return np.ndarray of bools
        
    def NumberOfObservations(self, sectors:Optional[List] = None, cycle:Optional[List] = None):  
        raise NotImplementedError
        #return np.ndarray of ints