import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.table import Column
from astropy.time import Time
from datetime import datetime

from file_exception import MyException
from file_init import Mike
import warnings


class Val:
    def __init__(self, file):
        '''
        Initialize binary fits cube by reading its header and data from a given file.
        '''

        self.file = file

        pass


    def validate_primary_header(self):
        '''
        Validate the primary header of the given FITS file by:
            1) checking if it complies to 2880 byte header standard
            2) checking primary header cards
            3) checking for duplicate cards
            4) checking that each card is 80 character long
            5) recording missing values in the header
        '''
        
        self.validate_header_size()
        self.validate_primary_header_cards()
        self.validate_header_cards()

        self.file.validated_header = True

        return
        

    def validate_header_size(self):
        '''
        Validate that the header size is a multiple of 2880 bytes.
        '''

        header_size = len(self.file.header.tostring())
        if 0 != header_size % 2880:
            raise MyException("File header does not conform to 2880 byte standard!")
        
        return


    def validate_primary_header_cards(self):
        '''
        Validate primary header cards 'SIMPLE, BITPIX, NAXIS, and END' exist and are valid.
        '''

        if self.file.header.cards[0].keyword != 'SIMPLE':
            raise MyException("The first keyword in the header must be SIMPLE.")

        required_keys = ['SIMPLE', 'BITPIX', 'NAXIS']
        for key in required_keys:
            if key not in self.file.header:
                raise MyException(f"Required keyword '{key}' missing in header.")
            
        if self.file.header.get('BITPIX') not in [8, 16, 32, 64, -32, -64]:
            raise MyException("BITPIX has an invalid value.")
        
        naxis = self.file.header.get('NAXIS')
        for i in range(1, naxis + 1):
            axis_key = f"NAXIS{i}"
            if axis_key not in self.file.header:
                raise MyException(f"Missing {axis_key} keyword.")
            if self.file.header[axis_key] < 0:
                raise MyException(f"{axis_key} must be non-negative.")
            
        with open(self.file.file_path, 'rb') as f:
            block_size = 2880
            header_bytes = b""
            while True:
                chunk = f.read(block_size)
                if not chunk:
                    break
                header_bytes += chunk
                if b'END' in chunk:
                    break

            if b'END' not in header_bytes:
                raise MyException("FITS header does not contain the END keyword.")
            
        return


    def validate_header_cards(self):
        '''
        Validate all other header cards to ensure 80 byte length and there are no duplicates.
        '''

        seen = set()
        for card in self.file.header.cards:
            key, value, comment = card
            if key in seen and key != "COMMENT" and key != "HISTORY":
                raise MyException(f"Duplicate keyword found: {key}")
            seen.add(key)

        types = []
        for card in self.file.header.cards:
            key, value, comment = card
            types.append(type(value))

            card_str = card.__str__()
            card_length = len(card_str)

            if card_length != 80:
                raise MyException(f"Card '{key}' is {card_length} characters long: {card_str}")
            
            if value is None or (isinstance(value, str) and not value.strip()):
                self.file.missing_values.append(key)

        if len(self.file.missing_values) > 0:
            warnings.warn(f"Values do not exist for {self.file.missing_values}", stacklevel=2)

        return


    def validate_data(self):
        '''
        Validate data of the given FITS file by:
            1) checking for negative or NaN values
            2) validate fields of data by provided data type
            3) convert datetime fields to datetime objects
        '''

        # try:
        #     with fits.open(self.file.file_path) as hdul:
        #         self.file.dataH = hdul[1].header
        #         self.file.data = Table(hdul[1].data)

        # except Exception as e:
        #     raise MyException(f"Error reading data from FITS file: {e}")

        self.check_values()

        for column in self.file.data.colnames:
            # Check that column conforms to correct type
            self.match_types(column)

            # Convert the time columns to datetime
            self.convert_to_datetime(column)
            
            # Check if the numeric columns have weird numbers
            self.check_numbers(column)

        self.file.validated_data = True

        return


    def check_values(self):
        # Check if the data has bad values (NaN or negative)
        Data = np.ma.array(self.file.data["DATA"])

        # Check each array in Data for NaN or 0s and create a masked array
        mask = np.isnan(Data) | (Data < 0)
        if np.any(mask):
            warnings.warn("🚫 Data contains negative or NaN values.", stacklevel=2)

        self.file.data["DATA"] = np.ma.masked_array(Data, mask=mask)

        return


    def match_types(self, column):
        '''
        Make sure all the column types match the data types of the values.
        '''

        column_data = self.file.data[column]

        # Find how many unique types there are in the column
        unique_types = set(type(value) for value in self.file.data[column])
        
        # If there are more than one type of value in a column then correct it
        if len(unique_types) > 1:
            # Try to convert the values to the most common type
            common_type = max(unique_types, key=lambda t: sum(isinstance(value, t) for value in self.file.data[column]))
            for i, value in enumerate(self.file.data[column]):
                if not isinstance(value, common_type):
                    try:
                        self.file.data[column][i] = common_type(value)
                    except Exception:
                        # Mask the messed up values if conversion fails
                        if not hasattr(self.data[column], 'mask'):
                            self.data[column] = np.ma.array(self.data[column])
                        self.data[column].mask[i] = True
            print(f"🚫 Column '{column}' contains mixed data types: {unique_types}.")
        else:
            # Otherwise convert the column to the common type
            common_type = list(unique_types)[0]
            try:
                if common_type is str:
                    # Automatic conversion to str doesn't work some of the time so we need to try UTF-8 first
                    if column_data.dtype.char == 'S':  # Check if it's a byte string
                        self.data[column] = np.char.decode(self.data[column], encoding='utf-8', errors='replace')
                    else:
                        self.file.data[column] = np.array(self.file.data[column], dtype=str)

                # Make sure the arrays in the column are all floats
                elif common_type is np.ndarray:
                    if not all(isinstance(float(item), float) for subarray in column_data for item in subarray):
                        # If not try to make them floats, otherwise raise an error
                        for subarray in column_data:
                            for ind, item in enumerate(subarray):
                                try:
                                    float(item)
                                except Exception:
                                    raise TypeError(f"🚫 Column '{column}' contains arrays with non-float elements: {item}")

                # Otherwise just convert the column to the common type
                else:
                    self.file.data[column] = self.file.data[column].astype(common_type)

            except Exception as e:
                # If that all fails, raise an error
                MyException(f"⚠️ Failed to change column '{column}' data type to {common_type}: {e}")

        return
    

    def convert_to_datetime(self, column):
        '''
        Convert the dates and times to astropy time ojects from the astropy time library
        Ensure other time-step columns are floats
        '''
        #convert date and time to astropy time objects, or fall back to float or string
        if any(keyword in column.upper() for keyword in ["DATE", "TIME"]) or column.upper() == "LST":
                try:
                    # Try converting with astropy Time
                    self.data[column] = Time(self.data[column])
                except (ValueError, TypeError):
                    try:
                        # Fallback: convert to float (e.g., durations)
                        self.data[column] = Column(self.data[column].astype(float), dtype='f8')
                    except ValueError:
                        try:
                            # Fallback: treat as string
                            self.data[column] = Column(self.data[column].astype(str), dtype='U')
                        except Exception:
                            raise MyException(f"⚠️ Failed to convert value '{self.data[column]}' in column '{column}' to astropy Time, float, or string.")
                        
        # Convert other time step columns to floats, fallback to string
        if any(keyword in column.upper() for keyword in [ "DURATION", "EXPOSURE", "MJD", "UTC", "UTSECS"]) or column.upper() == "LST":
                try:
                    # Convert to float (e.g., durations)
                    self.data[column] = Column(self.data[column].astype(float), dtype='f8')
                except ValueError:
                    try:
                        # Fallback: treat as string
                        self.data[column] = Column(self.data[column].astype(str), dtype='U')
                    except Exception:
                        raise MyException(f"⚠️ Failed to convert value '{self.data[column]}' in column '{column}' to astropy Time, float, or string.")
                        
        return


    def check_numbers(self, column):
        '''
        Check if certain columns have negative values and remove them.
        '''
        if isinstance(self.data[column], Time):
            return
        
        if self.data[column].dtype.kind in {'U', 'S', 'O'}:
            # Vectorized replacement of 'nan' or 'NaN' (case-insensitive) with np.nan in string columns
            def replace_nan(val):
                if isinstance(val, str) and val.strip().lower() == 'nan':
                    return np.nan
                return val
            self.data[column] = np.vectorize(replace_nan)(self.data[column])

        if column.upper() in ["DURATION", "EXPOSURE", "TSYS", "TCAL", "LST", "ELEVATION", "TAMBIENT", "PRESSURE", "HUMIDITY", "RESTFREQ", "FREQRES", "TRGTLONG", "MJD", "UTSECS" ]:
            if np.any(self.file.data[column] < 0):
                num_negatives = np.sum(self.file.data[column] < 0)
                print(f"Found {num_negatives} negative values in column '{column}'. out of {len(self.file.data[column])} total values.")
                self.file.data = self.file.data[self.file.data[column] >= 0]
        
        return
    

if __name__ == "__main__":
    '''
    Test function to implement validation.
    '''

    file = Mike("EpicONoFFHiRes/0136550.fits")
    v = Val(file)
    v.validate_primary_header()
    v.validate_data()
    
    print(file.validated_header)
    print(file.validated_data)
