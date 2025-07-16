# Third-party libraries
import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.table import Column
from astropy.time import Time


class Val:
    def __init__(self, filepath):
        self.filepath = filepath

        try:
            with fits.open(filepath) as hdul:
                self.header = hdul[0].header
                self.data = Table(hdul[1].data)
            # TODO add logging here
        except Exception as e:
            # TODO add logging here
            pass


    def validate_header_size(self):
        header_size = len(self.header.tostring())
        if 0 != header_size % 2880:
            # TODO add logging and error handling here
            pass
    

    def validate_primary_header_cards(self):
        if self.header.cards[0].keyword != "SIMPLE":
            # TODO add logging here
            pass

        required_keys = ['BITPIX', 'NAXIS']
        for key in required_keys:
            if key not in self.header:
                # TODO add logging and handling for missing keys here
                pass

        if self.header.get('BITPIX') not in [8, 16, 32, 64, -32, -64]:
            # TODO add handling here for endianness
            pass

        naxis = self.header.get('NAXIS')
        for i in range(1, naxis + 1):
            axis_key = f"NAXIS{i}"
            if axis_key not in self.header:
                # TODO add error handling here for missing axis
                pass
            if self.header[axis_key] < 0:
                # TODO add error handling here for negative axis keys
                pass

        with open(self.filepath, 'rb') as f:
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
                # TODO add error handling and logging for bad header here
                pass


    def validate_header_cards(self):
        seen = set()
        for card in self.header.cards:
            key, value, comment = card
            if key in seen and key not in ('COMMENT', 'HISTORY'):
                # TODO add logging for duplicate cards here
                pass
            
            seen.add(key)

        
        self.missing_values = []
        for card in self.header.cards:
            key, value, comment = card

            card_str = str(card)
            card_length = len(card_str)

            if card_length != 80:
                # TODO add logging for incorrect length cards
                pass

            if value is None or (isinstance(value, str) and not value.strip()):
                # TODO add logging for missing values for cards here
                pass

    
    def validate_primary_header(self):
        self.validate_header_size()
        self.validate_primary_header_cards()
        self.validate_header_cards()
        # File type validation no longer here --> moved to external pipeline workflow

        # TODO add logging here for completed header


    def check_values(self):
        data_array = np.ma.array(self.data['DATA'])

        mask = np.isnan(data_array)
        if np.any(mask):
            # TODO add logging to warn negative values
            pass

        self.data['DATA'] = np.ma.masked_array(data_array, mask = mask)


    def match_types(self, column):
        column_data_array = self.data[column]

        unique_types = set(type(value) for value in column_data_array)

        if len(unique_types) > 1:
            common_type = max(unique_types, key = lambda t: sum(isinstance(value, t) for value in column_data_array))

            for ind, i in enumerate(column_data_array):
                if not isinstance(i, common_type):
                    try:
                        column_data_array[ind] = common_type(i)
                    except:
                        if not hasattr(column_data_array, 'mask'):
                            column_data_array = np.ma.array(column_data_array)
                        column_data_array.mask[ind] = True

            # TODO add logging to state "column has mixed data types in X field"

        else:
            common_type = next(iter(unique_types))  # ✅ This works
            try:
                if common_type is str:
                    if column_data_array.dtype.char == 'S':
                        column_data_array = np.char.decode(column_data_array, encoding='utf-8', errors='replace')
                    else:
                        column_data_array = np.array(column_data_array[column], dtype = str)

                elif common_type is np.ndarry:
                    if not all(isinstance(float(item), float) for subarray in column_data_array for item in subarray):
                        for subarray in column_data_array:
                            for ind, i in enumerate(subarray):
                                try:
                                    float(i)
                                except:
                                    # TODO add logging to state that the array contains non float elements
                                    pass

                else:
                    self.data[column] = self.data[column].astype(common_type)

            except:
                # TODO add logging for failed conversions
                pass

            
    def convert_to_datetime(self, column):
        column_data_array = self.data[column]

        if any(keyword in column.upper() for keyword in ["TIME"]) or column.upper() == "LST":
            try:
                # Try converting with astropy Time
                type = column_data_array.dtype.type
                if type is str or type is np.str_:
                    for i, value in enumerate(column_data_array):
                        # Replace underscores with dashes and the first colon with 'T' for ISO format
                        new_value = value.replace('_', '-')

                        if 'T' not in new_value:
                            new_value = new_value.replace(':', 'T', 1)

                        self.data[column][i] = new_value
                self.data[column] = Time(column_data_array)

            except (ValueError, TypeError):
                try:
                    self.data[column] = Column(column_data_array.astype(float), dtype='f8')
                except ValueError:
                    try:
                        self.data[column] = Column(column_data_array.astype(str), dtype='U')
                    except:
                        # TODO add logging for all methods failed to convert to datetime
                        pass


        if any(keyword in column.upper() for keyword in ['DATE', 'TIME', 'DURATION', 'EXPOSURE', 'MJD', 'UTC', 'UTSECS']) or column.upper() == 'LST':
            try:
                self.data[column] = Column(column_data_array.astype(float), dtype='f8')
            except ValueError:
                try:
                    self.data[column] = Column(column_data_array.astype(str), dtype='U')
                except:
                    # TODO add logging for failed to convert times to float or str
                    pass


    def check_numbers(self, column):
        column_data_array = self.data[column]

        if isinstance(self.data[column], Time):
            return
        
        if column_data_array.dtype.kind in {'U', 'S', 'O'}:
            def replace_nan(val):
                if isinstance(val, str) and val.strip().lower() == 'nan':
                    return np.nan
                
                return val
            
            self.data[column] = np.vectorize(replace_nan)(column_data_array)

        if column.upper() in ["DURATION", "EXPOSURE", "TSYS", "TCAL", "LST", "ELEVATION", "TAMBIENT", "PRESSURE", "HUMIDITY", "RESTFREQ", "FREQRES", "TRGTLONG", "MJD", "UTSECS" ]:
            if np.any(column_data_array < 0):
                num_negatives = np.sum(column_data_array < 0)

                self.data = column_data_array[column_data_array >= 0]

                # TODO add logging for time array less than zero
        

    def validate_data(self):
        self.check_values()
        
        for column in self.data.colnames:
            self.match_types(column)
            self.convert_to_datetime(column)
            self.check_numbers(column)


    # def find_channels(self):
        

        
    def find_calibrations(self):
        ifnums = np.unique(self.data['IFNUM'])
        plnums = np.unique(self.data['PLNUM'])

        channel_count = ifnums * plnums

        data_start_ind = None
        post_cal_start_ind = None

        counter = 0

        cal_started = False
        pre_cal_complete = False

        for ind, i in enumerate(self.data):
            if i['CALSTATE'] == 1:
                cal_started = True
            
            if cal_started and i["CALSTATE"] == 0 and i["SWPVALID"] == 1 and not pre_cal_complete:
                data_start_ind = ind
                pre_cal_complete = True

            if pre_cal_complete and i["SWPVALID"] == 0 and self.data[ind - 1]["SWPVALID"] == 0:
                if post_cal_start_ind is None:
                    post_cal_start_ind = ind
            else:
                post_cal_start_ind = None

            if pre_cal_complete and i['CALSTATE'] == 0 and i['SWPVALID'] == 1:
                counter += 1

            if counter <= 3 * channel_count and i['SWPVALID'] == 1:
                pre_cal_complete = False

            if pre_cal_complete and i['SWPVALID'] == 0 and i['CALSTATE'] == 1:
                break

        indices = [data_start_ind, post_cal_start_ind]

        if self.header['OBSMODE'] == 'onoff':
            for ind, i in enumerate(self.data):
                target = 'onoff:off'

                if target in i['OBSMODE']:
                    offstart = ind
                    indices = [data_start_ind, offstart-1, offstart, post_cal_start_ind]

                    break

        
        index_str = ",".join(str(ind) if ind is not None else "None" for ind in indices)
        self.header['INDICES'] = (index_str, 'Indices: data_start, post_cal_start[, off_start-1, off_start]')

            
    # def get_frequency_range(self):


            

    
    # def sort_data(self):
    #     self.find_channels() --> needed in gain calibration for per channel gain deltas
    #     self.find_calibrations() --> still used here
    #     self.get_frequency_range() --> parses range from header but not needed here because it would just be rewritten to the header but we already parse it from there
    #     self.get_channel_range() --> same as above


    def save(self, output_path=None):
        if output_path is None:
            base, ext = os.path.splitext(self.filepath)
            output_path = f"{base}_validate{ext}"

        # Create Primary HDU with updated header
        primary_hdu = fits.PrimaryHDU(header=self.header)

        # Convert Table back to FITS HDU
        table_hdu = fits.BinTableHDU(data=self.data)

        # Combine HDUs
        hdulist = fits.HDUList([primary_hdu, table_hdu])

        # Write to file
        hdulist.writeto(output_path, overwrite=True)
        

    def validate(self):
        self.validate_primary_header()
        self.validate_data()

        self.save()

        # self.find_calibrations() --> this can be done as part of gain calibration

if __name__ == "__main__":
    filepath = "EpicONoFFHiRes/0135417_validate.fits"
    val = Val(filepath)
    val.validate_primary_header()