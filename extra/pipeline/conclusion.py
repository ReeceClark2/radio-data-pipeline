# Standard library
import copy
from operator import ilshift, index

# Third-party libraries
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from matplotlib.pyplot import axis

# Local application imports
import os
from matplotlib.pylab import f
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.table import Column
from astropy.time import Time
from scipy.stats import linregress
import rcr
from collections import defaultdict
from validate import Val
import csv
from dep_gain_calibrate import Gain_Calibrate

class Conclusion:
    def __init__(self, filepath, index, axis, slice_type, feeds=[]):
        self.filepath = filepath
        try:
            with fits.open(filepath) as hdul:
                self.header = hdul[0].header
                self.data = Table(hdul[1].data)
        except Exception as e:
            print("Error opening FITS file:", e)
            pass

        

        # Copy other attributes as needed
        self.data_start_ind = None
        self.post_cal_start_ind = None 
        self.offstart = None
        self.freqs = []
        self.labels = []

        #TO be saved
        self.params = [index, axis, slice_type, feeds]
        self.continuum = []
        self.flux_calibrated = []
        self.spectrum = []
        
        self.get_startend_freqs()
        self.find_calibrations()  # Find calibration indices
        self.split_slp_feed()



        #fluxily calibration

    def find_calibrations(self):
        # Attempt to read data indices from header cards
        for key, value in self.header.items():
            if key == 'DATAIND':
                self.data_start_ind = value
            elif key == 'POSTCIND':
                self.post_cal_start_ind = value
            elif key == 'ONOFFIND':
                self.offstart = value

    def get_startend_freqs(self):
        '''
        Get the start and stop frequencies for each channel in the data.
        param file: Radio_File class file
        returns: populates the file's freqs field with start and stop frequencies
        '''
        # Scour the header for the bandwidth and center frequencies
        for key, value in self.header.items():
            # Bandwidth
            if key == ("OBSBW"):
                band = value

            # Lowres center frequency
            elif key == ("OBSFREQ"):
                center = [value]

            # HIRES bands center frequencies
            elif key == ("HISTORY"):
                # If HIRES BANDS exist replace center with the HIRES center frequencies
                if value.startswith("HIRES bands"):
                    # Extract all integers from the string
                    value = value.replace(",", " ").strip()
                    value = value.split(" ")
                    center = []

                    # Split the value into individual words and numbers
                    for k in value:
                        k = str(k).strip()
                        try:
                            # If it's a float add it to the center list
                            k = float(k)
                            center.append(k)
                        except ValueError:
                            # Otherwise skip it
                            continue

        channels = len(np.unique([d['PLNUM'] for d in self.data]))



        for c in center:
            start_f = c - (band / 2)
            stop_f = c + (band / 2)
            
            for i in range(channels):
                self.freqs.append(np.array([start_f, stop_f]))
                # TODO add logger
        # New system: start and stop frequencies stored in a single Card with key "Startstop"
        # for key, value in self.header.items():
        #     if "Startstop" in self.header:
        #         # Expecting a string like "start,stop" or a list/tuple
        #         startstop = self.header["Startstop"]
        #         if isinstance(startstop, str):
        #             try:
        #                 start_f, stop_f = map(float, startstop.split(","))
        #             except Exception:
        #                 raise ValueError("Header 'Startstop' value is not in 'start,stop' format.")
        #         elif isinstance(startstop, (list, tuple)) and len(startstop) == 2:
        #             start_f, stop_f = float(startstop[0]), float(startstop[1])
        #         else:
        #             raise ValueError("Header 'Startstop' value is not a valid format.")
        #             channels = len(np.unique([d['PLNUM'] for d in self.data]))
        #         for i in range(channels):
        #             self.freqs.append(np.array([start_f, stop_f]))

    def user_cuts(self, indices, axis, slice_type, feeds=[]):
        """        Apply user-defined cuts to the data based on frequency or time intervals.
        indices: list of tuples defining the intervals to keep or cut
        axis: "spectrum" for frequency cuts, "continuum" for time cuts
        slice_type: "keep" to keep the specified intervals, "cut" to remove them
        feeds: list of feed numbers to apply the cuts to; if empty, applies to all feeds
        """  
        if feeds == []:
            feeds = np.arange(len(np.unique(d['IFNUM'] for d in self.data)))
        # Ensure indices are in pairs


        if axis == "spectrum":
            for i, c in enumerate(self.data):
                feednum = np.unique(c['IFNUM'])
                if feednum.size != 1:
                    print("Data is not split by feed. Please run split_slp_feed() first.")
                feednum = feednum[0]
                if feednum not in feeds:
                    continue

                data_column = c['DATA']
                length = data_column.shape[1]
                freqs = np.linspace(self.freqs[i][0], self.freqs[i][1], length)

                freq_mask = np.zeros_like(freqs, dtype=bool)
                for fmin, fmax in indices:
                    if slice_type == "keep":
                        freq_mask |= (freqs >= fmin) & (freqs <= fmax)
                    elif slice_type == "cut":
                        freq_mask |= (freqs >= fmin) & (freqs <= fmax)
                if slice_type == "cut":
                    freq_mask = ~freq_mask

                # Apply mask to each row in the DATA column
                sliced_data = np.array([row[freq_mask] for row in data_column])
                selected_freqs = freqs[freq_mask]

                # Update the DATA column with the sliced data
                c.replace_column('DATA', sliced_data)
                self.spectrum.append(selected_freqs)


        if axis == "continuum":
            for i, c in enumerate(self.data):
                feednum = np.unique(c['IFNUM'])
                if feednum.size != 1:
                    print("Data is not split by feed. Please run split_slp_feed() first.")
                feednum = feednum[0]
                
                if feednum not in feeds:
                    continue

                times = Time(c["DATE-OBS"], format="isot")
                t0 = Time(self.header["DATE"], format="isot")
                time_rel = (times - t0).sec  # Time delta in seconds

                # Check if any indices are out of bounds (before data start or after post cal start)
                data_start = time_rel[self.data_start_ind//len(self.data)] if self.data_start_ind is not None else 0
                post_cal_end = time_rel[self.post_cal_start_ind//len(self.data)] if self.post_cal_start_ind is not None else len(time_rel) - 1
                for tmin, tmax in indices:
                    if tmin < data_start:
                        print(f"Warning: tmin={tmin} is before data start index ({data_start})")
                    if tmax > post_cal_end:
                        print(f"Warning: tmax={tmax} is after post calibration start index ({post_cal_end})")

                new_table= []
                if slice_type == "keep":
                    for tmin, tmax in indices:
                        for j in range(len(time_rel)):
                            if time_rel[j] > tmin and time_rel[j] < tmax:
                                new_table.append(c[j])
                elif slice_type == "cut":
                    # Start with all rows, then filter out those within any (tmin, tmax) interval
                    mask = np.ones(len(time_rel), dtype=bool)
                    for tmin, tmax in indices:
                        for j in range(len(time_rel)):
                            if tmin < time_rel[j] < tmax:
                                mask[j] = False
                    new_table = [c[j] for j in range(len(time_rel)) if mask[j]]
                    # Convert new_table to an astropy Table only if it's not empty

                self.data[i] = Table(rows=new_table, names=c.colnames)
    
    
    def split_slp_feed(self):
        '''
        Split the data of the provided file by channel and feed.
        '''

        ifnums = np.unique(self.data["IFNUM"])
        plnums = np.unique(self.data["PLNUM"])

        data = []
        labels = []
        for i in ifnums:
            for j in plnums:
                subset_mask = (self.data["IFNUM"] == i) & (self.data["PLNUM"] == j)
                subset_data = self.data[subset_mask]
                data.append(subset_data)
                labels.append(f'Feed{i + 1},Channel{j + 1}')
        self.data = data
        self.labels = labels

        return

    
    def get_conversion_factor(self):
        """
        Calculate the conversion factor for flux calibration.
        This method should be implemented based on the specific requirements of the flux calibration process.
        """
        # Placeholder for conversion factor logic
        # This should include the actual implementation of conversion factor calculation
        with open("radio-data-pipeline/calibration_data/cal_obs_ids.csv", newline='') as csvfile:
            reader = csv.reader(csvfile)
            rows = list(reader)
            if len(rows) < 1 or len(rows[-1]) < 3:
                raise ValueError("CSV file does not contain enough data to determine conversion factor.")
            print (f"Using conversion factor from CSV: {float(rows[-1][2])}")
            return float(rows[-1][2])  # Convert string to float
        

    def Flux_Calibrate(self):
        """
        Perform flux calibration on the data.
        This method should be implemented based on the specific requirements of the flux calibration process.
        """
        # Placeholder for flux calibration logic
        # This should include the actual implementation of flux calibration
        conversion_factor = self.get_conversion_factor()/10000000
        for i in range(len(self.data)):
            for j, row in enumerate(self.data[i]['DATA']):
                self.data[i]['DATA'][j] = row * conversion_factor  # Assign the result back


    def make_spec(self):
        channels = len(self.data)
        for ind, i in enumerate(self.data):
            # Only get actual data, not the calibration data
            data_start = self.data_start_ind // channels
            data_end = self.post_cal_start_ind // channels
            data = self.data[ind][data_start:data_end]

            # Get the summed intensities for each frequency
            # Sum across time (axis=0) and divide by number of time samples
            result = np.sum(data['DATA'], axis=0) / len(data)

            # Reverse the list because the data is in reverse order
            result = result[::-1]

            # Check if spectrum already has frequency data from user_cuts
            if ind < len(self.spectrum) and len(self.spectrum[ind]) > 0:
                # Use existing frequency array from user_cuts
                freq_array = self.spectrum[ind]
                # Ensure frequency array matches result length
                if len(freq_array) != len(result):
                    # Create new frequency array if lengths don't match
                    freq_array = np.linspace(self.freqs[ind][1], self.freqs[ind][0], len(result))
                else:
                    # Reverse to match spectrum order
                    freq_array = freq_array[::-1]
            else:
                # Create frequency array from start/stop frequencies
                # Reverse order (high to low) to match reversed spectrum
                freq_array = np.linspace(self.freqs[ind][1], self.freqs[ind][0], len(result))

            # Ensure spectrum list is large enough
            while len(self.spectrum) <= ind:
                self.spectrum.append([])

            # Store as [frequency_array, intensity_array]
            self.spectrum[ind] = [freq_array, result]

        return
    
    
    def make_continuum(self):

        for data in self.data:
            t0 = Time(self.header["DATE"], format="isot")
            times = Time(data["DATE-OBS"], format='isot')
            time_rel = (times - t0).sec

            # Get the summed intensities for each time
            result = np.sum(data['DATA'], axis=1) / len(data['DATA'][0])

            self.continuum.append([time_rel, result])

    def save(self, output_path=None):
        if output_path is None:
            base, ext = os.path.splitext(self.filepath)
            output_path = f"{base}_flux_calibrated{ext}"

        # Create Primary HDU with updated header
        primary_hdu = fits.PrimaryHDU(header=self.header)

        # Option 1: Save each feed/channel as separate HDUs
        hdu_list = [primary_hdu]
        for i, table in enumerate(self.data):
            # Clean string columns to avoid encoding issues
            cleaned_table = table.copy()
            for col_name in cleaned_table.colnames:
                if cleaned_table[col_name].dtype.kind in ['U', 'S']:  # Unicode or byte string
                    # Convert problematic characters to ASCII-safe versions
                    cleaned_data = []
                    for item in cleaned_table[col_name]:
                        if isinstance(item, bytes):
                            try:
                                cleaned_item = item.decode('utf-8', errors='replace')
                            except:
                                cleaned_item = str(item)
                        else:
                            # Replace non-ASCII characters
                            cleaned_item = str(item).encode('ascii', errors='replace').decode('ascii')
                        cleaned_data.append(cleaned_item)
                    cleaned_table[col_name] = cleaned_data
            
            table_hdu = fits.BinTableHDU(data=cleaned_table, name=f'FEED_{i}')
            hdu_list.append(table_hdu)
        
        # Combine HDUs
        hdulist = fits.HDUList(hdu_list)

        # Write to file
        try:
            hdulist.writeto(output_path, overwrite=True)
        except UnicodeEncodeError as e:
            print(f"Unicode encoding error: {e}")
            print("Attempting to save without problematic string columns...")
            
            # Fallback: save without string columns that cause problems
            hdu_list_backup = [primary_hdu]
            for i, table in enumerate(self.data):
                # Keep only numerical columns
                numeric_table = table.copy()
                cols_to_remove = []
                for col_name in numeric_table.colnames:
                    if numeric_table[col_name].dtype.kind in ['U', 'S']:
                        cols_to_remove.append(col_name)
                
                for col in cols_to_remove:
                    numeric_table.remove_column(col)
                
                table_hdu = fits.BinTableHDU(data=numeric_table, name=f'FEED_{i}')
                hdu_list_backup.append(table_hdu)
            
            hdulist_backup = fits.HDUList(hdu_list_backup)
            hdulist_backup.writeto(output_path, overwrite=True)
            print(f"Saved file without string columns: {output_path}")

    def conclude(self):
        self.user_cuts(self.params[0], self.params[1], self.params[2], self.params[3])
        self.Flux_Calibrate()
        self.make_spec()
        self.make_continuum()
        self.save()

if __name__ == "__main__":
    filepath = "C:/Users/starb/Downloads/0136872.fits"
    val = Val(filepath)
    val.validate()
    file = Gain_Calibrate(filepath.replace(".fits", "_validate.fits"))
    file.gain_calibrate()
    flux = Conclusion(filepath.replace(".fits", "_validate_gain_calibrated.fits"), [[12,40], [60, 70], [100, 110], [1392, 1400], [1420,1500]], "continuum", "cut", feeds=[0,1 ])
    flux.conclude()
    # print (flux.continuum[0][1])
    # print ("spectrum")
    # print (flux.spectrum[0][1])

    import matplotlib.pyplot as plt

    # Plot Spectrums and Continuums in one image with two subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

    # Spectrum plot
    axs[0, 0].plot(flux.spectrum[0][0], flux.spectrum[0][1], label='Spectrum 1')
    axs[0, 0].plot(flux.spectrum[1][0], flux.spectrum[1][1], label='Spectrum 2')
    axs[0, 0].set_xlabel('Frequency')
    axs[0, 0].set_ylabel('Intensity')
    axs[0, 0].set_title('Spectrums')
    axs[0, 0].legend()

    axs[1, 0].plot(flux.spectrum[2][0], flux.spectrum[2][1], label='Spectrum 3')
    axs[1, 0].plot(flux.spectrum[3][0], flux.spectrum[3][1], label='Spectrum 4')
    axs[1, 0].set_xlabel('Frequency')
    axs[1, 0].set_ylabel('Intensity')
    axs[1, 0].set_title('Spectrums')
    axs[1, 0].legend()

    # Continuum plot
    axs[0, 1].plot(flux.continuum[0][0], flux.continuum[0][1], label='Continuum 1')
    axs[0, 1].plot(flux.continuum[1][0], flux.continuum[1][1], label='Continuum 2')
    axs[0, 1].set_xlabel('Time (s)')
    axs[0, 1].set_ylabel('Intensity')
    axs[0, 1].set_title('Continuums')
    axs[0, 1].legend()

    axs[1, 1].plot(flux.continuum[2][0], flux.continuum[2][1], label='Continuum 3')
    axs[1, 1].plot(flux.continuum[3][0], flux.continuum[3][1], label='Continuum 4')
    axs[1, 1].set_xlabel('Time (s)')
    axs[1, 1].set_ylabel('Intensity')
    axs[1, 1].set_title('Continuums')
    axs[1, 1].legend()

    plt.savefig("spectrum_continuum_combined.png")
