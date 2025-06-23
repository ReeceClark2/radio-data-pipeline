# Standard library
import copy
from operator import index

# Third-party libraries
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from matplotlib.pyplot import axis

# Local application imports
from file_exception import MyException
from file_init import Radio_File
from utils import average
from logger import Log_Collector


class Radio_Child_File(Radio_File):
    def __init__(self, radio_file, file, index, axis, slice_type, feeds=[]):
        self.params = [index, axis, slice_type, feeds]

        try:
            with fits.open(file) as hdul:
                self.header = hdul[0].header
                self.data = Table(hdul[1].data)
        except Exception as e:
            raise MyException(f"Error reading FITS file: {e}")
        
        self.split_slp_feed()

        # Copy other attributes as needed
        self.file_path = copy.deepcopy(radio_file.file_path)
        self.gain_start = copy.deepcopy(radio_file.gain_start)
        self.gain_end = copy.deepcopy(radio_file.gain_end)
        self.data_indices = copy.deepcopy(radio_file.data_indices)
        self.freqs = copy.deepcopy(radio_file.freqs)
        self.labels = copy.deepcopy(radio_file.labels)

        #TO be saved
        self.logger = Log_Collector(name=f"logger_{self.file_path}")
        self.params = [index, axis, slice_type, feeds]
        self.continuum = copy.deepcopy(radio_file.continuum)
        self.flux_calibrated = []
        self.spectrum = copy.deepcopy(radio_file.spectrum)
        # self.parent = file.file_path
        self.user_cuts(index, axis, slice_type, feeds)
        self.make_spec()


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
                    raise MyException("Data is not split by feed. Please run split_slp_feed() first.")
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
                self.spectrum[i][0] = selected_freqs


        if axis == "continuum":
            for i, c in enumerate(self.data):
                feednum = np.unique(c['IFNUM'])
                if feednum.size != 1:
                    raise MyException("Data is not split by feed. Please run split_slp_feed() first.")
                feednum = feednum[0]
                
                if feednum not in feeds:
                    continue

                times = Time(c["DATE-OBS"], format="isot")
                t0 = Time(self.header["DATE"], format="isot")
                time_rel = (times - t0).sec  # Time delta in seconds

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


    def make_spec(self):
        for ind, i in enumerate(self.data):
            # Only get actual data, not the calibration data
            data = self.data[ind][self.data_indices[ind][0]:self.data_indices[ind][1]]

            # Get the summed intensities for each frequency
            result = average(data, 0)

            # Reverse the list becaise the data is in reverse order
            result = result[::-1]

            # Get the frequency range for the channel based on the feed
            # if the spectrum field has not been fully populated yet create the frequency range from the original start and stop frequencies
            # Add the frequency and result to the file's spectrum field
            self.spectrum[ind] = ([np.array(self.spectrum[ind][0]), np.array(result)])
        
        return
    

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
