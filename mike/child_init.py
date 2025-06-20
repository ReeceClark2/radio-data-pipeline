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


class Radio_Child_File(Radio_File):
    def __init__(self, mike, file, index, axis, slice_type, feeds=[]):
        self.params = [index, axis, slice_type, feeds]

        try:
            with fits.open(file) as hdul:
                self.header = hdul[0].header
                self.data = Table(hdul[1].data)
        except Exception as e:
            raise MyException(f"Error reading FITS file: {e}")
        
        self.split_slp_feed()

        # self.data = copy.deepcopy(mike.data)
        # Copy other attributes as needed
        self.gain_start = copy.deepcopy(mike.gain_start)
        self.gain_end = copy.deepcopy(mike.gain_end)
        self.data_indicies = copy.deepcopy(mike.data_indicies)
        self.freqs = copy.deepcopy(mike.freqs)
        self.labels = copy.deepcopy(mike.labels)

        #TO be saved
        self.params = []
        self.continuum = copy.deepcopy(mike.continuum)
        self.spectrum = copy.deepcopy(mike.spectrum)
        # self.parent = file.file_path
        self.user_cuts(index, axis, slice_type, feeds)
        self.make_spec()
        print (self.params)



    def user_cuts(self, indicies, axis, slice_type, feeds=[]):
        """        Apply user-defined cuts to the data based on frequency or time intervals.
        indicies: list of tuples defining the intervals to keep or cut
        axis: "spectrum" for frequency cuts, "continuum" for time cuts
        type: "keep" to keep the specified intervals, "cut" to remove them
        feeds: list of feed numbers to apply the cuts to; if empty, applies to all feeds
        """        
        if feeds == []:
            feeds = np.arange(len(np.unique(d['IFNUM'] for d in self.data)))
        # Ensure indicies are in pairs


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
                for fmin, fmax in indicies:
                    if type == "keep":
                        freq_mask |= (freqs >= fmin) & (freqs <= fmax)
                    elif type == "cut":
                        freq_mask |= (freqs >= fmin) & (freqs <= fmax)
                if type == "cut":
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
                if type == "keep":
                    for tmin, tmax in indicies:
                        for j in range(len(time_rel)):
                            if time_rel[j] > tmin and time_rel[j] < tmax:
                                new_table.append(c[j])
                elif type == "cut":
                    # Start with all rows, then filter out those within any (tmin, tmax) interval
                    mask = np.ones(len(time_rel), dtype=bool)
                    for tmin, tmax in indicies:
                        for j in range(len(time_rel)):
                            if tmin < time_rel[j] < tmax:
                                mask[j] = False
                    new_table = [c[j] for j in range(len(time_rel)) if mask[j]]
                    # Convert new_table to an astropy Table only if it's not empty

                if new_table:
                    self.data[i] = Table(rows=new_table, names=c.colnames)
                else:
                    # If no rows matched, create an empty table with the same columns
                    self.data[i] = Table(names=c.colnames)


    def make_spec(self):
        for ind, i in enumerate(self.data):
            # Only get actual data, not the calibration data
            data = self.data[ind][self.data_indicies[ind][0]:self.data_indicies[ind][1]]

            # Get the summed intensities for each frequency
            result = average(data, 0)

            # Reverse the list becaise the data is in reverse order
            result = result[::-1]

            # Get the frequency range for the channel based on the feed
            # if the spectrum field has not been fully populated yet create the frequency range from the original start and stop frequencies
            frequencies = np.linspace(self.file.freqs[ind][0], self.file.freqs[ind][1], len(result))

            # Add the frequency and result to the file's spectrum field
            self.spectrum[ind] = ([np.array(frequencies), np.array(result)])
        
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
