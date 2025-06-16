import copy
from operator import index

from matplotlib.pyplot import axis
from file_init import Mike
import numpy as np
from astropy.table import Table
from file_exception import MyException
from astropy.time import Time

class Sully(Mike):
    def __init__(self, file, index, axis, type, feeds=[]):
        super().__init__(file.file_path)
        self.data = copy.deepcopy(file.data)
        # Copy other attributes as needed
        self.data_indicies = copy.deepcopy(file.data_indicies)
        self.freqs = copy.deepcopy(file.freqs)
        self.labels = copy.deepcopy(file.labels)

        #TO be saved
        self.params = []
        self.continuum = copy.deepcopy(file.continuum)
        self.spectrum = copy.deepcopy(file.spectrum)
        # self.parent = file.file_path
        self.user_cuts(index, axis, type, feeds)
        self.specMaker()



    def user_cuts(self, indices, axis, type, feeds=[]):
        """        Apply user-defined cuts to the data based on frequency or time intervals.
        indices: list of tuples defining the intervals to keep or cut
        axis: "spectrum" for frequency cuts, "continuum" for time cuts
        type: "keep" to keep the specified intervals, "cut" to remove them
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
                    for tmin, tmax in indices:
                        for j in range(len(time_rel)):
                            if time_rel[j] > tmin and time_rel[j] < tmax:
                                new_table.append(c[j])
                elif type == "cut":
                    # Start with all rows, then filter out those within any (tmin, tmax) interval
                    mask = np.ones(len(time_rel), dtype=bool)
                    for tmin, tmax in indices:
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

    def specMaker(self):
        for ind, i in enumerate(self.data):
            #get feed and polarization numbers
            feednum = np.unique(self.data[ind]["IFNUM"])[0]
            polnum = np.unique(self.data[ind]["PLNUM"])[0]
            # only get actual data, not the calibration data
            data = self.data[ind][self.data_indicies[ind][0]:self.data_indicies[ind][1]]

            # get the summed intensities for each frequency
            intensities = data['DATA']
            count = intensities.shape[0]
            result = np.sum(intensities, axis = 0) / count
            #reverse the list becaise the data is in reverse order
            result = result[::-1]

            #get the frequency range for the channel based on the feed
            # if the spectrum field has not been fully populated yet create the frequency range from the original start and stop frequencies
            frequencies = self.spectrum[ind][0]

            # Add the frequency and result to the file's spectrum field
            self.spectrum[ind] = ([np.array(frequencies), np.array(result)])
        return
