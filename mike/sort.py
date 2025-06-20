# Third-party libraries
import numpy as np
from astropy.table import Table

# Local application imports
from file_init import Radio_File


class Sort:
    def __init__(self, file):
        '''
        Initialize file and confirm header and data have been validated.
        '''

        self.file = file

        # if self.file.validated_header == False:
        #     raise MyException('FITS header has not been validated!')
        # if self.file.validated_data == False:
        #     raise MyException('FITS data has not been validated!')


    def split_slp_feed(self):
        '''
        Split the data of the provided file by channel and feed.
        '''

        ifnums = np.unique(self.file.data["IFNUM"])
        plnums = np.unique(self.file.data["PLNUM"])

        self.file.continuum.extend([None for i in range(len(ifnums) * len(plnums))])
        self.file.spectrum.extend([None for i in range(len(ifnums) * len(plnums))])

        data = []
        labels = []
        for i in ifnums:
            for j in plnums:
                subset_mask = (self.file.data["IFNUM"] == i) & (self.file.data["PLNUM"] == j)
                subset_data = self.file.data[subset_mask]
                data.append(subset_data)
                labels.append(f'Feed{i + 1},Channel{j + 1}')
        self.file.data = data
        self.file.labels = labels

        return 


    def divide_section(self, section):
        '''
        Locate indicies within some channel of the data.

        Params:
        section: astropy FITS table

        Returns:
        array: 1st index of data, 1st index of post calibration
        '''

        # Initialize data start index and post calibration start index
        data_start_ind = None
        post_cal_start_ind = None

        # Create counter for limiting 'false data'
        counter = 0

        # Initialize calibration not started and pre calibration not completed
        cal_started = False
        pre_cal_complete = False

        for ind, point in enumerate(section):
            if point["CALSTATE"] == 1:
                cal_started = True

            if cal_started and point["CALSTATE"] == 0 and point["SWPVALID"] == 1 and not pre_cal_complete:
                data_start_ind = ind
                pre_cal_complete = True

            if pre_cal_complete and point["SWPVALID"] == 0 and section[ind-1]["SWPVALID"] == 0:
                if post_cal_start_ind is None:
                    post_cal_start_ind = ind
            else:
                post_cal_start_ind = None

            if pre_cal_complete and point["CALSTATE"] == 0 and point["SWPVALID"] == 1:
                counter += 1

            if counter <= 3 and point["CALSTATE"] == 1:
                pre_cal_complete = False

            if pre_cal_complete and point["SWPVALID"] == 0 and point["CALSTATE"] == 1:
                break


        if data_start_ind is None:
            for ind, point in enumerate(section):
                if point["SWPVALID"] == 1 and data_start_ind is None:
                    data_start_ind = ind
                    continue  

                if data_start_ind is not None and point["SWPVALID"] == 0:
                    post_cal_start_ind = ind
                    break


        if self.file.header["OBSMODE"] == "onoff":
            for ind, point in enumerate(section):
                target = 'onoff:off'
                if target in point["OBSMODE"]: 
                    offstart = ind
                    indicies = np.array([data_start_ind, offstart-1, offstart, post_cal_start_ind])   
                    break
                            
        else: 
             indicies = np.array([data_start_ind, post_cal_start_ind])

        return indicies


    def sort_data(self):
        '''
        Iterate through each channel of data and append to data_indicies.
        '''

        for i in self.file.data:
            result = self.divide_section(i)
            self.file.data_indicies.append(result)

        return
    

    def get_startend_freqs(self):
        '''
        Get the start and stop frequencies for each channel in the data.
        param file: Radio_File class file
        returns: populates the file's freqs field with start and stop frequencies
        '''
        freqs = []
        
        # Scour the header for the bandwidth and center frequencies
        for key, value in self.file.header.items():
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

        channels = len(np.unique([d['PLNUM'] for d in self.file.data]))
        for c in center:
            start_f = c - (band / 2)
            stop_f = c + (band / 2)
            
            for _ in range(channels):
                freqs.append(np.array([start_f, stop_f]))
        self.file.freqs = freqs

    
    def get_startstop_channels(self):
        '''
        Get the start and stop channels for each channel in the data.
        param file: Radio_File class file
        returns: cuts out the channels in the data that are not in the start and stop channels
        '''

        for i in range(len(self.file.data)):
            # Search through the header 
            for key, value in self.file.header.items():
                if key == ("HISTORY"):
                    if value.startswith("START,STOP"):
                        # Extract all integers from the string
                        value = value.replace(",", " ").strip()
                        value = value.split(" ")
                        # Split the value into individual words and numbers

                        channels = []
                        for k in value:
                            k = str(k).strip()
                            try:
                                # If it's an integer add it to the channels list
                                k = int(k)
                                channels.append(k)
                            except ValueError:
                                # If it can't become a integer, skip it
                                continue
                        
                        # Remove all string type characters from the channels list
                        channel1 = int(channels[0])
                        channel2 = int(channels[1])
            
            # Cut the DATA column to only include the channels in the start and stop channels
            t = Table(self.file.data[i])
            t.replace_column('DATA', np.array([row[channel1:channel2] for row in t['DATA']]))
            self.file.data[i] = t

        return
    

    def sort(self):
        self.split_slp_feed()
        self.sort_data()
        self.get_startend_freqs()
        self.get_startstop_channels()

        return


    def section_debug(self):
        '''
        Iterate through data and check the pre calibration, data, and post calibration of each.
        '''

        for ind, i in enumerate(self.file.data):
            current_indicies = self.file.data_indicies[ind]
            data_start_ind = current_indicies[0]
            post_cal_start_ind = current_indicies[1]

            pre_cal = i[:data_start_ind]
            data = i[data_start_ind:post_cal_start_ind]
            post_cal = i[post_cal_start_ind:]

            print(data_start_ind, post_cal_start_ind)
            print(len(pre_cal),len(data),len(post_cal))
            print (data_start_ind, post_cal_start_ind)


if __name__ == "__main__":
    file = Radio_File("C:/Users/starb/Downloads/0136870.fits")

    np.set_printoptions(threshold=100000)
    keep_indicies = [[1394,1400], [1401, 1402]]  # Example indicies to keep

    s = Sort(file)
    s.sort()
    s.section_debug()

