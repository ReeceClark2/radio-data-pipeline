# Third-party libraries
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import traceback

# Local application imports
from data_processing.file_init import Radio_File
from data_processing.flux_calibration import Flux_Cal
from data_processing.gain_calibration import Gain_Cal
from data_processing.sort import Sort
from data_processing.spectrum import Spectrum
from data_processing.val import Val
from data_processing.weather import Weather


if __name__ == "__main__":
    try:
        files = ["C:/Users/starb/Downloads/0136993.fits", "C:/Users/starb/Downloads/0137023.fits", "C:/Users/starb/Downloads/0137053.fits"]

        for i in files:
            print(i)
            file = Radio_File(i)

            v = Val(file)
            v.validate_primary_header()
            v.validate_data()

            s = Sort(file)
            s.sort()

            w = Weather(file)
            w.weather_correction()

            gc = Gain_Cal(file)
            gc.gain_cal()

            spec = Spectrum(file)
            spec.make_spec()

            for entry in file.logger.get_log_entries():
                print(f"{entry['level']}: {entry['message']}")

            file.logger.save(file.file_path)


        # fig, axs = plt.subplots(2, 2, figsize=(12, 8))

        # # Plot first two spectra (top-left)
        # for i in range(2):
        #     x, y = file.spectrum[i]
        #     axs[0, 0].plot(x, y, label=file.labels[i], linewidth=0.5)
        # axs[0, 0].set_xlabel('Frequency (MHz)')
        # axs[0, 0].set_ylabel('Intensity')
        # axs[0, 0].set_title('Spectra 1 & 2')
        # axs[0, 0].legend()

        # # Plot last two spectra (bottom-left)
        # for i in range(2, 4):
        #     x, y = file.spectrum[i]
        #     axs[1, 0].plot(x, y, label=file.labels[i], linewidth=0.5)
        # axs[1, 0].set_xlabel('Frequency (MHz)')
        # axs[1, 0].set_ylabel('Intensity')
        # axs[1, 0].set_title('Spectra 3 & 4')
        # axs[1, 0].legend()

        # # Plot first two continuum (top-right)
        # for i in range(2):
        #     x, y = file.continuum[i]
        #     axs[0, 1].plot(x, y, label=file.labels[i], linewidth=0.5)
        # axs[0, 1].set_xlabel('Time (Seconds)')
        # axs[0, 1].set_ylabel('Intensity (Jy)')
        # axs[0, 1].set_title('Continuum 1 & 2')
        # axs[0, 1].legend()

        # # Plot last two continuum (bottom-right)
        # for i in range(2, 4):
        #     x, y = file.continuum[i]
        #     axs[1, 1].plot(x, y, label=file.labels[i], linewidth=0.5)
        # axs[1, 1].set_xlabel('Time (Seconds)')
        # axs[1, 1].set_ylabel('Intensity (Jy)')
        # axs[1, 1].set_title('Continuum 3 & 4')
        # axs[1, 1].legend()

        # file = Radio_File("C:/Users/starb/Downloads/0135383.fits")
        # cal_file = Radio_File("C:/Users/starb/Downloads/0135383.fits")

        # v = Val(file)
        # v.validate_primary_header()
        # v.validate_data()

        # s = Sort(file)
        # s.sort()

        # w = Weather(file)
        # w.weather_correction()

        # gc = Gain_Cal(file)
        # gc.gain_cal()

        # v = Val(cal_file)
        # v.validate_primary_header()
        # v.validate_data()

        # s = Sort(cal_file)
        # s.sort()

        # w = Weather(cal_file)
        # w.weather_correction()

        # gc = Gain_Cal(cal_file)
        # gc.gain_cal()

        # # fc = Flux_Cal(file, cal_file)
        # # fc.flux_cal()

        # spec = Spectrum(file)
        # spec.make_spec()


        # # Plot first two spectra (top-left)
        # for i in range(2):
        #     x, y = file.spectrum[i]
        #     axs[0, 0].plot(x, y, label=file.labels[i], linewidth=0.5)
        # axs[0, 0].set_xlabel('Frequency (MHz)')
        # axs[0, 0].set_ylabel('Intensity')
        # axs[0, 0].set_title('Spectra 1 & 2')
        # axs[0, 0].legend()

        # # Plot last two spectra (bottom-left)
        # for i in range(2, 4):
        #     x, y = file.spectrum[i]
        #     axs[1, 0].plot(x, y, label=file.labels[i], linewidth=0.5)
        # axs[1, 0].set_xlabel('Frequency (MHz)')
        # axs[1, 0].set_ylabel('Intensity')
        # axs[1, 0].set_title('Spectra 3 & 4')
        # axs[1, 0].legend()

        # # Plot first two continuum (top-right)
        # for i in range(2):
        #     x, y = file.continuum[i]
        #     axs[0, 1].plot(x, y, label=file.labels[i], linewidth=0.5)
        # axs[0, 1].set_xlabel('Time (Seconds)')
        # axs[0, 1].set_ylabel('Intensity (Jy)')
        # axs[0, 1].set_title('Continuum 1 & 2')
        # axs[0, 1].legend()

        # # Plot last two continuum (bottom-right)
        # for i in range(2, 4):
        #     x, y = file.continuum[i]
        #     axs[1, 1].plot(x, y, label=file.labels[i], linewidth=0.5)
        # axs[1, 1].set_xlabel('Time (Seconds)')
        # axs[1, 1].set_ylabel('Intensity (Jy)')
        # axs[1, 1].set_title('Continuum 3 & 4')
        # axs[1, 1].legend()


        # plt.tight_layout()
        # plt.savefig('GB Ours', dpi=400)

    except Exception as e:
        if file is not None:
            file.logger.error("An exception occurred: " + str(e))
            file.logger.error(traceback.format_exc())
            file.logger.save(file.file_path)
        else:
            # Fallback: print and log to generic file if logger isn't available
            with open("crash_log.txt", "w") as f:
                f.write("Unhandled exception:\n")
                traceback.print_exc(file=f)
            print("An exception occurred before logger was initialized. See crash_log.txt for details.")


