
import numpy as np
import datetime as dt
def matlab2datetime(matlab_datenum):
        # equivalent of datestr when converting from serial number to date
       
        day = dt.datetime.fromordinal(int(matlab_datenum))
        dayfrac = dt.timedelta(days=float(matlab_datenum)%1) - dt.timedelta(days = 366)
        return day + dayfrac