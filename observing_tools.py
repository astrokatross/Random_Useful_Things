import os
import sys
import logging
 
from datetime import datetime, timedelta
 
import numpy as np
 
from astropy.coordinates import SkyCoord, Angle, AltAz, EarthLocation
import astropy.units as u
from astropy.table import Table
from astropy.time import Time, TimeDelta
 
mwa_location = EarthLocation(lat=-26.7031*u.degree, lon=116.671*u.degree, height=377.827*u.m) #mwa coordinates hard-coded here
 
atca_location = EarthLocation(lon = Angle(149.5501388, unit = u.degree),
                              lat = Angle('-30 18 46.385', unit = u.degree),
                              height = 236.87*u.m)
 
alma_location = EarthLocation(lat = Angle('-23 01 22.42', unit = u.degree), lon = Angle('-67 45 17.44', unit = u.degree), height = 5017*u.m)
 
gmrt_location = EarthLocation(lat = Angle('19 05 47.46', unit = u.degree), lon = Angle('74 02 59.07', unit = u.degree), height = 650.0*u.m)
 
parkes_location = EarthLocation(lon = Angle('148 15 54.636', unit = u.degree), lat = Angle('-32 59 48.636', unit = u.degree), height = 414.8*u.m)

siding_spring_location = EarthLocation.of_site('Siding Spring Observatory')

# vla_location = EarthLocation(lat = Angle('34 04 43.46', unit = u.degree), lon = Angle('74 02 59.07', unit = u.degree)()
 
vla_location = EarthLocation.of_site('vla')

class source:
    def __init__(self, ra, dec):
        if isinstance(ra, str):
            if ' ' in ra or 'h' in ra or ':' in ra:
                self.ra = Angle(ra, unit = u.h)
        elif '.' in str(ra):
            self.ra = Angle(ra, unit = u.deg)
        if isinstance(dec, str):
             
            if ' ' in dec or 'd' in dec or ':' in dec:
                self.dec = Angle(dec, unit = u.deg)
        elif '.' in str(dec):
            self.dec = Angle(dec, unit = u.deg)
 
        self.coords = SkyCoord(self.ra, self.dec, frame = 'icrs')
         
 
 
class observation_tools(source):#(mwa_voevent_parse, mwa_trigger_observation_parameters):
 
    def __init__(self, ra, dec, location, time = None):
 
        source.__init__(self, ra, dec)
         
        self.location = location #EarthLocation
         
        if time == None:
            self.time = Time(Time.now(), location = self.location)
        else:
            self.time = Time(time, format = 'iso', location = self.location)
             
                
        self.elevation = self.calculate_elevation()
        #self.is_up, self.min_time, self.max_time = self.determine_observing_window() #whether it is up within 12 hours of trigger event, the first time it is up after the trigger event, and the last time within 12 hours that it's up
 
 
    def calculate_elevation(self, time = None):
        """
        calculate the elevation of the target source
        """
 
        if time == None:
            time = self.time
             
        source_location = self.coords
 
        #calculate altitude and azimuth angles:
        source_Alt_Az = source_location.transform_to(AltAz(obstime = time, location = self.location))
         
        return source_Alt_Az.alt #astropy Latitude object
 
    def determine_observing_window(self, time = None, altitude_threshold = 45.0*u.degree):
 
        """
        determine the earliest and latest time that the trigger source location is above altitude_threshold (in degrees)
        within 12 hours of trigger time
        default alt threshold is 45 deg above horizon
        might want to consider case when source is only up for couple of minutes after trigger
        """
 
        if time is None:
            time = self.time #default value for time argument
 
        #utc_offset = 8.0*u.hour #AWST
        start_time = time.utc
        delta_time = np.linspace(0.0, 12.0, 200)*u.hour
        times = start_time + delta_time #an array of times spanning 12 hours after the trigger time
        source_altitudes = self.calculate_elevation(time = times) #will be a list of altitudes
        above_indexes = np.where(source_altitudes > altitude_threshold)[0] #where the altitudes are above altitude_threshold
        above_altitudes = source_altitudes[above_indexes]
        above_times = times[above_indexes]
 
        if len(above_times) == 0:
            return False, Time(-1, format = 'jd'), Time(-1, format = 'jd')
        else:
            #return True, when it is first up after the event, and the last time it is up 12 hours after event
            return True, above_times[0], above_times[-1]
 
 
