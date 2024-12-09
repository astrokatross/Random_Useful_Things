import astropy
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.time import Time

from observing_tools import *

from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle
#from astral import Astral


def deg2HMS(ra='', dec='', round=False):
  RA, DEC, rs, ds = '', '', '', ''
  if dec:
    if str(dec)[0] == '-':
      ds, dec = '-', abs(dec)
    deg = int(dec)
    decM = abs(int((dec-deg)*60))
    if round:
      decS = int((abs((dec-deg)*60)-decM)*60)
    else:
      decS = (abs((dec-deg)*60)-decM)*60
    DEC = '{0}{1} {2} {3}'.format(ds, deg, decM, decS)
  
  if ra:
    if str(ra)[0] == '-':
      rs, ra = '-', abs(ra)
    raH = int(ra/15)
    raM = int(((ra/15)-raH)*60)
    if round:
      raS = int(((((ra/15)-raH)*60)-raM)*60)
    else:
      raS = ((((ra/15)-raH)*60)-raM)*60
    RA = '{0}{1} {2} {3}'.format(rs, raH, raM, raS)
  
  if ra and dec:
    return (RA, DEC)
  else:
    return RA or DEC

# print(deg2HMS(ra=48.53268756992491, dec=-34.441721004670704))

flora_src = observation_tools("17h 28m 12.13s", "-46d 08m 01s", parkes_location)
# source1 = observation_tools('17h 02m 38.976s', '21d 34m 35.91s', siding_spring_location)
# source2 = observation_tools('3h 14m 7.845016781979055s', '-34d 26m 30.195616814534674s', siding_spring_location)
# source3 = observation_tools('02h 48m 38s', '-32d 13m 36s', siding_spring_location)
# source4 = observation_tools('03h 30m 23s', '-07d 40m 52s', siding_spring_location)

# calibrator = observation_tools('19h 39m 25s', '-63d 42m 45s', siding_spring_location)
# calibrator = observation_tools('17h 02m 38.98s', '+21d 34m 35.91s', vla_location)

d = Time('2024-03-23 00:00:00', format = 'iso') #2019-04-30 12:00:00 AEST
times = d + np.arange(0.0, 36.0, 0.01)*u.hour
source1_alts = np.asarray(flora_src.calculate_elevation(time = times))
# source2_alts = np.asarray(source2.calculate_elevation(time = times))
# source3_alts = np.asarray(source3.calculate_elevation(time = times))
# source4_alts = np.asarray(source4.calculate_elevation(time = times))

# anu_alts = np.asarray(calibrator.calculate_elevation(time = times))

# t = Table.read('/data/random_useful_code/useful-code/target_list.csv', format = 'ascii.csv')


# standards = [observation_tools(ra, dec, atca_location) for ra, dec in t['ra', 'dec']]
# standard_names = t['Name']
# standard_alts = [np.asarray(standard.calculate_elevation(time = times)) for standard in standards]


plt.plot((times - d).to(u.hour), source1_alts, label = 'Src2', c = 'C6', linewidth =2.5)
# plt.plot((times - d).to(u.hour), source2_alts, label = 'Src1', c = 'C4', linewidth =2.5)
# # plt.plot((times - d).to(u.hour), source3_alts, label = '024838', c = 'C9', linewidth =2.5)
# plt.plot((times - d).to(u.hour), source4_alts, label = '033023', c = 'orange', linewidth =2.5)

# plt.plot((times - d).to(u.hour), anu_alts, label = '1934', c = 'r', linewidth =2.5)


# i = 0
# for standard_alt,standard_name in zip(standard_alts, standard_names):
#     plt.plot((times - d).to(u.hour), standard_alt, label = standard_name, c = 'k', linewidth =2.5)
#     # if i < 10:
    #     if t[i]['P'] == 'P':
    #         plt.plot((times -d).to(u.hour), standard_alt, label = standard_name, linewidth = 2.5)
    #     else:
    #         plt.plot((times -d).to(u.hour), standard_alt, label = standard_name)
    # else:
    #     if t[i]['P'] == 'P':
    #         plt.plot((times -d).to(u.hour), standard_alt, '--', label = standard_name, linewidth = 2.5)
    #     else:
    #         plt.plot((times -d).to(u.hour), standard_alt, '--', label = standard_name)

    # i+=1
# plt.axvline(14, c = 'k')
# plt.axvline(0, c = 'k')
        
plt.legend()
plt.xlim(0,24)
plt.ylim(12, 90)
# plt.show()wifes_srcs_223933
# plt.legend()
# plt.axvline((mwa_sunset - d).to(u.hour).value, c= 'C1')
# plt.axvline((mwa_sunrise - d).to(u.hour).value, c = 'C1')
# plt.axvline((anu_sunset - d).to(u.hour).value, c = 'C0')
# plt.axvline((anu_sunrise - d).to(u.hour).value, c = 'C0')
# plt.axvline((mwa_uptimes[0] - d).to(u.hour).value, c = 'C0')
# plt.axvline((mwa_uptimes[-1] - d).to(u.hour).value, c = 'C0')

# plt.axhline(45.0)
plt.axhline(30.0)
plt.axvline(1.)
plt.xlabel('Time from {} UT (hour)'.format((d).iso))
plt.ylabel('Source elevation')
# plt.xlim(2, 60)
plt.savefig("/data/random_useful_code/useful-code/flora_src.png")
# plt.show()

#print((0.03*u.hour).to(u.minute))

# J032237 	03h 22m 37s	-48d 20m 16s
