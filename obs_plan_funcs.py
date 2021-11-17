#importing useful packages
import numpy as np
import matplotlib.pyplot as plt
#---
from astropy import units as u
from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import get_sun, get_moon, SkyCoord, EarthLocation, AltAz
#---
from astroplan import plots, observability_table
from astroplan import Observer, FixedTarget
from astroplan import AltitudeConstraint, AirmassConstraint, AtNightConstraint
from astroplan import is_observable, is_always_observable, months_observable
#---
from astroquery.jplhorizons import Horizons
#---
import os
import pytz


def obs_times(obs_start_date,obs_end_date,night_start_time,night_end_time,utc_offset):
    obs_start_time = str(obs_start_date+' '+night_start_time)
    obs_start_time = Time(obs_start_time) - utc_offset 

    obs_end_time = str(obs_end_date+' '+night_end_time) 
    obs_end_time = Time(obs_end_time) - utc_offset  

    midnight = str(obs_end_date+' '+'00:00:00')   
    midnight = Time(midnight) - utc_offset     

    return obs_start_time, obs_end_time, midnight

def obs_stuff(longitude, latitude, elevation, obs_start_time, obs_end_time):
    # Specify the location of your observatory:
    obs_location = EarthLocation.from_geodetic(longitude, latitude, elevation)

    # Make an observer at your observatory:
    observatory = Observer(location=obs_location)

    # making time range of osbervation period
    obs_time_range = Time([obs_start_time, obs_end_time])

    return obs_location, observatory, obs_time_range


def file_reader(targets_file):
    target_names = []
    with open(targets_file, "r") as file:
        for line in file: 
            line = line.strip("\n") 
            target_names.append(line) 
    file.close()
    return target_names


def jpl_querier(target_names,obs_location_code,epoch_dates):
    target_info = []
    for target_name in target_names:
        obj = Horizons(id=target_name,location=obs_location_code,epochs=epoch_dates)
        obj_data = obj.ephemerides(skip_daylight=True)
        JPLtableData = obj_data['targetname','datetime_jd','datetime_str','RA','DEC','airmass','AZ','EL']
        target_info.append(JPLtableData)
    return target_info


def coord_obs(target_info):
    fixed_targets_lst = []
    skycoord_targets_lst = []

    for target in target_info:
        target_name = target['targetname']
        my_ra = target['RA']
        my_dec = target['DEC']
        jd_times = target['datetime_jd']
        obs_time_jd = Time(jd_times,format='jd',scale='utc')

        # making SkyCoord object
        c = SkyCoord(ra=my_ra,dec=my_dec,unit=(u.hourangle,u.deg),obstime=obs_time_jd)
        skycoord_targets_lst.append(c)
        
        # making FixedObject object
        fixed_target = FixedTarget(coord=c,name=target_name)
        fixed_targets_lst.append(fixed_target)
    
    return skycoord_targets_lst, fixed_targets_lst


def altaz_info(skycoord_targets_lst,midnight,obs_location):
    # defining midnight for calculation purposes
    delta_midnight = np.linspace(-12, 12, len(skycoord_targets_lst[0])) * u.hour 
    times = midnight + delta_midnight

    # creating an Alt/Az frame for calculation purposes
    altazframe = AltAz(obstime=times, location=obs_location)

    # finding the Alt/Az information for the Sun and Moon
    sunaltazs = get_sun(times).transform_to(altazframe)
    moonaltazs = get_moon(times).transform_to(altazframe)

    # finding the Alt/Az information for your targets
    target_altazs = []
    for target in skycoord_targets_lst:
        target_altazs.append(target.transform_to(altazframe))

    return delta_midnight, sunaltazs, moonaltazs, target_altazs


def zerolines(lim_alt):
    # making a zero-line for the limiting alitude of your telescope
    zeroLine_x1=np.arange(-12,90,10)
    zeroLine_y1=np.full(len(zeroLine_x1),lim_alt)

    # making a warning zero-line for the limiting alitude of your telescope
    zeroLine_y2=np.full(len(zeroLine_x1),lim_alt+10)

    return zeroLine_x1, zeroLine_y1, zeroLine_y2


def plot_alt(lim_alt,delta_midnight,sunaltazs,moonaltazs,target_altazs,target_names,plots_save_path,obs_dates_display,obs_dates_save,j):
    zeroLine_x1, zeroLine_y1, zeroLine_y2 = zerolines(lim_alt)
    plt.figure(figsize=(20,10))
    plt.plot(zeroLine_x1,zeroLine_y1,"-",linewidth=7,color="deeppink",label="limiting altitude")
    plt.plot(zeroLine_x1,zeroLine_y2,"-",linewidth=2,color="deeppink",label="warning altitude")
    plt.plot(delta_midnight, moonaltazs.alt, color=[0.75]*3, ls='--', label='Moon')

    # plotting the alitudes of the objects
    for i in range(len(target_altazs)):  
    #     plt.scatter(delta_midnight,target_altazs[i].alt,label='{}'.format(target_names[i]), lw=0, s=8) 
        plt.plot(delta_midnight,target_altazs[i].alt,label='{}'.format(target_names[i])) 

    # defining twilight/nighttimes
    # plt.fill_between(delta_midnight, 0, 90, sunaltazs.alt < -0*u.deg, color='0.5', zorder=0)  
    # plt.fill_between(delta_midnight, 0, 90, sunaltazs.alt < -18*u.deg, color='k', zorder=0)   
    plt.fill_between(delta_midnight, 0, 90, sunaltazs.alt < -0*u.deg, color='sandybrown', zorder=0, alpha = 0.5)  
    plt.fill_between(delta_midnight, 0, 90, sunaltazs.alt < -6*u.deg, color='cornflowerblue', zorder=0, alpha = 0.5) 
    plt.fill_between(delta_midnight, 0, 90, sunaltazs.alt < -12*u.deg, color='darkblue', zorder=0, alpha = 0.5)  
    plt.fill_between(delta_midnight, 0, 90, sunaltazs.alt < -18*u.deg, color='k', zorder=0, alpha = 0.9)  

    plt.legend()
    plt.xlim(-3, 6)          #can change this as appropriate, it is the number of hours before and after midnight (ie on your x-axis)
    plt.ylim(lim_alt-5, 90)  #could range from 0-90 degrees
    plt.xlabel('Hours from Midnight')  
    plt.ylabel('Altitude [deg]')  
    plt.title("Object Visibility on {}".format(obs_dates_display))
    plt.grid()
    plt.savefig(plots_save_path/'TargetVisibility_{}-{}.jpeg'.format(obs_dates_save,j),dpi=900)
    plt.show()

def plot_target_position(fixed_targets_lst,observatory, target_position_plot_times, plots_save_path, obs_dates_display, obs_dates_save,j):
    plt.figure(figsize=(10,10))
    for target in fixed_targets_lst:
        plots.plot_sky(target, observatory, target_position_plot_times)
    plt.title("Target Positions on {}".format(obs_dates_display))
    # plt.legend()
    plt.savefig(plots_save_path/'TargetPositions_{}-{}.jpeg'.format(obs_dates_save,j),dpi=900)
    plt.show()