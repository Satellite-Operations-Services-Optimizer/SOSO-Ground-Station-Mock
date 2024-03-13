## Imports
import json
import os
import numpy as np
import matplotlib.pyplot as plt
from math import atan, degrees, atan2
from math import radians, cos, sin, sqrt
from skyfield.api import load, EarthSatellite, Topos
from datetime import timedelta, datetime
from skyfield.sgp4lib import EarthSatellite
from math import radians, cos, sin, asin, sqrt
from skyfield import api
from skyfield.positionlib import ICRF, Geocentric
from skyfield.constants import (AU_M, ERAD, DEG2RAD,
                                IERS_2010_INVERSE_EARTH_FLATTENING, tau)
from skyfield.units import Angle

from numpy import einsum, sqrt, arctan2, pi, cos, sin
# from mock_SAT_ver2 import load_satellites, ephemeris_power_schedule_start_end, location_satellite, eclipse_sunlight, image_validation

## Step 1: Load the TLE files
ts = load.timescale() # Create timescale object for TLE computation
satellites = [] # Empty list of satellite objects (SOSO-1, SOSO-2, etc)

# Create empty lists for the latitudes and longitudes of each satellite
latitudes = [[] for _ in range(5)]
longitudes = [[] for _ in range(5)]

for i in range(1, 6): # Iterate over numbers 1 to 5
    try: # Used to catch and handle exceptions in the code
        with open(f'SOSO-{i}_TLE.json') as f: # For-loop and f-string used to open the TLE files for SOSO-1, SOSO-2, etc.
            data = json.load(f) # Load the JSON data from the file
            name = data['name']
            line1 = data['line1']
            line2 = data['line2']
        satellite = EarthSatellite(line1, line2, name, ts) # Create new satellite object where line 1 = tle[1], line 2 = tle[2], title = tle[0], and ts for timescale
        satellites.append(satellite) # Add each satellite to the empty list
    except IndexError: # Handles TLE files without title line or missing lines
        print(f"Error: TLE file for satellite {i} is not formatted correctly.") # Output of error

## Step 2: (Maintenance) Is the satellite in eclipse or in sunlight?
eph = load('de421.bsp')  # Load the JPL ephemeris DE421
sun = eph['sun']  # Get the 'sun' object from the ephemeris
earth = eph['earth']  # Get the 'earth' object from the ephemeris

## Step 3: Select time interval for satellite and ground station accesses.

# Ask the user for the start and end times
start_time_str = input("Enter the start time (YYYY-MM-DD HH:MM:SS): ")
end_time_str = input("Enter the end time (YYYY-MM-DD HH:MM:SS): ")

# Convert the input strings to datetime objects
start_time = datetime.strptime(start_time_str, "%Y-%m-%d %H:%M:%S")
end_time = datetime.strptime(end_time_str, "%Y-%m-%d %H:%M:%S")

# Convert the datetime objects to skyfield Time objects
start_time_skyfield = ts.utc(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute, start_time.second)
end_time_skyfield = ts.utc(end_time.year, end_time.month, end_time.day, end_time.hour, end_time.minute, end_time.second)

# Convert the datetime objects to skyfield Time objects
start_time_skyfield = ts.utc(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute, start_time.second)
end_time_skyfield = ts.utc(end_time.year, end_time.month, end_time.day, end_time.hour, end_time.minute, end_time.second)

## Reference: https://stackoverflow.com/questions/54969534/skyfield-visible-area-underneath-earthsatellite

def reverse_terra(xyz_au, gast, iterations=3):
    """Convert a geocentric (x,y,z) at time `t` to latitude and longitude.
    Returns a tuple of latitude, longitude, and elevation whose units
    are radians and meters.  Based on Dr. T.S. Kelso's quite helpful
    article "Orbital Coordinate Systems, Part III":
    https://www.celestrak.com/columns/v02n03/
    """
    x, y, z = xyz_au
    R = sqrt(x*x + y*y)

    lon = (arctan2(y, x) - 15 * DEG2RAD * gast - pi) % tau - pi
    lat = arctan2(z, R)

    a = ERAD / AU_M
    f = 1.0 / IERS_2010_INVERSE_EARTH_FLATTENING
    e2 = 2.0*f - f*f
    i = 0
    C = 1.0
    while i < iterations:
        i += 1
        C = 1.0 / sqrt(1.0 - e2 * (sin(lat) ** 2.0))
        lat = arctan2(z + a * C * e2 * sin(lat), R)
    elevation_m = ((R / cos(lat)) - a * C) * AU_M
    earth_R = (a*C)*AU_M
    return lat, lon, elevation_m, earth_R

def subpoint(self, iterations):
    """Return the latitude an longitude directly beneath this position.

    Returns a :class:`~skyfield.toposlib.Topos` whose ``longitude``
    and ``latitude`` are those of the point on the Earth's surface
    directly beneath this position (according to the center of the
    earth), and whose ``elevation`` is the height of this position
    above the Earth's center.
    """
    if self.center != 399:  # TODO: should an __init__() check this?
        raise ValueError("you can only ask for the geographic subpoint"
                            " of a position measured from Earth's center")
    t = self.t
    xyz_au = einsum('ij...,j...->i...', t.M, self.position.au)
    lat, lon, elevation_m, self.earth_R = reverse_terra(xyz_au, t.gast, iterations)

    from skyfield.toposlib import Topos
    return Topos(latitude=Angle(radians=lat),
                    longitude=Angle(radians=lon),
                    elevation_m=elevation_m)

def earth_radius(self):
    return self.earth_R

def satellite_visible_area(earth_radius, satellite_elevation):
    """Returns the visible area from a satellite in square meters.

    Formula is in the form is 2piR^2h/R+h where:
        R = earth radius
        h = satellite elevation from center of earth
    """
    return ((2 * pi * ( earth_radius ** 2 ) * 
            ( earth_radius + satellite_elevation)) /
            (earth_radius + earth_radius + satellite_elevation))
    
def satellite_visible_square_area(earth_radius, satellite_elevation):
    """Returns the visible square area from a satellite in square kilometers.

    Formula is in the form is s^2 where:
        s = 2 * sqrt(2Rh + h^2)
        R = earth radius
        h = satellite elevation from center of earth
    """
    s = sqrt(2 * earth_radius * satellite_elevation + satellite_elevation ** 2)
    return s ** 2

# Filter satellites based on your criteria (e.g., name)
selected_satellites = [satellite for satellite in satellites if satellite.name.startswith('SOSO')]

# Iterate over selected satellites
for satellite in selected_satellites:
    try:
        # Get position of the satellite at the start time
        geocentric = satellite.at(start_time_skyfield)
        geocentric.subpoint = subpoint.__get__(geocentric, Geocentric)
        geocentric.earth_radius = earth_radius.__get__(geocentric, Geocentric)
        
        # Get subpoint details
        geodetic_sub = geocentric.subpoint(3)
        print(f'Satellite: {satellite.name}')
        print('Geodetic latitude:', geodetic_sub.latitude)
        print('Geodetic longitude:', geodetic_sub.longitude)
        print('Geodetic elevation (km)', int(geodetic_sub.elevation.km))
        print('Geodetic earth radius (km)', int(geocentric.earth_radius() / 1000))

        geocentric_sub = geocentric.subpoint(0)
        
        # Calculate the ground distance (in meters) from the satellite
        ground_distance = sqrt(geocentric_sub.elevation.m ** 2 + geocentric_sub.elevation.m ** 2)

        # Calculate the angle of elevation (in radians)
        angle_of_elevation_rad = atan(geocentric_sub.elevation.m / ground_distance)

        # Convert the angle to degrees
        angle_of_elevation_deg = degrees(angle_of_elevation_rad)
        
        print('Geocentric latitude:', geocentric_sub.latitude)
        print('Geocentric longitude:', geocentric_sub.longitude)
        print('Geocentric elevation / Altitude (km)', int(geocentric_sub.elevation.km))
        print('Geocentric earth radius (km)', int(geocentric.earth_radius() / 1000))
        print('Angle of Elevation (degrees)', angle_of_elevation_deg)
        print('Visible area (km^2) SPHERICAL CAP', satellite_visible_area((geocentric.earth_radius() / 1000), geocentric_sub.elevation.km))
        print('Visible area (km^2) SQUARE AREA', satellite_visible_square_area((geocentric.earth_radius() / 1000), geocentric_sub.elevation.km))
        
        print('---')
    except Exception as e:
        print(f"Error processing satellite {satellite.name}: {e}")

# Haversine formula to calculate the distance between two points on the Earth's surface
def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0  # Radius of the Earth in kilometers
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * atan2(sqrt(a), sqrt(1-a))
    distance = R * c
    return distance

# Function to calculate the distance to the horizon from the satellite's elevation
def distance_to_horizon(elevation):
    R = 6371.0  # Radius of the Earth in kilometers
    # Distance to the horizon in kilometers
    horizon_distance = sqrt((R + elevation)**2 - R**2)
    return horizon_distance

# Function to calculate the approximate field of view area
def field_of_view_area(elevation):
    horizon_distance = distance_to_horizon(elevation)
    # Approximate the field of view as a circular segment with the horizon distance as the radius
    area = pi * (horizon_distance ** 2)
    return area

# Iterate over selected satellites
for satellite in selected_satellites:
    try:
        # Get position of the satellite at the start time
        geocentric = satellite.at(start_time_skyfield)
        geocentric.subpoint = subpoint.__get__(geocentric, Geocentric)
        geocentric.earth_radius = earth_radius.__get__(geocentric, Geocentric)
        
        # Get subpoint details
        geodetic_sub = geocentric.subpoint(3)
        
        # Calculate the field of view area
        fov_area = field_of_view_area(geodetic_sub.elevation.km)
        
        # Output the calculated field of view area
        print(f'Satellite: {satellite.name}')
        print(f'Field of view area (km^2): {fov_area}')
        
        print('---')
    except Exception as e:
        print(f"Error processing satellite {satellite.name}: {e}")
