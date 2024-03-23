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
        
#############################################################################################
########################### Plotting the Field of View Area ################################
#############################################################################################

# Assuming the viewing angle is 30 degrees (half of the full FoV of 60 degrees)
viewing_angle_deg = 30
viewing_angle_rad = radians(viewing_angle_deg)

# Earth's radius in kilometers
earth_radius_km = 6378.137

for i, satellite in enumerate(satellites):  # Loop over the satellites in the list.
    current_time_skyfield = start_time_skyfield
    while current_time_skyfield.tt < end_time_skyfield.tt:  # Compare Julian dates
        position_current = satellite.at(current_time_skyfield).position.km  # Plain (x, y, z) coordinates at the current time (Center of Earth)
        subpoint_current = satellite.at(current_time_skyfield).subpoint()
        latitude_current = subpoint_current.latitude.degrees  # Latitude at the current time
        longitude_current = subpoint_current.longitude.degrees  # Longitude at the current time
        
        # Get the positions of the Earth, Sun, and satellite
        earth_pos = earth.at(current_time_skyfield).position.km
        sun_pos = sun.at(current_time_skyfield).position.km

        # Append the latitude and longitude to their respective lists
        latitudes[i].append(latitude_current)
        longitudes[i].append(longitude_current)

        # Calculate altitude from position data
        semi_major_axis_km = satellite.model.a * 6378.137  # Get the semi-major axis in kilometers
        altitude_current = semi_major_axis_km - 6378.137  # The altitude is the semi-major axis minus the Earth's radius

        # Calculate FOV
        # fov_current = degrees(2 * atan(12742 / (2 * (altitude_current + 6371))))
        # No need to calculate FOV since it is given to us.
        
        # Calculate the distance from the satellite to the horizon
        # This calculation assumes a simple spherical Earth model
        distance_to_horizon = sqrt((altitude_current)**2)

        # Calculate the FoV width and height on the ground using the viewing angle
        # The FoV dimensions are calculated as twice the distance from the satellite to the edge of the FoV
        # This uses the tangent of the viewing angle and the distance to the horizon
        fov_width = 2 * distance_to_horizon * np.tan(viewing_angle_rad)
        fov_height = fov_width  # Assuming a square FoV for simplicity

        print(f"SOSO-{i + 1} at Current Time:")
        print(f"  Position: {position_current} km")
        print(f"  Latitude: {latitude_current} degrees")
        print(f"  Longitude: {longitude_current} degrees")
        print(f"  Altitude: {altitude_current} km")
        print(f"  Satellite {satellite.name}: FoV width: {fov_width:.2f} km, FoV height: {fov_height:.2f} km")
        # print(f"  Field of View: {fov_current} degrees")
        print(f"  Time: {current_time_skyfield.utc_strftime('%Y %b %d %H:%M:%S')}")

        current_time_skyfield = ts.utc(current_time_skyfield.utc_datetime() + timedelta(minutes=1)) # Print all variables every minute from start and end times.
        
# Define image types and their properties
image_types = {
    "Spotlight": {"height": 10, "width": 10, "transfer_time": 120},
    "Medium": {"height": 40, "width": 20, "transfer_time": 45},
    "Low": {"height": 40, "width": 20, "transfer_time": 20}
}

# Define the viewing angle
viewing_angle = 30  # degrees

# Iterate over the image orders
for i in range(1, 50):
    with open(f'SampleOrders/Order_{i}.json') as f:
        order = json.load(f)

    # Get the image type
    image_type = image_types[order["ImageType"]]

    # Calculate the FOV square area
    for satellite in satellites:
        fov_current = fov_width # km

        # Check if the image order falls within the FOV
        if image_type["height"] <= fov_current and image_type["width"] <= fov_current:
            # Calculate the start and end times
            start_time = datetime.strptime(order["ImageStartTime"], "%Y-%m-%dT%H:%M:%S")
            end_time = datetime.strptime(order["ImageEndTime"], "%Y-%m-%dT%H:%M:%S")
            
            # Convert the datetime objects to skyfield Time objects
            start_time_skyfield = ts.utc(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute, start_time.second)
            end_time_skyfield = ts.utc(end_time.year, end_time.month, end_time.day, end_time.hour, end_time.minute, end_time.second)

            # Get the satellite's position at the start and end times
            start_position = satellite.at(start_time_skyfield).position.km
            end_position = satellite.at(end_time_skyfield).position.km

            # Check if the satellite stays within the FOV during the transfer time (Capture Time)
            transfer_end_time = start_time + timedelta(seconds=image_type["transfer_time"])
            if start_time <= transfer_end_time <= end_time:
                if start_position[0] <= end_position[0] and start_position[1] <= end_position[1]: # x-axis (start_position[0]) and y-axis (start_position[1]) of satellite
                    print(f"Order {i} is acceptable for {satellite.name}.")
                    print(f"Transfer completed at time {transfer_end_time}, lat {order['Latitude']}, lon {order['Longitude']}")
                    
                    # Plot the satellite's position
                    plt.figure(figsize=(10, 10))
                    plt.plot([order['Longitude']], [order['Latitude']], 'ro')
                    plt.xlim(order['Longitude'] - 1, order['Longitude'] + 1)
                    plt.ylim(order['Latitude'] - 1, order['Latitude'] + 1)
                    plt.grid(True)
                    plt.title(f"Satellite {satellite.name} Position for Order {i}")
                    plt.xlabel("Longitude")
                    plt.ylabel("Latitude")
                    plt.show()
                else:
                    print(f"Order {i} is unacceptable for {satellite.name}.")
            else:
                print(f"Order {i} is unacceptable for {satellite.name}.")
        else:
            print(f"Order {i} is unacceptable for {satellite.name}.")