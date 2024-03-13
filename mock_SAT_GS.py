# ------------------------------------------------------- #
# Combination of mock_SAT.py and mock_GS.py to demonstrate flow of program for mapping.
# Testing Procedures for SAT:
# - Latitude and Longitude of each satellite. (Method: Hand calculation)
# - Altitude of satellite. (Method: Hand calculation)
# - FOV of satellite based on formula used from Earth coverage data. (Method: Hand calculation + 2D projection) 
# - Image Validation based on CSA requirements (Method: 2D projection)
# Testing Procedures for GS:
# - Access points between satellites. (Method: STK simulation)
# - Preconfiguration of GS for delta-time. (Method: STK simulation)
# ------------------------------------------------------- #

## Imports
import json
import os
import numpy as np
import matplotlib.pyplot as plt
from math import atan, degrees
from skyfield.api import load, EarthSatellite, Topos
from datetime import timedelta, datetime
from skyfield.sgp4lib import EarthSatellite
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

# Power Management Example
P_sunlit = 500 # in Watts during Sunlight
# 200-800 Watts for research sat.
# 1000-1500 Watts for commercial sat.
P_eclipse = P_sunlit * 0.4 # in Watts during Eclipse (assuming 40% of power is used)

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

## Step 4: Get live data of the satellites' position (i.e. x, y, z coordinates, latitude, longitude)
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

        print(f"SOSO-{i + 1} at Current Time:")
        print(f"  Position: {position_current} km")
        print(f"  Latitude: {latitude_current} degrees")
        print(f"  Longitude: {longitude_current} degrees")
        print(f"  Altitude: {altitude_current} km")
        # print(f"  Field of View: {fov_current} degrees")
        print(f"  Time: {current_time_skyfield.utc_strftime('%Y %b %d %H:%M:%S')}")

        current_time_skyfield = ts.utc(current_time_skyfield.utc_datetime() + timedelta(minutes=1)) # Print all variables every minute from start and end times.

for i, satellite in enumerate(satellites):  # Loop over the satellites in the list.
    current_time_skyfield = start_time_skyfield
    while current_time_skyfield.tt < end_time_skyfield.tt:  # Compare Julian dates
        sat_pos_current = satellite.at(current_time_skyfield).position.km  # Get the position of the satellite relative to Earth at the current time

        is_sunlit_current = satellite.at(current_time_skyfield).is_sunlit(eph) # Check if satellite is sunlit at current time

        if is_sunlit_current:
            print(f"SOSO-{i + 1} at Current Time: {is_sunlit_current}")
            print(f"  The satellite is in sunlight. Power is unrestricted = {P_sunlit} W.")
            print(f"  Time: {current_time_skyfield.utc_strftime('%Y %b %d %H:%M:%S')}")
        else:
            print(f"SOSO-{i + 1} at Current Time: {is_sunlit_current}")
            print(f"  The satellite is in eclipse. All activities must use less power than a predefined limit = {P_eclipse} W.")
            print(f"  Time: {current_time_skyfield.utc_strftime('%Y %b %d %H:%M:%S')}")

        current_time_skyfield = ts.utc(current_time_skyfield.utc_datetime() + timedelta(minutes=1))

## Step 5: Load the Image Orders & Validation Methodology

# Spotlight image = 10 km height x 10 km width square area, 120 second transfer time, 512 MB.
# Medium image = 40 km height x 20 km width square area, 45 second transfer time, 256 MB.
# Low image = 40 km height x 20 km width square area, 20 second transfer time, 128 MB.

# Given: Satellite FOV -> Viewing Angle = 30 degrees, Full View = 60 degrees.
# If FOV square area of 577 km height x 577 km width propagated on ground track doesn't cover square area of image in image order list, then unacceptable.
# Variables to consider in image order: Latitude, Longitude, Image Start Time, Image End Time, Revisit Time (True or False).
# Image order square area should be within the satellite FOV square area throughout the transfer time during between the start and end time of the image.

# False Image Orders: Order 37 for, 

# # Define the image types and their dimensions
# image_types = {
#     "Spotlight": [10, 10],
#     "Medium": [40, 20],
#     "Low": [40, 20]
# }

# # Iterate over the image orders
# for i in range(1, 51):
#     with open(f'SampleOrders/Order_{i}.json') as f:
#         order = json.load(f)
#     image_type = order["ImageType"]
#     lat, lon = order["Latitude"], order["Longitude"]
#     start_time = datetime.strptime(order["ImageStartTime"], "%Y-%m-%dT%H:%M:%S")
#     end_time = datetime.strptime(order["ImageEndTime"], "%Y-%m-%dT%H:%M:%S")

#     # Calculate the square area for the image type
#     image_area = image_types[image_type]

#     # Iterate over the satellites
#     for satellite in satellites:
#         # Calculate the satellite's FOV square area
#         satellite_area = [577, 577]  # Replace with actual calculation

#         # Note 2: Dependent on altitude. 
        
#         # Check if the image area falls within the satellite's FOV
#         if image_area[0] <= satellite_area[0] and image_area[1] <= satellite_area[1]:
#             # Calculate the transfer time
#             transfer_time = end_time - start_time
#             # Note 1: Change for 120, 45, 20 sec for transfer time.

#             # Check if the satellite stays within the FOV during the transfer time
#             t = ts.utc(start_time.year, start_time.month, start_time.day, start_time.hour, start_time.minute, start_time.second)
#             geocentric = satellite.at(t)
#             subpoint = geocentric.subpoint()
#             start_lat, start_lon = subpoint.latitude.degrees, subpoint.longitude.degrees

#             t = ts.utc(end_time.year, end_time.month, end_time.day, end_time.hour, end_time.minute, end_time.second)
#             geocentric = satellite.at(t)
#             subpoint = geocentric.subpoint()
#             end_lat, end_lon = subpoint.latitude.degrees, subpoint.longitude.degrees

#             # If the satellite stays within the FOV, the image is acceptable
#             if start_lat <= lat <= end_lat and start_lon <= lon <= end_lon:
#                 print(f"Image order {i} is unacceptable for satellite {satellite.name}")
#             else:
#                 print(f"Image order {i} is acceptable for satellite {satellite.name}")
#         else:
#             print(f"Image order {i} is acceptable for satellite {satellite.name}")

#         # Plot the FOV and image area
#         fig, ax = plt.subplots()
#         ax.add_patch(plt.Rectangle((lon - image_area[1] / 2, lat - image_area[0] / 2), image_area[1], image_area[0], fill=None, edgecolor='r'))
#         ax.add_patch(plt.Rectangle((start_lon - satellite_area[1] / 2, start_lat - satellite_area[0] / 2), satellite_area[1], satellite_area[0], fill=None, edgecolor='b'))
#         ax.set_xlim([min(lon - image_area[1] / 2, start_lon - satellite_area[1] / 2) - 10, max(lon + image_area[1] / 2, start_lon + satellite_area[1] / 2) + 10])
#         ax.set_ylim([min(lat - image_area[0] / 2, start_lat - satellite_area[0] / 2) - 10, max(lat + image_area[0] / 2, start_lat + satellite_area[0] / 2) + 10])
#         plt.show()
        
#         # Note 3: add the lat and lon for the time for when image is acceptable.

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
        semi_major_axis_km = satellite.model.a * 6378.137  # Get the semi-major axis in kilometers
        altitude_current = semi_major_axis_km - 6378.137  # The altitude is the semi-major axis minus the Earth's radius
        fov_current = degrees(2 * atan(12742 / (2 * (altitude_current + 6371))))  # Calculate FOV

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
            
###################################################################################################################################################################################################
################################ GROUND STATION CODE
###################################################################################################################################################################################################

## Step 4: Ground Station Class
class GroundStation:
    def __init__(self, name, latitude, longitude, height, mask, uplink_rate, downlink_rate):
        self.name = name
        self.latitude = latitude
        self.longitude = longitude
        self.height = height
        self.mask = mask
        self.uplink_rate = uplink_rate
        self.downlink_rate = downlink_rate

# Downlink & Uplink Transfer Times
def simulate_uplink(self, data_size):
    # Simulate uplink data transfer
    # Memory Scrubs, Maintainence, Outage Requests
    transfer_time = data_size / self.uplink_rate
    return transfer_time

def simulate_downlink(self, data_size):
    # Simulate downlink data transfer
    # Image Orders
    transfer_time = data_size / self.downlink_rate
    return transfer_time

# Create ground stations
inuvik_station = GroundStation("Inuvik", 68.3195, -133.549, 102.5, 5, 40000, 100000000)
prince_albert_station = GroundStation("Prince Albert", 53.2124, -105.934, 490.3, 5, 40000, 100000000)
gatineau_station = GroundStation("Gatineau", 45.5846, -75.8083, 240.1, 5, 40000, 100000000)

# Step 5: Loop over the ground stations (Preconfiguration, Elevation constraint, and Access Points)
for ground_station in [inuvik_station, prince_albert_station, gatineau_station]:
    ground_station_topos = Topos(ground_station.latitude, ground_station.longitude)

    # Get the position of the satellite relative to the ground station at the start and end times
    relative_position_start = (satellite - ground_station_topos).at(start_time_skyfield)
    relative_position_end = (satellite - ground_station_topos).at(end_time_skyfield)

    # Get the altitude (elevation) of the satellite from the perspective of the ground station
    elevation_angle_start = relative_position_start.altaz()[0]
    elevation_angle_end = relative_position_end.altaz()[0]

    # Check if the satellite is above the ground station's mask angle at both start and end times
    if elevation_angle_start.degrees > ground_station.mask and elevation_angle_end.degrees > ground_station.mask:
        print(f"The satellite is visible from {ground_station.name} between {start_time} and {end_time} with elevations of {elevation_angle_start.degrees} to {elevation_angle_end.degrees} degrees.")
        # Preconfigure the ground station for 5 minutes before communication
        t0 = start_time
        t1 = t0 + timedelta(seconds=5)  # 5 minutes
        print(f"Preconfiguring {ground_station.name} for communication with the satellite...")
        print(f"Start time: {t0}")
        print(f"End time: {t1}")
    else:
        print(f"The satellite is not visible from {ground_station.name} between {start_time} and {end_time}.")

# Create a list to store access points for each ground station
access_points = {}

# Define a function to check if the satellite is in contact with the ground station
def is_in_contact(satellite, ground_station, time):
    ground_station_topos = Topos(ground_station.latitude, ground_station.longitude)
    relative_position = (satellite - ground_station_topos).at(time)
    elevation_angle = relative_position.altaz()[0]
    return elevation_angle.degrees > ground_station.mask

# Iterate through the time interval minute by minute
current_time = start_time
while current_time < end_time:
    current_time_skyfield = ts.utc(current_time.year, current_time.month, current_time.day, current_time.hour, current_time.minute, current_time.second)
    
    for ground_station in [inuvik_station, prince_albert_station, gatineau_station]:
        if is_in_contact(satellite, ground_station, current_time_skyfield):
            if ground_station.name not in access_points:
                access_points[ground_station.name] = {
                    "Access Timestamp Start": [],
                    "Access Timestamp End": []
                }
            access_points[ground_station.name]["Access Timestamp Start"].append(current_time.strftime('%Y-%m-%d %H:%M:%S'))
            while current_time < end_time:
                current_time += timedelta(minutes=1)
                current_time_skyfield = ts.utc(current_time.year, current_time.month, current_time.day, current_time.hour, current_time.minute, current_time.second)
                if not is_in_contact(satellite, ground_station, current_time_skyfield):
                    break
            access_points[ground_station.name]["Access Timestamp End"].append(current_time.strftime('%Y-%m-%d %H:%M:%S'))
        else:
            current_time += timedelta(minutes=1)

# Print access points for each ground station
for ground_station, timestamps in access_points.items():
    print(f"{ground_station}:")
    print(f"Access Timestamp Start: {timestamps['Access Timestamp Start']}")
    print(f"Access Timestamp End: {timestamps['Access Timestamp End']}")

# Create a dictionary to keep track of satellite access
satellite_access = {ground_station.name: None for ground_station in [inuvik_station, prince_albert_station, gatineau_station]}

# Iterate through the time interval minute by minute
current_time = start_time
while current_time < end_time:
    current_time_skyfield = ts.utc(current_time.year, current_time.month, current_time.day, current_time.hour, current_time.minute, current_time.second)
    
    for ground_station in [inuvik_station, prince_albert_station, gatineau_station]:
        if is_in_contact(satellite, ground_station, current_time_skyfield):
            # Check if the ground station is not in use by any other satellite at the current time
            if satellite_access[ground_station.name] is None:
                if ground_station.name not in access_points:
                    access_points[ground_station.name] = {
                        "Access Timestamp Start": [],
                        "Access Timestamp End": [],
                        "Satellite Name": []  # Initialize the 'Satellite Name' list
                    }
                else:
                    # Initialize the 'Satellite Name' list if it doesn't exist
                    if "Satellite Name" not in access_points[ground_station.name]:
                        access_points[ground_station.name]["Satellite Name"] = []
                while current_time < end_time:
                    current_time += timedelta(minutes=1)
                    current_time_skyfield = ts.utc(current_time.year, current_time.month, current_time.day, current_time.hour, current_time.minute, current_time.second)
                    if not is_in_contact(satellite, ground_station, current_time_skyfield):
                        break
                access_points[ground_station.name]["Access Timestamp End"].append(current_time.strftime('%Y-%m-%d %H:%M:%S'))
                access_points[ground_station.name]["Satellite Name"].append(satellite.name)  # Store the satellite name accessing the ground station
                satellite_access[ground_station.name] = satellite.name  # Mark the ground station as in use by the current satellite
            else:
                current_time += timedelta(minutes=1)
        else:
            current_time += timedelta(minutes=1)
            # Reset ground station access if the satellite moves out of range
            for gs in satellite_access.keys():
                if satellite_access[gs] == satellite.name:
                    satellite_access[gs] = None

# Print access points for each ground station
for ground_station, data in access_points.items():
    print(f"{ground_station}:")
    for i in range(len(data["Access Timestamp Start"])):
        print(f"Access Timestamp Start: {data['Access Timestamp Start'][i]}")
        print(f"Access Timestamp End: {data['Access Timestamp End'][i]}")
        print(f"Satellite Name: {data['Satellite Name'][i]}")
        print()

###################################################################################################################################################################################################
################################ FUNCTIONS
###################################################################################################################################################################################################

# Function 1: Read TLE files for SOSO constellation
def load_satellites():
    ts = load.timescale()
    satellites = []

    for i in range(1, 6):
        try:
            with open(f'SOSO-{i}_TLE.json') as f:
                data = json.load(f)
                name = data['name']
                line1 = data['line1']
                line2 = data['line2']
            satellite = EarthSatellite(line1, line2, name, ts)
            satellites.append(satellite)
        except IndexError:
            print(f"Error: TLE file for satellite {i} is not formatted correctly.")

    # return ts, satellites
    pass

# Function 2: Ephemeris Data, Power Draw, and User-defined Start and End Times for Schedule
def ephemeris_power_schedule_start_end():
    ## Step 2: (Maintenance) Is the satellite in eclipse or in sunlight?
    eph = load('de421.bsp')  # Load the JPL ephemeris DE421
    sun = eph['sun']  # Get the 'sun' object from the ephemeris
    earth = eph['earth']  # Get the 'earth' object from the ephemeris

    # Power Management Example
    P_sunlit = 500 # in Watts during Sunlight
    # 200-800 Watts for research sat.
    # 1000-1500 Watts for commercial sat.
    P_eclipse = P_sunlit * 0.4 # in Watts during Eclipse (assuming 40% of power is used)

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
    pass

# Function 3: Location of the Satellites
def location_satellite():
    ## Step 4: Get live data of the satellites' position (i.e. x, y, z coordinates, latitude, longitude)
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

            print(f"SOSO-{i + 1} at Current Time:")
            print(f"  Position: {position_current} km")
            print(f"  Latitude: {latitude_current} degrees")
            print(f"  Longitude: {longitude_current} degrees")
            print(f"  Altitude: {altitude_current} km")
            # print(f"  Field of View: {fov_current} degrees")
            print(f"  Time: {current_time_skyfield.utc_strftime('%Y %b %d %H:%M:%S')}")

            current_time_skyfield = ts.utc(current_time_skyfield.utc_datetime() + timedelta(minutes=1)) # Print all variables every minute from start and end times.
    pass

# Function 4: Power Management
def eclipse_sunlight():
    for i, satellite in enumerate(satellites):  # Loop over the satellites in the list.
        current_time_skyfield = start_time_skyfield
        while current_time_skyfield.tt < end_time_skyfield.tt:  # Compare Julian dates
            sat_pos_current = satellite.at(current_time_skyfield).position.km  # Get the position of the satellite relative to Earth at the current time

            is_sunlit_current = satellite.at(current_time_skyfield).is_sunlit(eph) # Check if satellite is sunlit at current time

            if is_sunlit_current:
                print(f"SOSO-{i + 1} at Current Time: {is_sunlit_current}")
                print(f"  The satellite is in sunlight. Power is unrestricted = {P_sunlit} W.")
                print(f"  Time: {current_time_skyfield.utc_strftime('%Y %b %d %H:%M:%S')}")
            else:
                print(f"SOSO-{i + 1} at Current Time: {is_sunlit_current}")
                print(f"  The satellite is in eclipse. All activities must use less power than a predefined limit = {P_eclipse} W.")
                print(f"  Time: {current_time_skyfield.utc_strftime('%Y %b %d %H:%M:%S')}")

            current_time_skyfield = ts.utc(current_time_skyfield.utc_datetime() + timedelta(minutes=1))
    pass

# Function 5: Image Validation
def image_validation():
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
            semi_major_axis_km = satellite.model.a * 6378.137  # Get the semi-major axis in kilometers
            altitude_current = semi_major_axis_km - 6378.137  # The altitude is the semi-major axis minus the Earth's radius
            fov_current = degrees(2 * atan(12742 / (2 * (altitude_current + 6371))))  # Calculate FOV

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

                # Check if the satellite stays within the FOV during the transfer time
                transfer_end_time = start_time + timedelta(seconds=image_type["transfer_time"])
                if start_time <= transfer_end_time <= end_time:
                    if start_position[0] <= end_position[0] and start_position[1] <= end_position[1]:
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
    pass

###################################################################################################################################################################################################
################################ FRONT-END
###################################################################################################################################################################################################

# Convert lists to numpy arrays after collecting all data
latitudes_soso1 = np.array(latitudes[0])
longitudes_soso1 = np.array(longitudes[0])

latitudes_soso2 = np.array(latitudes[1])
longitudes_soso2 = np.array(longitudes[1])

latitudes_soso3 = np.array(latitudes[2])
longitudes_soso3 = np.array(longitudes[2])

latitudes_soso4 = np.array(latitudes[3])
longitudes_soso4 = np.array(longitudes[3])

latitudes_soso5 = np.array(latitudes[4])
longitudes_soso5 = np.array(longitudes[4])

## [FRONT-END] Cesium Front-end (Satellite View) -> Parsing data into json file for js (cesium_SAT.js)
# Create a dictionary to hold your data
data = {
    'soso1': {
        'latitudes': latitudes_soso1.tolist(),
        'longitudes': longitudes_soso1.tolist()
    },
    'soso2': {
        'latitudes': latitudes_soso2.tolist(),
        'longitudes': longitudes_soso2.tolist()
    },
    'soso3': {
        'latitudes': latitudes_soso3.tolist(),
        'longitudes': longitudes_soso3.tolist()
    },
    'soso4': {
        'latitudes': latitudes_soso4.tolist(),
        'longitudes': longitudes_soso4.tolist()
    },
    'soso5': {
        'latitudes': latitudes_soso5.tolist(),
        'longitudes': longitudes_soso5.tolist()
    },
    # ... continue this for all your satellites ...
}

# Write the data to a JSON file
with open('satellite_data.json', 'w') as f:
    json.dump(data, f)
    
# Create a new figure
fig = plt.figure(1)

# Add a subplot with a projection of 'mollweide'
ax = fig.add_subplot(111, projection='mollweide')

# Convert degrees to radians as required by matplotlib for mollweide projection
latitudes_soso1_rad = np.deg2rad(latitudes_soso1)
longitudes_soso1_rad = np.deg2rad(longitudes_soso1)

latitudes_soso2_rad = np.deg2rad(latitudes_soso2)
longitudes_soso2_rad = np.deg2rad(longitudes_soso2)

latitudes_soso3_rad = np.deg2rad(latitudes_soso3)
longitudes_soso3_rad = np.deg2rad(longitudes_soso3)

latitudes_soso4_rad = np.deg2rad(latitudes_soso4)
longitudes_soso4_rad = np.deg2rad(longitudes_soso4)

latitudes_soso5_rad = np.deg2rad(latitudes_soso5)
longitudes_soso5_rad = np.deg2rad(longitudes_soso5)

# Plot the trajectories of the satellites
ax.plot(longitudes_soso1_rad, latitudes_soso1_rad, color='red', label='SOSO-1')
ax.plot(longitudes_soso2_rad, latitudes_soso2_rad, color='yellow', label='SOSO-2')
ax.plot(longitudes_soso3_rad, latitudes_soso3_rad, color='pink', label='SOSO-3')
ax.plot(longitudes_soso4_rad, latitudes_soso4_rad, color='blue', label='SOSO-4')
ax.plot(longitudes_soso5_rad, latitudes_soso5_rad, color='green', label='SOSO-5')

# Add a legend
ax.legend()

# Show the plot
plt.show()

plt.figure(2)
plt.title('Top Down View of Satellite FOV and Image Order')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()