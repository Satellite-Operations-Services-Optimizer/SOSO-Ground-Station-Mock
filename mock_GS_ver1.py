## Version 1:
## Imports
from skyfield.api import load, EarthSatellite, Topos
from datetime import timedelta, datetime
import matplotlib.pyplot as plt
import json
from czml import czml
import numpy as np
from math import atan, degrees


## Step 1: Ground Station Class
class GroundStation:
    def __init__(self, name, latitude, longitude, height, mask, uplink_rate, downlink_rate):
        self.name = name
        self.latitude = latitude
        self.longitude = longitude
        self.height = height
        self.mask = mask
        self.uplink_rate = uplink_rate
        self.downlink_rate = downlink_rate

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
                access_points[ground_station.name]["Access Timestamp Start"].append(current_time.strftime('%Y-%m-%d %H:%M:%S'))
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
        
# Downlink & Uplink Transfer Times
def simulate_uplink(self, data_size):
    # Simulate uplink data transfer
    transfer_time = data_size / self.uplink_rate
    return transfer_time

def simulate_downlink(self, data_size):
    # Simulate downlink data transfer
    transfer_time = data_size / self.downlink_rate
    return transfer_time

# Sample for Satellites
class EarthSatellite:
    def __init__(self, line1, line2, name, ts):
        self.line1 = line1
        self.line2 = line2
        self.name = name
        self.ts = ts
        self.data_size = 0  # Initialize data size

    def set_data_size(self, size):
        # Set the data size for the satellite (e.g., in bits)
        self.data_size = size

    def simulate_uplink(self, ground_station):
        if self.data_size == 0:
            return 0  # No data to transfer
        return ground_station.simulate_uplink(self.data_size)

    def simulate_downlink(self, ground_station):
        if self.data_size == 0:
            return 0  # No data to transfer
        return ground_station.simulate_downlink(self.data_size)
    
## Example 1: For Downlink & Uplink Transfer
# Create a satellite
satellite = EarthSatellite(line1, line2, name, ts)

# Set the data size to simulate
data_size = 10000000  # For example, 10 megabits

# Simulate uplink and downlink
uplink_time = satellite.simulate_uplink(inuvik_station)
downlink_time = satellite.simulate_downlink(inuvik_station)

print(f"Uplink time: {uplink_time:.2f} seconds")
print(f"Downlink time: {downlink_time:.2f} seconds")