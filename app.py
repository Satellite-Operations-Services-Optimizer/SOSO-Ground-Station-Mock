from flask import Flask, request, jsonify
from fastapi.encoders import jsonable_encoder
from datetime import datetime
from Models import satellite_schedule, ground_station_request

app = Flask(__name__)

@app.route('/image-requests', methods=['POST'])
def add_image_request():
    # Get the request URL and print it
    request_url = request.url
    print(f"Received POST request at URL: {request_url}")
    
    data = request.get_json()
    try:
        # Validate request data
        latitude = float(data['Latitude'])
        longitude = float(data['Longitude'])
        priority = int(data['Priority'])
        image_type = data['ImageType']
        image_start_time = data['ImageStartTime']
        image_end_time = data['ImageEndTime']
        delivery_time = data['DeliveryTime']
        revisit_time = data['RevisitTime']

        # Store the image request
        new_request = {
            'Latitude': latitude,
            'Longitude': longitude,
            'Priority': priority,
            'ImageType': image_type,
            'ImageStartTime': image_start_time,
            'ImageEndTime': image_end_time,
            'DeliveryTime': delivery_time,
            'RevisitTime': revisit_time,
        }
        image_requests.append(new_request)

        return jsonify({"message": "Image request added successfully."}), 200
    except (KeyError, ValueError) as e:
        return jsonify({"error": "Invalid request data."}), 400
    
@app.route('/ActivityRequest/general', methods=['POST'])
def sendActivityRequest():
    requestBody = request.get_json()
    
    return jsonify(requestBody), 200

@app.route('/satellite_schedule', methods=['POST'])
def recieveSatelliteSchedule():
    data = request.get_json()
    
    try:
        # Validate request data
        sat_schedule = satellite_schedule(**data)
        recieved_schedule = jsonable_encoder(sat_schedule)
        print("satellite schedule recieved: ", recieved_schedule )
        return recieved_schedule , 200
    except (KeyError, ValueError) as e:
        print(e)
        return jsonify({"error": "Invalid request data."}), 400
       
        
    

@app.route('/ground_station_schedule', methods=['POST'])
def recieveGroundStationSchedule():
    data = request.get_json()
    
    try:
        # Validate request data
        sat_schedule = ground_station_request(**data)
        recieved_schedule = jsonable_encoder(sat_schedule)
        
        print("ground station schedule recieved: ", recieved_schedule )
        return recieved_schedule , 200
    except (KeyError, ValueError) as e:
        print(e)
        return jsonify({"error": "Invalid request data."}), 400
    
if __name__ == '__main__':
    app.run(debug=True)