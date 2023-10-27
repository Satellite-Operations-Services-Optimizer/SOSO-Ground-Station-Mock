from flask import Flask, request, jsonify

app = Flask(__name__)

# Mock Ground Station data (image requests)
image_requests = []

# Endpoint to add a new image request
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

# Endpoint to retrieve a list of image requests
@app.route('/image-requests', methods=['GET'])
def get_image_requests():
    return jsonify({"image_requests": image_requests})

if __name__ == '__main__':
    app.run(debug=True)
