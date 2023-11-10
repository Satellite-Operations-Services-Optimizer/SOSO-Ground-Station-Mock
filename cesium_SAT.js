const fs = require('fs');

// Read the JSON file
let rawdata = fs.readFileSync('satellite_data.json');

// Parse the data to a JavaScript object
let data = JSON.parse(rawdata);

// Now you can access the latitude and longitude data in your JavaScript code
var latitudes_soso1 = data.soso1.latitudes;
var longitudes_soso1 = data.soso1.longitudes;

// Create a Cesium Viewer instance
var viewer = new Cesium.Viewer('cesiumContainer');

// Add the trajectories of the satellites to the viewer
for (var i = 0; i < latitudes_soso1.length; i++) {
    viewer.entities.add({
        position : Cesium.Cartesian3.fromDegrees(longitudes_soso1[i], latitudes_soso1[i]),
        point : {
            pixelSize : 5,
            color : Cesium.Color.RED,
            outlineColor : Cesium.Color.WHITE,
            outlineWidth : 2
        }
    });
}

for (var i = 0; i < latitudes_soso2.length; i++) {
    viewer.entities.add({
        position : Cesium.Cartesian3.fromDegrees(longitudes_soso2[i], latitudes_soso2[i]),
        point : {
            pixelSize : 5,
            color : Cesium.Color.YELLOW,
            outlineColor : Cesium.Color.WHITE,
            outlineWidth : 2
        }
    });
}

for (var i = 0; i < latitudes_soso3.length; i++) {
    viewer.entities.add({
        position : Cesium.Cartesian3.fromDegrees(longitudes_soso3[i], latitudes_soso3[i]),
        point : {
            pixelSize : 5,
            color : Cesium.Color.PINK,
            outlineColor : Cesium.Color.WHITE,
            outlineWidth : 2
        }
    });
}

// ... continue this for all your satellites ...
