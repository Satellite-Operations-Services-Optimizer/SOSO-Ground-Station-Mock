# SOSO-Ground-Station-Mock

Repo for ground stations.

# How to run the program

## Basic setup

In order to run the program, make sure you have python and flask installed. <br>
To check whether flask is installed, run the following command, and look for Flask in the list.

```
pip list
```

## Running the program

Run the command below to run the program.

```
flask run
```

# Current Endpoints

List of current endpoints available

## Image Requests

```
/image-requests
```

Returns a confirmation of reception, along with a copy of the data received.

## Activity Requests

```
/ActivityRequest/general
```

Returns a confirmation of reception, along with a copy of the code.
