import numpy as np

def great_circle_distance(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance between two points on the Earth's surface.

    Parameters:
    lat1, lon1: Latitude and longitude of the first point in decimal degrees.
    lat2, lon2: Latitude and longitude of the second point in decimal degrees.

    Returns:
    Distance between the two points in kilometers.
    """
    # Convert latitude and longitude from degrees to radians
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)

    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    # Radius of the Earth in kilometers (mean radius)
    R = 6371.0

    # Calculate the distance
    distance = R * c

    return distance