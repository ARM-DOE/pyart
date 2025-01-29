"""
Data and functions common to all types of NEXRAD files.

"""
# The functions in this module are intended to be used in other
# nexrad related modules. The functions are not and should not be
# exported into the pyart.io namespace.

import pandas as pd

def get_nexrad_location(station):
    """
    Return the latitude, longitude and altitude of a NEXRAD station.

    Parameters
    ----------
    station : str
        Four letter NEXRAD station ICAO name.

    Returns
    -------
    lat, lon, alt : float
        Latitude (in degrees), longitude (in degrees), and altitude
        (in meters above mean sea level) of the NEXRAD station.

    """
    loc = NEXRAD_LOCATIONS[station.upper()]

    # Convert from feet to meters for elevation units
    loc["elev"] = loc["elev"] * 0.3048

    return loc["lat"], loc["lon"], loc["elev"]


# Locations of NEXRAD locations can be retrieved from NOAA's
# Historical Observing Metadata Repository (HOMR).
# The NEXRAD station location information is available
# at the following link:
# https://www.ncei.noaa.gov/access/homr/

# Read the nexrad-stations.txt file directly from the URL
nexrad_station_file = "https://www.ncei.noaa.gov/access/homr/file/nexrad-stations.txt;jsessionid=EC5C4FA497A46EFAF53FD6068F870D3E"

# Define the column specifications
colspecs = [
    (0, 8),  # NCDCID
    (9, 13),  # ICAO
    (14, 19),  # WBAN
    (20, 50),  # NAME
    (51, 71),  # COUNTRY
    (72, 74),  # ST
    (75, 105),  # COUNTY
    (106, 115),  # LAT
    (116, 126),  # LON
    (127, 133),  # ELEV
    (134, 139),  # UTC
    (140, 190)   # STNTYPE
]

# Define the column names
colnames = [
    "NCDCID",
    "ICAO",
    "WBAN",
    "NAME", 
    "COUNTRY", 
    "ST", 
    "COUNTY", 
    "LAT", 
    "LON", 
    "ELEV", 
    "UTC", 
    "STNTYPE",
]

# Read the fixed-width file
df_nexrad_loc = pd.read_fwf(
    nexrad_station_file, colspecs=colspecs, names=colnames, skiprows=2
)


# Create the dictionary
NEXRAD_LOCATIONS = {
    row["ICAO"]: {
        "lat": float(row["LAT"]),
        "lon": float(row["LON"]),
        "elev": int(row["ELEV"]),
    }
    for _, row in df_nexrad_loc.iterrows()
}
