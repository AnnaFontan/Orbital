''' astroConstants (from Horizon: 16.09.2022)
This functions returns the astrodynamic-related physical constants * to be updated *

List of identifiers:
    1   Universal gravity constant (G) [km^3/(kg*s^2)]
    2   Astronomical Unit (AU) [km]
    3   Speed of light in the vacuum [km/s]
    4   Standard free fall (the acceleration due to gravity on the Earth's surface at sea level) [m/s^2]

    * Space for other constants *

    Planetary constant [km^3/s^2]: (mu = mass * G)
    10  Sun
    11  Mercury 
    12  Venus
    13  Earth
    14  Mars system
    15  Jupiter system
    16  Saturn system
    17  Uranus system
    18  Neptune system
    19  Moon

    Mean radius [km]:
    20  Sun
    21  Mercury
    22  Venus
    23  Earth
    24  Mars
    25  Jupiter
    26  Saturn
    27  Uranus
    28  Neptune
    29  Moon

    Second zonal harmonics J2 [-]:
    31  Mercury
    32  Venus
    33  Earth
    34  Mars
    35  Jupiter
    36  Saturn
    37  Uranus
    38  Neptune
    39  Moon
'''
def astroConstants(input):

    if (input == 1):
        return 6.67430e-20
    elif (input == 2):
        return 149597870.7
    elif (input == 3):
        return 299792.458
    elif (input == 4):
        return 9.80665

    elif (input == 10):
        return 1.327124400181e11
    elif (input == 11):
        return 22031.868551
    elif (input == 12):
        return 324858.592000
    elif (input == 13):
        return 398600.435507
    elif (input == 14):
        return 42828.375816
    elif (input == 15):
        return 126712764.100000
    elif (input == 16):
        return 37940584.841800
    elif (input == 17):
        return 5794556.400000
    elif (input == 18):
        return 6836527.100580
    elif (input == 19):
        return 4902.800118

    elif (input == 20):
        return 6.9551e5
    elif (input == 21):
        return 2439.4 
    elif (input == 22):
        return 6051.8 
    elif (input == 23):
        return 6371.0084
    elif (input == 24):
        return 3389.50
    elif (input == 25):
        return 69911 
    elif (input == 26):
        return 58232
    elif (input == 27):
        return 25362
    elif (input == 28):
        return 24622
    elif (input == 29):
        return 1737.4

    elif (input == 31):
        return 60e-6
    elif (input == 32):
        return 4.458e-6 
    elif (input == 33):
        return 1.08263e-3
    elif (input == 34):
        return 1.96045e-3
    elif (input == 35):
        return 14.736e-3
    elif (input == 36):
        return 16.298e-3
    elif (input == 37):
        return 3.34343e-3
    elif (input == 38):
        return 3.411e-3
    elif (input == 39):
        return 202.7e-6
    
    else:
        print('Number not yet saved!')