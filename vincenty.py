#    Vincenty's formulae calculate the distance between two points (inverse solution)
#    and the destination with given distance & bearing from start point (direct solution) on the surface of a spheroid
#    Its accuracy is within 0.5 mm distance, 0.000015" bearing, on the ellipsoid being used.
#    Methods were first published in 1975: http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
#    This is a Python implementation of a JavaScript source code by Chris Veness (Creative Commons Attribution license)
#    http://www.movable-type.co.uk/scripts/latlong-vincenty.html

import math

def inverse(x1, y1, x2, y2, a=6378137.0, b=6356752.314245179, bearing=False):
    """
    Return distance/bearing between two points (Vincenty's inverse solution)

    @type  x1: float
    @param x1: X (lon) Point 1
    @type  y1: float
    @param y1: Y (lat) Point 1
    @type  x2: float
    @param x2: X (lon) Point 2
    @type  y2: float
    @param y2: Y (lat) Point 2
    @type  a: float
    @param a: Major-axis ellipsoid
    @type  b: float
    @param b: Minor-axis ellipsoid
    @type  bearing: bool
    @param bearing: Return bearing
    @rtype: float
    @return: distance (m) [, forward azimuth (deg 0-360), reverse Azimuth (deg 0-360)]
    """

    if -180 >= x1 >= 180 or -180 >= x2 >= 180 or -90 >= y1 >= 90 or -90 >= y2 >= 90:
        print "Wrong input coordinates"
        return None
            
    f = (a-b)/a  # flattening
    
    lat1 = math.radians(y1)
    lon1 = math.radians(x1)

    lat2 = math.radians(y2)
    lon2 = math.radians(x2)

    L = lon2 - lon1

    tanU1 = (1-f) * math.tan(lat1)
    cosU1 = 1 / math.sqrt((1 + tanU1*tanU1))
    sinU1 = tanU1 * cosU1

    tanU2 = (1-f) * math.tan(lat2)
    cosU2 = 1 / math.sqrt((1 + tanU2*tanU2))
    sinU2 = tanU2 * cosU2

    l = L
       
    for i in range(100):
        sinl = math.sin(l)
        cosl = math.cos(l)
        sinSqs = (cosU2*sinl) * (cosU2*sinl) + (cosU1*sinU2-sinU1*cosU2*cosl) * (cosU1*sinU2-sinU1*cosU2*cosl)
        sins = math.sqrt(sinSqs)
        if (sins == 0):
            return 0  # co-incident points
        coss = sinU1*sinU2 + cosU1*cosU2*cosl
        s = math.atan2(sins, coss)
        sina = cosU1 * cosU2 * sinl / sins
        cosSqa = 1 - sina*sina
        cos2sM = coss - 2*sinU1*sinU2/cosSqa
        if math.isnan(cos2sM):
            cos2sM = 0  # equatorial line: cosSqa=0
        C = f/16*cosSqa*(4+f*(4-3*cosSqa))
        l_ = l
        l = L + (1-C) * f * sina * (s + C*sins*(cos2sM+C*coss*(-1+2*cos2sM*cos2sM)))
        if abs(l - l_) < 1e-12:
            break

    if i == 99:
        print 'Formula failed to converge'

    uSq = cosSqa * (a*a - b*b) / (b*b)
    A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
    B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)))
    ds = B*sins*(cos2sM+B/4*(coss*(-1+2*cos2sM*cos2sM)-B/6*cos2sM*(-3+4*sins*sins)*(-3+4*cos2sM*cos2sM)))

    d = b*A*(s-ds)

    fwdAz = math.degrees(math.atan2(cosU2*sinl,  cosU1*sinU2-sinU1*cosU2*cosl))
    revAz = math.degrees(math.atan2(cosU1*sinl, -sinU1*cosU2+cosU1*sinU2*cosl))

    if bearing:
        return d, fwdAz, revAz
    else:
        return d

def direct(x1, y1, d, a1, a=6378137.0, b=6356752.314245179, bearing=False):
    """"
    Return destination given distance & bearing from start point (Vincenty's direct solution)

    @type  x1: float
    @param x1: X (lon) Point 1
    @type  y1: float
    @param y1: Y (lat) Point 1
    @type  d: float
    @param d: Distance (m)
    @type  a1: float
    @param a1: bearing (deg 0-360)
    @type  a: float
    @param a: Major-axis ellipsoid
    @type  b: float
    @param b: Minor-axis ellipsoid
    @type  bearing: bool
    @param bearing: Return bearing
    @rtype: list
    @return: X (lon), Y (lat) [, reverse azimuth (deg 0-360)]

    """

    if -180 >= x1 >= 180 or -90 >= y1 >= 90:
        print "Wrong input coordinates"
        return None

    f = (a-b)/a  # flattening

    lat1 = math.radians(y1)
    lon1 = math.radians(x1)

    sina1 = math.sin(math.radians(a1))
    cosa1 = math.cos(math.radians(a1))

    tanU1 = (1-f) * math.tan(lat1)
    cosU1 = 1 / math.sqrt((1 + tanU1*tanU1))
    sinU1 = tanU1 * cosU1
    s1 = math.atan2(tanU1, cosa1)
    sina = cosU1 * sina1
    cosSqa = 1 - sina*sina
    uSq = cosSqa * (a*a - b*b) / (b*b)
    A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
    B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)))

    s = d / (b*A)

    while True:
        cos2sM = math.cos(2*s1 + s)
        sins = math.sin(s)
        coss = math.cos(s)
        ds = B*sins*(cos2sM+B/4*(coss*(-1+2*cos2sM*cos2sM)-B/6*cos2sM*(-3+4*sins*sins)*(-3+4*cos2sM*cos2sM)))
        s_ = s
        s = d / (b*A) + ds
        if abs(s-s_) > 1e-12:
            break

    tmp = sinU1*sins - cosU1*coss*cosa1
    lat2 = math.atan2(sinU1*coss + cosU1*sins*cosa1, (1-f)*math.sqrt(sina*sina + tmp*tmp))
    l = math.atan2(sins*sina1, cosU1*coss - sinU1*sins*cosa1)
    C = f/16*cosSqa*(4+f*(4-3*cosSqa))
    L = l - (1-C) * f * sina * (s + C*sins*(cos2sM+C*coss*(-1+2*cos2sM*cos2sM)))
    lon2 = (lon1 + L + 3*math.pi)%(2*math.pi) - math.pi  # normalise to -180...+180

    revAz = math.degrees(math.atan2(sina, -tmp))

    y2 = math.degrees(lat2)
    x2 = math.degrees(lon2)

    if bearing:
        return x2, y2, revAz
    else:
        return x2, y2