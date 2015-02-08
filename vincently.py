import arcpy
import math
import sys

def get_distance(point1, point2, in_srs="WGS 1984", bearing=False):
    """Return distance, initial, and final bearing between two points (Vincenty method)

    Vincenty's solution for the distance between points on an ellipsoidal earth model
    is accurate to within 0.5 mm distance, 0.000015" bearing, on the ellipsoid being used.
    This code is a python implementation of code found here
    http://www.movable-type.co.uk/scripts/latlong-vincenty.html

    @type  point1: arcpy.Point
    @param point1: Point(lon,lat)
    @type  point2: arcpy.Point
    @param point2 Point(lon,lat)
    @type  in_sra: str
    @param in_srs: Spatial reference system name or path to prj file, default "WGS 1984"
    @rtype:   number
    @return:  distance (m)
    """

    if in_srs == "WGS 1984":
        a = 6378137.0  # major semi-axis of the ellipsoid
        b = 6356752.314245179  # minor semi-axis of the ellipsoid
    else:
        try:
            srs = arcpy.SpatialReference(in_srs)
            a = srs.semiMajorAxis
            b = srs.semiMinorAxis
        except:
            print "Spatial reference system in not well defined"
            sys.exit(2)
            
    f = (a-b)/a  # flattening
    
    lat1 = math.radians(point1.Y)
    lon1 = math.radians(point1.X)

    lat2 = math.radians(point2.Y)
    lon2 = math.radians(point2.X)

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

    fwdAz = math.atan2(cosU2*sinl,  cosU1*sinU2-sinU1*cosU2*cosl)
    revAz = math.atan2(cosU1*sinl, -sinU1*cosU2+cosU1*sinU2*cosl)

    if bearing:
        return d, fwdAz, revAz
    else:
        return d

def make_buffer(point, m, in_srs="WGS 1984"):
    """Return a polygon geometry representing a buffer around a point with an approximate radius of [input value]

    Expect a point with lat/lon coordinates and buffer radius in meter
    Buffer radius is converted into degree using Vincenty method in north-south and east-west direction
    Final buffer is the mean radius (in degree) of north-south and east-west radius
    It represents the input radius best, closest to the equator.

    @type  point: arcpy.Point
    @param point: Point(lon,lat)
    @type  m: number
    @param m2 buffer radius (meter)
    @type  in_sra: str
    @param in_srs: Spatial reference system name or path to prj file, default "WGS 1984"
    @rtype:   arcpy.geometry
    @return:  Polygon geometry
    """

    diff = 0.45  # proxy value, about 50000 m equator

    lat = point.Y
    lon = point.X

    point_v = arcpy.Point(lon, lat + diff)
    point_h = arcpy.Point(lon + diff, lat)

    d_v = get_distance(point, point_v, in_srs, False)
    d_h = get_distance(point, point_h, in_srs, False)

    deg_v = diff * m / d_v
    deg_h = diff * m / d_h

    deg = (deg_v + deg_h)/2

    try:
        srs = arcpy.SpatialReference(in_srs)
    except:
        print "Spatial reference system in not well defined"
        sys.exit(2)
    p_geom = arcpy.PointGeometry(point, srs)
    buffer_geom = p_geom.buffer(deg)
    
    return buffer_geom

def get_area(point, height, width, in_srs="WGS 1984"):
    """Return area of height and width in degree at a given point

    Expect a point with lat/lon coordinates and height and width in degree
    Converts heights and width form degree to meters at the given location Vincenty method
    Calculate area in square meters

    @type  point: arcpy.Point
    @param point: Point(lon,lat)
    @type  height: number
    @param height: height (degree)
    @type  width: number
    @param width: width (degree)
    @type  in_sra: str
    @param in_srs: Spatial reference system name or path to prj file, default "WGS 1984"
    @rtype:   number
    @return:  Area in square meter
    """

    lat = point.Y
    lon = point.X

    point_v = arcpy.Point(lon, lat + width)
    point_h = arcpy.Point(lon + height, lat)
    
    d_v = get_distance(point, point_v, in_srs, False)
    d_h = get_distance(point, point_h, in_srs, False)

    a = d_v * d_h

    return a

