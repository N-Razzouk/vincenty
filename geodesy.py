import arcpy
from geographiclib.geodesic import Geodesic

def p(lat, lon):
    return {'lat': lat, 'lon': lon}

def get_axis(in_srs):
    """
    Return semi-major axis and semi-minor axis of geographic coordinate system

    @type in_sra: str
    @param in_srs: Spatial reference system name or path to prj file, default "WGS 1984"
    @rtype: list
    @return: semi-major axis, semi-minor axis
    """

    if in_srs == "WGS 1984":
        a = 6378137.0
        b = 6356752.314245179
    else:
        srs = arcpy.SpatialReference(in_srs)
        if srs == srs.GCS:
            a = srs.semiMajorAxis
            b = srs.semiMinorAxis
        else:
            print "Invalide spatial reference system"
            return None
    return a, b

def make_buffer(point, m, in_srs="WGS 1984", complex=False, n=50):
    """
    Return a polygon geometry representing a buffer around a point with an approximate radius of [input value]

    Expect a point with lon/lat coordinates and buffer radius in meter
    For simple (not complex) method, buffer is formed by a circle with a radius derived from the average distance (in degree)
    from center point to points a bearing of 0 and 90 degree in distance m.
    Results for this method are most accurate close to the Equator. With further distance to the Equator,
    buffer in degrees would look more like an ellipsoid.
    For complex solution, buffer is formed by n vertexes in equal distance and spacing around center point.


    @type  point: arcpy.Point
    @param point: Point(lon,lat)
    @type  m: float
    @param m: buffer radius (meter)
    @type  in_sra: str
    @param in_srs: Spatial reference system name or path to prj file, default "WGS 1984"
    @type  complex: bool
    @param complex: Return complex buffer geometry
    @type  complex: int
    @param complex: Number of vertexes for complex buffer
    @rtype:   arcpy.geometry
    @return:  Polygon geometry
    """

    x = point.X
    y = point.Y

    try:
        a, b = get_axis(in_srs)
    except TypeError:
        return None

    f = (a-b)/a  # flattening
    geo = Geodesic(a, f)

    srs = arcpy.SpatialReference(in_srs)

    if complex:
        arr = arcpy.Array()
        for i in range(0, 360, 360/int(n)):
            direct = geo.Direct(y, x, i, m)
            arr.append(arcpy.Point(direct['lon2'], direct['lat2']))
        buffer_geom = arcpy.Polygon(arr, srs)
    else:
        x1, y1 = geo.Direct(y, x, 0, m)
        x2, y2 = geo.Direct(y, x, 90, m)

        dy = y1 - y
        dx = x2 - x

        d = (dx + dy)/2

        p_geom = arcpy.PointGeometry(point, srs)
        buffer_geom = p_geom.buffer(d)

    return buffer_geom

def get_area(point, height, width, in_srs="WGS 1984"):
    """Return area of height and width in degree at a given point

    Expect a point with lon/lat coordinates and height and width in degree
    Converts heights and width form degree to meters at the given location Vincenty method
    Calculate area in square meters

    @type  point: arcpy.Point
    @param point: Point(lon,lat)
    @type  height: float
    @param height: height (degree)
    @type  width: float
    @param width: width (degree)
    @type  in_sra: str
    @param in_srs: Spatial reference system name or path to prj file, default "WGS 1984"
    @rtype: float
    @return: Area (m2)
    """

    x = point.X
    y = point.Y

    x1 = x + width
    y1 = y

    x2 = x
    y2 = y + height

    try:
        a, b = get_axis(in_srs)
    except TypeError:
        return None

    f = (a-b)/a  # flattening

    geo = Geodesic(a, f)

    a = geo.Area([p(y, x), p(y1, x1), p(y2, x1), p(y2, x2)])

    #inverse1 = geo.Inverse(y, x, y1, x1)
    #inverse2 = geo.Inverse(y, x, y2, x2)

    #d1 = inverse1['s12']
    #d2 = inverse2['s12']

    #a = d1 * d2

    return a

