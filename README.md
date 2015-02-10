# Vincenty
Calculate metric distance, buffer and area for lat/lon coordinates using Vincently method.

Python implementation of java script code published [here] (http://www.movable-type.co.uk/scripts/latlong-vincenty.html):

## Dependencies
Needs arcpy

## Usage

```python
import arcpy
import geodesy

point1 = arcpy.Point(5,6)  # default WGS 1984
point2 = arcpy.Point(6,6)  # default WGS 1984

distance =  geodesy.py.get_distance(point1, point2)  # returns distance between points in meter
```

```python
import arcpy
import geodesy

point = arcpy.Point(5,6)  # default WGS 1984
radius = 50000  # radius in meter

buffer_geom =  geodesy.get_buffer(point,radius)  # returns a polygon geometry (circle)
```

```python
import arcpy
import geodesy

point = arcpy.Point(5,6)  # default WGS 1984
radius = 50000  # radius in meter

buffer_geom =  geodesy.make_buffer(point,radius,complex=True)  # returns a polygon geometry (default 50 vertexes)
```



```python
import arcpy
import geodesy

point = arcpy.Point(5,6)  # default WGS 1984
pixel_size =  0.00027777778  # This is the pixel size of a Hanson pixel

pixel_area =  geodesy.get_area(point, pixel_size, pixel_size)  # returns area of a pixel at given location
```





