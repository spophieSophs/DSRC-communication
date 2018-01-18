import json, os, geojsonio
import geopandas as gpd
from math import sqrt, log10, sin, cos, pi
from shapely.geometry import shape, MultiLineString, MultiPoint, Point, LineString, box
from geopy.distance import vincenty
import shapely
from shapely.ops import transform
from functools import partial
import pyproj
import scipy.stats as stat
from rtree import index
"""This script includes all the necessary functions of building geometry and receiving power model.
   NOTE-1: update the radio propagation model, fixed the bugs of double counting the transmitting power."""

def caldistance(Tx, Rx):
    """This function calculates distance between two nodes, x-y, in meters."""
    distance = sqrt((float(Tx[0]) - float(Rx[0]))**2 + (float(Tx[1])-float(Rx[1]))**2)
    return distance

def checkdistance(dist, dist_th):
    """This function checks if two nodes are in the transmission range, with input of distance and relevant threshold."""
    return bool (dist < dist_th)

def meetsROIareaxml(position, polygon):
    """This function check if the vehicle position(x, y) is in the ROI area, building geometry-lust.poly.xml"""
    """
        Input:
        position:       the vehicle position in sumo network, (x, y)
        polygon:        the region of interest area, the coordinates of left-bottom point and right-up point[(x, y), (x, y)]
    Output:
        bool:           whether the position is in the ROI area"""
    lebt = polygon[0]
    riup = polygon[1]

    if lebt[0] < position[0] < riup[0] and lebt[1] < position[1] < riup [1]:
        return True
    else:
        return False 

def meetsROIarea(position, polygon):
    """This function checks if the vehicle position is in the ROI area"""
    """
    Input:
        position:       the vehicle position in sumo network, (lat, lon)
        polygon:        the region of interest area, the coordinates of left-bottom point and right-up point[(lat, lon), (lat, lon)]
    Output:
        bool:           whether the position is in the ROI area
    """
    lebt = polygon[0]
    riup = polygon[1]
    if lebt[0] < position[0] < riup[0] and lebt[1] < position[1] < riup [1]:
        return True
    else:
        return False 

def LOSNLOSxml(js, idx, TxPoint, RxPoint):
    """This function reads the newly generated geojson file(origin: lust.poly.xml, coordinates:x-y), and returns the number of walls and crossing
     distance of Tx and Rx transmission. If n == 0, the transmission is in LOS.  NOTE: Fast version, the rtree inserting in maincode."""
    """input:
        js:                 python object, after reading the geojson file, unit: meter  
        idx:                rtree index after insert the building bounding box   
        TxPoint             a tuple of x/y, instance:  TxPoint = (6826.12, 7137.92)
        RxPoint             a tuple of x/y, instance:  RxPoint = (6923.42, 7153.99)
       output:
        d:                  total length of the obstacle's intersection, in meters
        n:                  the number of borders of the building which intersect the line of sight 
    """

    # the distance and wall numbers of intersection
    distance = []
    wall_num = []
    # get the LineString of TxPoint and RxPoint
    Line1 = LineString([TxPoint, RxPoint])
    # get the bounding box of TxPoint and RxPoint 
    (minx, miny, maxx, maxy) = shape(Line1).bounds
    Box1 = box(minx, miny, maxx, maxy)
    # rough intersection of buildings'bounding box and Tx/RxPoints bounding box
    for pos in idx.intersection(Box1.bounds):
    # for pos in idx.intersection(Box1.bounds):
        feature = js['features'][pos]
        polygon = shape(feature['geometry'])
        # get the exact intersection line or geometry collection of Tx/RxPoint and buildings
        intersect_line = polygon.intersection(Line1)
        # If there is an intersection in any builiding, thus a LineString/MultiLineString exists
        if intersect_line:
            # If it is a MultiLineString Object, with crossing multiparts in the building
            if isinstance(intersect_line, shapely.geometry.multilinestring.MultiLineString):
                # calculate the crossing distance in the building & calculate the number of crossing walls
                walls = 2*len(intersect_line.geoms)
                distance.append(intersect_line.length)
                wall_num.append(walls)
            else:
                # else, just a LineString
                inter_point = list(intersect_line.coords)
                # calculate the crossing distance in the building, in meters
                inter_distance = caldistance(inter_point[0], inter_point[1])
                # calculate the crossing walls in the building
                walls = len(inter_point)
                distance.append(inter_distance)
                wall_num.append(walls)
    d = sum(distance)
    n = sum(wall_num)   
    return d, n

def LOSNLOSxy(inputfile, TxPoint, RxPoint):
    """This function reads the building polygons file(geojson file, from xml, unit: meter), and returns the number of walls and crossing distance of Tx 
    and Rx transmission. If n == 0, the transmission is in LOS. NOTE: Slow version, each time read all builidings."""
    """input:
        inputfile           region of interest, geojson file, unit: meter
        TxPoint             a tuple of x/y, instance:  TxPoint = (6826.12, 7137.92)
        RxPoint             a tuple of x/y, instance:  RxPoint = (6923.42, 7153.99)
       output:
        d:                  total length of the obstacle's intersection, in meter
        n:                  the number of borders of the building which intersect the line of sight 
    """
    # the transmission line of transmitting node and receiving node
    Line1 = LineString([TxPoint, RxPoint])
    with open(inputfile, 'r') as f:
        js = json.load(f)  
    
    # # debug output: the number of obstacle buildings
    # i = 0
    # the distance and wall numbers of intersection
    distance = []
    wall_num = []

    for feature in js['features']:
        polygon = shape(feature['geometry'])
        intersect_line = polygon.intersection(Line1)
        if intersect_line:
            # i += 1
            # print("intersection building:", i)
            # print(intersect_line)
            if isinstance(intersect_line, shapely.geometry.multilinestring.MultiLineString):
                intersect_line = polygon.intersection(Line1)
                # calculate the crossing distance in the building and the number of crossing walls
                walls = 2*len(intersect_line.geoms)
                distance.append(intersect_line.length)
                wall_num.append(walls)
            else:
                inter_point = list(intersect_line.coords)
                # calculate the crossing distance in the building, in meters
                inter_distance = caldistance(inter_point[0], inter_point[1])
                # print("intersect_distance:", inter_distance)
                # calculate the crossing walls in the building
                walls = len(inter_point)
                # print("crossing", walls, "walls")
                distance.append(inter_distance)
                wall_num.append(walls)
    d = sum(distance)
    n = sum(wall_num)   
    return d, n

def LOSNLOS(inputfile, TxPoint, RxPoint):
    """This function reads the building polygons file(geojson file), and returns the number of walls and crossing distance of Tx and Rx transmission.
    If n == 0, the transmission is in LOS. NOTE: origin version, lat/lon, read all buildings everytime."""
    """input:
        inputfile           region of interest, geojson file, unit: degree
        TxPoint             a tuple of lon/lat     TxPoint = (6.1034, 49.61828)
        RxPoint             a tuple of lon/lat     RxPoint = (6.1042, 49.61826)
       output:
        d:                  total length of the obstacle's intersection
        n:                  the number of borders of the building which intersect the line of sight 
    """
    # the transmission line of transmitting node and receiving node
    Line1 = LineString([TxPoint, RxPoint])
    with open(inputfile, 'r') as f:
        js = json.load(f)  
    
    # # the number of obstacle buildings
    # i = 0
    # the distance and wall numbers of intersection
    distance = []
    wall_num = []

    for feature in js['features']:
        polygon = shape(feature['geometry'])
        intersect_line = polygon.intersection(Line1)
        if intersect_line:
            # i += 1
            # print("intersection building:", i)
            # print(intersect_line)
            if isinstance(intersect_line, shapely.geometry.multilinestring.MultiLineString):
                intersect_line = polygon.intersection(Line1)
                # convert the distance of degrees to meters
                project = partial(pyproj.transform, pyproj.Proj(init="EPSG:4326"), pyproj.Proj(init="EPSG:32633"))
                line2 = transform(project, intersect_line)
                # calculate the crossing distance in the building
                # print("intersect_distance", line2.length, "meters")
                # calculate the number of crossing walls
                walls = 2*len(intersect_line.geoms)
                # print("crossing", walls, "walls")
                distance.append(line2.length)
                wall_num.append(walls)
            else:
                inter_point = list(intersect_line.coords)
                # calculate the crossing distance in the building, in meters
                inter_distance = vincenty(inter_point[0], inter_point[1]).meters
                # print("intersect_distance:", inter_distance)
                # calculate the crossing walls in the building
                walls = len(inter_point)
                # print("crossing", walls, "walls")
                distance.append(inter_distance)
                wall_num.append(walls)
    d = sum(distance)
    n = sum(wall_num)   
    return d, n

def NonLineofSight(Pt,txHeight, rxHeight, Gt, Gr, dist, d, n):
    """This function calculates the receiving power of NLOS transmission, based on two ray interference model and obstacle shadowing models."""
    """Input: 
        Pt:                     the transmit power, in mW  
        txHeight, rxHeight:     the antenna height of tx and rx nodes
        Gt, Gr:                 the antenna gain of tx and rx nodes, in dBi
        dist:                   distance between Tx and Rx
        d:                      total length of the obstacle's intersection
        n:                      the number of borders of the building which intersect the line of sight 
       Output:
        Pr:                     received power, in dBm
    """
    # Relative permitivity of the ground, depends on the scenario.
    er = 1.02
    # speed of light, in meters
    c = 299792458
    # CCH central frequency, in GHz
    f = 5.890
    # wave length, in meters
    lamda = c/f*10**(-9)
    # LOS distance
    dLos = sqrt((txHeight-rxHeight)**2 + dist**2)
    # Ground reflected distance
    dRef = sqrt((txHeight+rxHeight)**2 + dist**2)
    # Sine and cosine of incident angle Theta
    sinTheta = (txHeight+rxHeight)/dRef
    cosTheta = dLos/dRef
    # Phase difference
    phDif = 2 * pi * (dLos-dRef) / lamda 
    # Power in W
    # Pt = (10**(Pt/10))/1000
    # Convert antenna gains from dB
    Gt = 10**(Gt/10)
    Gr = 10**(Gr/10)
    # Reflection coefficient, horizontal polarization
    GM = (sinTheta- sqrt(er-cosTheta**2)) / (sinTheta+ sqrt(er-cosTheta**2))
    # Two ray interference path loss component, unit in dB
    att_2ray = pow(4* pi*(dist/lamda)*1/(sqrt((pow((1 + GM* cos(phDif)),2) + pow(GM,2)*pow(sin(phDif),2)))),2)
    # no antenna influences
    AvePr = Pt / att_2ray
    ########################################################
    # comment out Gt.Gr influences
    # # Receving power in W after two ray interference model
    # AvePr = Pt * Gt * Gr / att_2ray 
    ########################################################
    # y: serves as a rough approximation of the internal structure of a building, in dB per meter
    y = 0.4
    # beta: represents the attenuation a transmission experecnes due to the exterior wall of a building, in dB per meter
    beta = 9
    # obstacle: attenuation of obstacle shadowing 
    obstacle = beta*n + y*d
    # received power in dBm, after obstacle model
    Pr = 10*(log10(AvePr/obstacle))
    return Pr

def LineofSight(Pt, txHeight, rxHeight, Gt, Gr, dist):
    """This function calculates the receving power of LOS transmission, based on two ray interference model and nakagami-m fading model"""
    """Input:
        Pt:                     transmitting power in mW
        txHeight, rxHeight:     the antenna height of tx and rx nodes
        Gt, Gr:                 the antenna gain of tx and rx nodes, in dBi
        dist:                   distance between Tx and Rx

    Output:
        Pr:                     received power in dBm 
    """
    # Relative permitivity of the ground, depends on the scenario.
    er = 1.02
    # speed of light, in meters
    c = 299792458
    # CCH central frequency, in GHz
    f = 5.890
    # wave length, in meters
    lamda = c/f*10**(-9)
    # LOS distance
    dLos = sqrt((txHeight-rxHeight)**2 + dist**2)
    # Ground reflected distance
    dRef = sqrt((txHeight+rxHeight)**2 + dist**2)
    # Sine and cosine of incident angle Theta
    sinTheta = (txHeight+rxHeight)/dRef
    cosTheta = dLos/dRef
    # Phase difference
    phDif = 2 * pi * (dLos-dRef) / lamda 
    # Power in W
    # Pt = (10**(Pt/10))/1000
    # Convert antenna gains from dB
    Gt = 10**(Gt/10)
    Gr = 10**(Gr/10)
    # Reflection coefficient, horizontal polarization
    GM = (sinTheta- sqrt(er-cosTheta**2)) / (sinTheta+ sqrt(er-cosTheta**2))
    # Two ray interference path loss component, unit in dB
    att_2ray = pow(4* pi*(dist/lamda)*1/(sqrt((pow((1 + GM* cos(phDif)),2) + pow(GM,2)*pow(sin(phDif),2)))),2)
    # no antenna gain 
    AvePr = Pt / att_2ray
    # ####################################################
    # comment out antenna gain influences
    # Receving power in W after two ray interference model
    # AvePr = Pt * Gt * Gr / att_2ray 
    ######################################################
    #############################################################################
    # this section is about Nakagami-fading model  
    # # Nakagami-m parameter, the distance threshold from veins code, in meters
    # Dist_thre = 80
    # # nakagami-m fading, m-close & m-far
    # M_close = 1.5
    # M_far = 0.75
    # if (dist <= Dist_thre):
    #     m = M_close
    # else:
    #     m = M_far
    # # Nakagami-m distribution, shape-m, location-AvePr/m 
    # NakAtt = stat.gamma.rvs(m, AvePr/m)
    # # received power in dBm, after Nakagami-m fading
    # Prec = 10*(log10(AvePr/NakAtt)) + 30
    ############################################################################
    Prec = 10 * log10(AvePr)
    # Prec = 10*log10(AvePr) + 30
    return Prec
