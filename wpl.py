import math

def main():


  arrayCentre = [147.215965, -41.547531]  # Lon, Lat
  blp, CA = blp_array_location(arrayCentre, rotation = 180.0)

  arrayCentre = [147.215965, -41.547531]  # Lon, Lat
  arc, az = map_2points(147.215414, -41.547638, 147.215606, -41.547849, degrees=True)
  az -= 90

  lon, lat = ll_arc_distance(arrayCentre[0], arrayCentre[1], 10.0 / 6371000.0, az, degrees=True)
  arrayCentre = [lon, lat]
  blp, CA = blp_array_location(arrayCentre, rotation = 180.0)

  kml_output(blp, CA)





def blp_array_location(arrayCentre, rotation = 0.0, antennaSpacing = 2.725, antennaLength = 2.814, padding = 1.0):
  RE = 6371000.0 # Radius of the earth in meters
  sqrt2 = math.sqrt(2)

  blp = [ [], [], [] ]
  compactionArea = []

  # Define the bottom left corners:
  # First array
  lon = arrayCentre[0]
  lat = arrayCentre[1]
  x = (antennaSpacing + antennaLength / (2 * sqrt2))
  y = (antennaSpacing / 2 - antennaLength / (2 * sqrt2))
  lon, lat = ll_arc_distance(lon, lat, x/RE, 270.0 + rotation, degrees = True)
  lon, lat = ll_arc_distance(lon, lat, y/RE, 0.0 + rotation, degrees = True)
  blp[0].append([lon, lat])

  # Second array
  lon = arrayCentre[0]
  lat = arrayCentre[1]
  x = (antennaSpacing * 5 / 2 + antennaLength / (2 * sqrt2))
  y = (antennaSpacing * 5 / 2 + antennaLength / (2 * sqrt2))
  lon, lat = ll_arc_distance(lon, lat, x/RE, 270.0 + rotation, degrees = True)
  lon, lat = ll_arc_distance(lon, lat, y/RE, 180.0 + rotation, degrees = True)
  blp[1].append([lon, lat])

  # Third array:
  lon = arrayCentre[0]
  lat = arrayCentre[1]
  x = (antennaSpacing / 2 - antennaLength / (2 * sqrt2))
  y = (antennaSpacing * 5 / 2 + antennaLength / (2 * sqrt2))
  lon, lat = ll_arc_distance(lon, lat, x/RE, 90.0 + rotation, degrees = True)
  lon, lat = ll_arc_distance(lon, lat, y/RE, 180.0 + rotation, degrees = True)
  blp[2].append([lon, lat])


  # Trace the three arrays from the bottom left corner
  for i in range(3):
    # Top left corner:
    lon = blp[i][-1][0]
    lat = blp[i][-1][1]
    x = 0
    y = (2 * antennaSpacing + antennaLength / sqrt2)
    lon, lat = ll_arc_distance(lon, lat, x/RE, 90.0 + rotation, degrees = True)
    lon, lat = ll_arc_distance(lon, lat, y/RE, 0.0 + rotation, degrees = True)
    blp[i].append([lon, lat])

    # Top right corner:
    lon = blp[i][-1][0]
    lat = blp[i][-1][1]
    x = (2 * antennaSpacing + antennaLength / sqrt2)
    y = 0
    lon, lat = ll_arc_distance(lon, lat, x/RE, 90.0 + rotation, degrees = True)
    lon, lat = ll_arc_distance(lon, lat, y/RE, 0.0 + rotation, degrees = True)
    blp[i].append([lon, lat])

    # Bottom right corner:
    lon = blp[i][-1][0]
    lat = blp[i][-1][1]
    x = 0
    y = (2 * antennaSpacing + antennaLength / sqrt2)
    lon, lat = ll_arc_distance(lon, lat, x/RE, 90.0 + rotation, degrees = True)
    lon, lat = ll_arc_distance(lon, lat, y/RE, 180.0 + rotation, degrees = True)
    blp[i].append([lon, lat])





  # Start at bottom left, move clockwise
  #Identify point (1)
  lon = arrayCentre[0]
  lat = arrayCentre[1]
  x = antennaSpacing * 5 / 2 + antennaLength / (2 * sqrt2) + padding
  y = antennaSpacing * 5 / 2 + antennaLength / (2 * sqrt2) + padding
  lon, lat = ll_arc_distance(lon, lat, x/RE, 270.0 + rotation, degrees = True)
  lon, lat = ll_arc_distance(lon, lat, y/RE, 180.0 + rotation, degrees = True)
  compactionArea.append([lon,lat])

  # (1) -> (2)
  lon = compactionArea[-1][0]
  lat = compactionArea[-1][1]
  x = 0.0
  y = antennaSpacing * 2 + antennaLength / sqrt2 + padding * 2.0
  lon, lat = ll_arc_distance(lon, lat, x/RE, 90.0 + rotation, degrees = True)
  lon, lat = ll_arc_distance(lon, lat, y/RE, 0.0 + rotation, degrees = True)
  compactionArea.append([lon,lat])

  # (2) -> (3)
  lon = compactionArea[-1][0]
  lat = compactionArea[-1][1]
  x = antennaSpacing * 3 / 2.0
  y = 0.0
  lon, lat = ll_arc_distance(lon, lat, x/RE, 90.0 + rotation, degrees = True)
  lon, lat = ll_arc_distance(lon, lat, y/RE, 0.0 + rotation, degrees = True)
  compactionArea.append([lon,lat])

  # (3) -> (4)
  lon = compactionArea[-1][0]
  lat = compactionArea[-1][1]
  x = 0.0
  y = antennaSpacing * 3
  lon, lat = ll_arc_distance(lon, lat, x/RE, 90.0 + rotation, degrees = True)
  lon, lat = ll_arc_distance(lon, lat, y/RE, 0.0 + rotation, degrees = True)
  compactionArea.append([lon,lat])

  # (4) -> (5)
  lon = compactionArea[-1][0]
  lat = compactionArea[-1][1]
  x = antennaSpacing * 2.0 + antennaLength / sqrt2 + padding * 2
  y = 0.0
  lon, lat = ll_arc_distance(lon, lat, x/RE, 90.0 + rotation, degrees = True)
  lon, lat = ll_arc_distance(lon, lat, y/RE, 0.0 + rotation, degrees = True)
  compactionArea.append([lon,lat])

  # (5) -> (6)
  lon = compactionArea[-1][0]
  lat = compactionArea[-1][1]
  x = 0.0
  y = antennaSpacing * 3
  lon, lat = ll_arc_distance(lon, lat, x/RE, 90.0 + rotation, degrees = True)
  lon, lat = ll_arc_distance(lon, lat, y/RE, 180.0 + rotation, degrees = True)
  compactionArea.append([lon,lat])

  # (6) -> (7)
  lon = compactionArea[-1][0]
  lat = compactionArea[-1][1]
  x = antennaSpacing * 3 / 2.0
  y = 0.0
  lon, lat = ll_arc_distance(lon, lat, x/RE, 90.0 + rotation, degrees = True)
  lon, lat = ll_arc_distance(lon, lat, y/RE, 0.0 + rotation, degrees = True)
  compactionArea.append([lon,lat])

  # (7) -> (8)
  lon = compactionArea[-1][0]
  lat = compactionArea[-1][1]
  x = 0.0
  y = antennaSpacing * 2 + antennaLength / sqrt2 + padding * 2.0
  lon, lat = ll_arc_distance(lon, lat, x/RE, 90.0 + rotation, degrees = True)
  lon, lat = ll_arc_distance(lon, lat, y/RE, 180.0 + rotation, degrees = True)
  compactionArea.append([lon,lat])






  return blp, compactionArea


def kml_output(blp, CA):
  d = 1.8 # m

  fileName = "BLP.kml"

  KML_STRING=""
  KML_STRING+="<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
  KML_STRING+=" <kml xmlns=\"http://www.opengis.net/kml/2.2\">\n"
  KML_STRING+="    <Document>\n"
  KML_STRING+="  <name>Paths</name>\n"
  KML_STRING+="  <description>\n"
  KML_STRING+="    "+fileName+"\n"
  KML_STRING+="  </description>\n"
  KML_STRING+="  <Style id=\"yellowLineGreenPoly\">\n"
  KML_STRING+="    <LineStyle>\n"
  KML_STRING+="      <color>7f00ffff</color>\n"
  KML_STRING+="      <width>4</width>\n"
  KML_STRING+="    </LineStyle>\n"
  KML_STRING+="    <PolyStyle>\n"
  KML_STRING+="      <color>7f00ff00</color>\n"
  KML_STRING+="    </PolyStyle>\n"
  KML_STRING+="  </Style>\n"
  KML_STRING+="  <Placemark>\n"
  KML_STRING+="  <MultiGeometry>\n"
  KML_STRING+="    <name>BLP</name>\n"
  KML_STRING+="    <description>BLP</description>\n"
  KML_STRING+="    <styleUrl>#yellowLineGreenPoly</styleUrl>\n"

  for arr in blp:
    KML_STRING+="    <Polygon>\n"
    KML_STRING+="      <extrude>0</extrude>\n"
    KML_STRING+="      <tessellate>1</tessellate>\n"
    KML_STRING+="      <altitudeMode>relativeToGround</altitudeMode>\n"
    KML_STRING+="      <outerBoundaryIs>"
    KML_STRING+="      <LinearRing>"
    KML_STRING+="      <coordinates>\n"

    for p in arr:
      KML_STRING+="       " + str(p[0]) + "," + str(p[1]) + "," + str(d) + "\n"
    KML_STRING+="       " + str(arr[0][0]) + "," + str(arr[0][1]) + "," + str(d) + "\n"

    KML_STRING+="        </coordinates>\n"
    KML_STRING+="        </LinearRing>\n"
    KML_STRING+="        </outerBoundaryIs>\n"
    KML_STRING+="      </Polygon>\n"



  KML_STRING+="    <Polygon>\n"
  KML_STRING+="      <extrude>0</extrude>\n"
  KML_STRING+="      <tessellate>1</tessellate>\n"
  KML_STRING+="      <altitudeMode>relativeToGround</altitudeMode>\n"
  KML_STRING+="      <outerBoundaryIs>"
  KML_STRING+="      <LinearRing>"
  KML_STRING+="      <coordinates>\n"
  for p in CA:
    KML_STRING+="       " + str(p[0]) + "," + str(p[1]) + "," + str(d-0.5) + "\n"
  KML_STRING+="       " + str(CA[0][0]) + "," + str(CA[0][1]) + "," + str(d-0.5) + "\n"
  KML_STRING+="        </coordinates>\n"
  KML_STRING+="        </LinearRing>\n"
  KML_STRING+="        </outerBoundaryIs>\n"
  KML_STRING+="      </Polygon>\n"

  KML_STRING+="    </MultiGeometry>\n"
  KML_STRING+="    </Placemark>\n"
  KML_STRING+="   </Document>\n"
  KML_STRING+=" </kml>\n"

  try:
    kmlFileName = fileName.replace(".txt",".kml")
    with open(kmlFileName,"w") as f:
      f.write(KML_STRING)
      f.write("\n")
  except:
    raise








###############################################################################
#
#  Function: ll_arc_distance
#
#  Calculates the longitude and latitude on a sphere some arc distance from a
#  given starting point at some azimuth from North. Arc distance is always in
#  radians, other parameters in radians unless degrees is specified.
#
#  Parameters:
#    lon - Initial longitude
#    lat - Initial latitude
#    arc - Arc distance to move (always radians)
#    az - azimuthal direction to move (0 is north)
#
#  Keywords:
#    degrees - If set, lon/lat/az are used and returned in degrees.
#
#  Returns:
#    lon,lat - The final longitude and latitude
#
###############################################################################
import math
import numpy as np
def ll_arc_distance(lon,lat,arc,az,degrees=False):
  if arc==0: return lon,lat
  k=math.pi/180.0 if degrees else 1.0
  cdist=math.cos(arc)  #Arc is always radians
  sdist=math.sin(arc)
  clat=math.cos(lat*k)
  slat=math.sin(lat*k)
  caz=math.cos(k*az)
  phi=np.arcsin(slat*cdist + clat*sdist*caz)
  lam=(k*lon)+np.arctan2(sdist*math.sin(k*az),clat*cdist-slat*sdist*caz)
  while lam<=(-1*np.pi): lam+=2*np.pi
  while lam>np.pi: lam-=2*np.pi
  return lam/k,phi/k



###############################################################################
#
# Function: map_2points
#
# Calculates the azimuth and arclength between two points on a sphere.
#
# Parameters:
#   lon0 - Initial longitude
#   lat0 - Initial latitude
#   lon1 - Final longitude
#   lat1 - Final latitude
#
# Keywords:
#   degrees - specify that inputs are in degrees and outputs should be in
#             degrees. Note that arclength is always in radians.
#             (default: False)
#
# Returns:
#   arc - Arc length between the two points
#   az - Azimuthal direction from initial to final point
#
#
###############################################################################
def map_2points(lon0,lat0,lon1,lat1,degrees=False):
  """Determine the arc and azimuth between two points on a sphere"""
  k=1.0
  if degrees: k=math.radians(1.0) #convert to radians
  if abs(k*lat0)>(math.pi/2.) or abs(k*lat1)>(math.pi/2.):
    raise ValueError('Latitude out of range')
  clat0=math.cos(k*lat0)
  slat0=math.sin(k*lat0)
  clat1=math.cos(k*lat1)
  slat1=math.sin(k*lat1)
  cdlon=math.cos(k*lon1-k*lon0)
  sdlon=math.sin(k*lon1-k*lon0)
  cosc=(slat0*slat1)+(clat0*clat1*cdlon)
  if cosc<-1.0: cosc=-1.0
  if cosc>1.0: cosc=1.0                 #Incase of rounding errors
  sinc=math.sqrt(1.0-(cosc*cosc))
  caz=((clat0*slat1)-(slat0*clat1*cdlon))/sinc
  saz=(sdlon*clat1)/sinc                #IDL checks for tiny angle to decide if point is antipodal
  arc=math.acos(cosc)
  az=math.atan2(saz,caz)
  return arc,az/k
###############################################################################
###############################################################################




main()
