import arcpy, time, math, os
import sys
reload(sys)
sys.setdefaultencoding("latin-1")
start_time = time.time()

# http://stackoverflow.com/questions/343865/how-to-convert-from-utm-to-latlng-in-python-or-javascript
'''
datumTable: [
        { eqRad: 6378137.0, flat: 298.2572236 },    // WGS 84
        { eqRad: 6378137.0, flat: 298.2572236 },    // NAD 83
        { eqRad: 6378137.0, flat: 298.2572215 },    // GRS 80
        { eqRad: 6378135.0, flat: 298.2597208 },    // WGS 72
        { eqRad: 6378160.0, flat: 298.2497323 },    // Austrailian 1965
        { eqRad: 6378245.0, flat: 298.2997381 },    // Krasovsky 1940
        { eqRad: 6378206.4, flat: 294.9786982 },    // North American 1927
        { eqRad: 6378388.0, flat: 296.9993621 },    // International 1924
        { eqRad: 6378388.0, flat: 296.9993621 },    // Hayford 1909
        { eqRad: 6378249.1, flat: 293.4660167 },    // Clarke 1880
        { eqRad: 6378206.4, flat: 294.9786982 },    // Clarke 1866
        { eqRad: 6377563.4, flat: 299.3247788 },    // Airy 1830
        { eqRad: 6377397.2, flat: 299.1527052 },    // Bessel 1841
        { eqRad: 6377276.3, flat: 300.8021499 }     // Everest 1830
    ],
'''
def utmToLatLng(zone, easting, northing, datum, northernHemisphere=True):
	if not northernHemisphere:
		northing = 10000000 - northing
	datumDict = {'NAD83': {'eqRad': 6378137.0, 'flat': 298.2572236}, 'NAD27': {'eqRad': 6378206.4, 'flat': 294.9786982}};
	if (datum in datumDict):
		a = datumDict[datum]['eqRad']
		flat = datumDict[datum]['flat']
	else:
		a = datumDict['NAD83']['eqRad']
		flat = datumDict['NAD83']['flat']
		
	# f = 1/flat; e1sq = 2 * f - f * f; e = math.sqrt(e1sq)
	f = 1/flat
	e1sq = 2 * f - f * f
	e = math.sqrt(e1sq)
	#a = 6378137
	#e = 0.081819191
	#e1sq = 0.006739497
	k0 = 0.9996
		
	arc = northing / k0
	mu = arc / (a * (1 - math.pow(e, 2) / 4.0 - 3 * math.pow(e, 4) / 64.0 - 5 * math.pow(e, 6) / 256.0))

	ei = (1 - math.pow((1 - e * e), (1 / 2.0))) / (1 + math.pow((1 - e * e), (1 / 2.0)))

	ca = 3 * ei / 2 - 27 * math.pow(ei, 3) / 32.0

	cb = 21 * math.pow(ei, 2) / 16 - 55 * math.pow(ei, 4) / 32
	cc = 151 * math.pow(ei, 3) / 96
	cd = 1097 * math.pow(ei, 4) / 512
	phi1 = mu + ca * math.sin(2 * mu) + cb * math.sin(4 * mu) + cc * math.sin(6 * mu) + cd * math.sin(8 * mu)

	n0 = a / math.pow((1 - math.pow((e * math.sin(phi1)), 2)), (1 / 2.0))

	r0 = a * (1 - e * e) / math.pow((1 - math.pow((e * math.sin(phi1)), 2)), (3 / 2.0))
	fact1 = n0 * math.tan(phi1) / r0

	_a1 = 500000 - easting
	dd0 = _a1 / (n0 * k0)
	fact2 = dd0 * dd0 / 2

	t0 = math.pow(math.tan(phi1), 2)
	Q0 = e1sq * math.pow(math.cos(phi1), 2)
	fact3 = (5 + 3 * t0 + 10 * Q0 - 4 * Q0 * Q0 - 9 * e1sq) * math.pow(dd0, 4) / 24

	fact4 = (61 + 90 * t0 + 298 * Q0 + 45 * t0 * t0 - 252 * e1sq - 3 * Q0 * Q0) * math.pow(dd0, 6) / 720

	lof1 = _a1 / (n0 * k0)
	lof2 = (1 + 2 * t0 + Q0) * math.pow(dd0, 3) / 6.0
	lof3 = (5 - 2 * Q0 + 28 * t0 - 3 * math.pow(Q0, 2) + 8 * e1sq + 24 * math.pow(t0, 2)) * math.pow(dd0, 5) / 120
	_a2 = (lof1 - lof2 + lof3) / math.cos(phi1)
	_a3 = _a2 * 180 / math.pi

	latitude = 180 * (phi1 - fact1 * (fact2 + fact3 + fact4)) / math.pi

	if not northernHemisphere:
		latitude = -latitude

	longitude = ((zone > 0) and (6 * zone - 183.0) or 3.0) - _a3

	return (latitude, longitude)

def createPolygonFeatureClass(featureName, featureData, featureFieldList, featureInsertCursorFields):
	print "Create " + featureName + " feature class"
	featureNameNAD83 = featureName + "_NAD83"
	featureNameNAD83Path = arcpy.env.workspace + "\\"  + featureNameNAD83
	arcpy.CreateFeatureclass_management(arcpy.env.workspace, featureNameNAD83, "POLYGON", "", "DISABLED", "DISABLED", "", "", "0", "0", "0")
	# Process: Define Projection
	arcpy.DefineProjection_management(featureNameNAD83Path, "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]")
	# Process: Add Fields	
	for featrueField in featureFieldList:
		arcpy.AddField_management(featureNameNAD83Path, featrueField[0], featrueField[1], featrueField[2], featrueField[3], featrueField[4], featrueField[5], featrueField[6], featrueField[7], featrueField[8])
	# Process: Append the records
	cntr = 1
	try:
		with arcpy.da.InsertCursor(featureNameNAD83, featureInsertCursorFields) as cur:
			for rowValue in featureData:
				cur.insertRow(rowValue)
				cntr = cntr + 1
	except Exception as e:
		print "\tError: " + featureName + ": " + e.message
	# Change the projection to web mercator
	arcpy.Project_management(featureNameNAD83Path, arcpy.env.workspace + "\\" + featureName, "PROJCS['WGS_1984_Web_Mercator_Auxiliary_Sphere',GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Mercator_Auxiliary_Sphere'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',0.0],PARAMETER['Standard_Parallel_1',0.0],PARAMETER['Auxiliary_Sphere_Type',0.0],UNIT['Meter',1.0]]", "NAD_1983_To_WGS_1984_5", "GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]")
	arcpy.Delete_management(featureNameNAD83Path, "FeatureClass")
	print "Finish " + featureName + " feature class."

INPUT_PATH = "input"
OUTPUT_PATH = "output"
if arcpy.Exists(OUTPUT_PATH + "\\DEMGrids.gdb"):
	os.system("rmdir " + OUTPUT_PATH + "\\DEMGrids.gdb /s /q")
os.system("del " + OUTPUT_PATH + "\\*DEMGrids*.*")
arcpy.CreateFileGDB_management(OUTPUT_PATH, "DEMGrids", "9.3")
arcpy.env.workspace = OUTPUT_PATH + "\\DEMGrids.gdb"

inputFile = open('input/index.html')
input = inputFile.read()
htmlSnippets = input.split('<TR>')
del htmlSnippets[0]
del htmlSnippets[0]

def createPolygon(htmlSnippet):
	items = htmlSnippet.split('<TD VALIGN="TOP" class="style1"  >')
	del items[0]
	tile = items[0].strip().split(">")[2].split("<")[0]
	city = items[1].strip()
	city = '' if city == '&nbsp;</TD>' else city.split(">")[1].split("<")[0]
	zone = int(items[2].strip().split(">")[1].split("<")[0])
	minUTMEasting = int(items[3].strip().split(">")[1].split("<")[0])
	maxUTMEasting = int(items[4].strip().split(">")[1].split("<")[0])
	minUTMNorthing = int(items[5].strip().split(">")[1].split("<")[0])
	maxUTMNorthing = int(items[6].strip().split(">")[1].split("<")[0])
	latlngWestSouth = utmToLatLng(zone, minUTMEasting, minUTMNorthing, 'NAD83')
	latlngWestNorth = utmToLatLng(zone, minUTMEasting, maxUTMNorthing, 'NAD83')
	latlngEastNorth = utmToLatLng(zone, maxUTMEasting, maxUTMNorthing, 'NAD83')
	latlngEastSouth = utmToLatLng(zone, maxUTMEasting, minUTMNorthing, 'NAD83')
	pnt = arcpy.Point()
	ary = arcpy.Array()
	pnt.X = latlngWestSouth[1]
	pnt.Y = latlngWestSouth[0]
	ary.add(pnt)
	pnt.X = latlngWestNorth[1]
	pnt.Y = latlngWestNorth[0]
	ary.add(pnt)
	pnt.X = latlngEastNorth[1]
	pnt.Y = latlngEastNorth[0]
	ary.add(pnt)
	pnt.X = latlngEastSouth[1]
	pnt.Y = latlngEastSouth[0]
	ary.add(pnt)
	polygon = arcpy.Polygon(ary)
	return [polygon, tile, city, str(zone)]
featureData = map(createPolygon, htmlSnippets)
featureName = "DEMGrids"
featureFieldList = [["TILE", "TEXT", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", ""], ["CITY", "TEXT", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", ""], ["UTMZONE", "TEXT", "", "", "", "", "NON_NULLABLE", "NON_REQUIRED", ""]]
featureInsertCursorFields = ("SHAPE@", "TILE", "CITY", "UTMZONE")
createPolygonFeatureClass(featureName, featureData, featureFieldList, featureInsertCursorFields)

elapsed_time = time.time() - start_time
print elapsed_time