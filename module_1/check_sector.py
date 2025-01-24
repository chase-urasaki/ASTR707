from tess_stars2px import tess_stars2px_function_entry
import argparse
from astroquery.mast import Catalogs

p = argparse.ArgumentParser()
p.add_argument("objectname",nargs=1)

args = p.parse_args()
objectname = args.objectname
objectname = str(objectname[0])
radSearch = 5./3600.
catalogData = Catalogs.query_object(objectname, radius = radSearch, catalog = "TIC")
ticid = int(catalogData[0]['ID'])
print('TIC counterpart = ',ticid)
ra = catalogData[0]['ra']
dec = catalogData[0]['dec']

outID, outEclipLong, outEclipLat,outSec,outCam,outCcd,outColPix,outRowPix,scinfo = tess_stars2px_function_entry(ticid,ra,dec)
if (len(outID) > 0):
    print('Sector','Camera','CCD','Column','Row')
    for u,v,w,p,q in zip(outSec,outCam,outCcd,outColPix,outRowPix):
        print(u,v,w,p,q)
else:
    print('TESS has not observed this star')


