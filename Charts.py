from astropy.io.votable import parse,writeto,parse_single_table,from_table
from astropy.coordinates import SkyCoord
from astropy import units as u 
import util as ut
import importlib

importlib.reload(ut)


Katalog1=parse_single_table("Test_Source3.xml")
Katalog=Katalog1.to_table()

KatalogDict={}

print(281)
for x in range(281):#len(Katalog["Name"])):
	output_dict={}
	output_dict["Coordinates"]=SkyCoord(ra=Katalog["RA"][x],dec=Katalog["DEC"][x],unit="deg")
	output_dict["Coord_Error"]=Katalog["Coord_Error"][x]
	output_dict["Coord_Error_Source"]=Katalog["Coord_Error_Source"][x]
	output_dict["XMM_RA"]=Katalog["XMM_RA"][x]
	output_dict["XMM_DEC"]=Katalog["XMM_DEC"][x]
	output_dict["XMM_Error"]=Katalog["XMM_Error"][x]
	output_dict["Integral_RA"]=Katalog["Integral_RA"][x]
	output_dict["Integral_DEC"]=Katalog["Integral_DEC"][x]
	output_dict["BAT_RA"]=Katalog["BAT_RA"][x]
	output_dict["BAT_DEC"]=Katalog["BAT_DEC"][x]
	output_dict["Chandra_RA"]=Katalog["Chandra_RA"][x]
	output_dict["Chandra_DEC"]=Katalog["Chandra_DEC"][x]
	output_dict["Chandra_Error"]=Katalog["Chandra_Error"][x]
	output_dict["XRT_RA"]=Katalog["XRT_RA"][x]
	output_dict["XRT_DEC"]=Katalog["XRT_DEC"][x]
	output_dict["XRT_Error"]=Katalog["XRT_Error"][x]
	output_dict["SimbadRA"]=Katalog["SimbadRA"][x]
	output_dict["SimbadDEC"]=Katalog["SimbadDEC"][x]
	output_dict["twoMassRA"]=Katalog["_2MASS_RA"][x]
	output_dict["twoMassDEC"]=Katalog["_2MASS_DEC"][x]
	output_dict["twoMassError"]=Katalog["_2MASS_Error"][x]
	output_dict["GAIA_RA"]=Katalog["GAIA_RA"][x]
	output_dict["GAIA_DEC"]=Katalog["GAIA_DEC"][x]
	output_dict["GAIA_Error"]=Katalog["GAIA_Error"][x]

	output_dict["RA"]=Katalog["RA"][x]
	output_dict["DEC"]=Katalog["DEC"][x]
	KatalogDict[Katalog["Name"][x]]=output_dict


uu=1
for key in KatalogDict:
	print(uu)
	try:
		ut.mk_finding_chart(KatalogDict[key]["SimbadRA"].data,KatalogDict[key]["SimbadDEC"].data,KatalogDict[key]["Coordinates"],KatalogDict[key]["Coord_Error"].data,KatalogDict[key]["Coord_Error_Source"],key,KatalogDict[key]["XMM_RA"].data,
			KatalogDict[key]["XMM_DEC"].data,KatalogDict[key]["XMM_Error"].data,KatalogDict[key]["XRT_Error"].data,KatalogDict[key]["XRT_RA"].data,KatalogDict[key]["XRT_DEC"].data,KatalogDict[key]["BAT_RA"].data,KatalogDict[key]["BAT_DEC"].data,KatalogDict[key]["Integral_RA"].data,KatalogDict[key]["Integral_DEC"].data,
			KatalogDict[key]["Chandra_RA"].data,KatalogDict[key]["Chandra_DEC"].data,KatalogDict[key]["Chandra_Error"].data,KatalogDict[key]["twoMassRA"].data,KatalogDict[key]["twoMassDEC"].data,KatalogDict[key]["twoMassError"].data,
			KatalogDict[key]["GAIA_RA"].data,KatalogDict[key]["GAIA_DEC"].data,KatalogDict[key]["GAIA_Error"].data)
	except (IndexError,ValueError):
		pass
	uu+=1

