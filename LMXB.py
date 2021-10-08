import matplotlib
import astropy.units as u 
import numpy as np 
from astropy.table import QTable 
import astroquery
import requests
from astroquery.vizier import Vizier
from astroquery.heasarc import Heasarc
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
from astropy.coordinates import Angle
from astropy.io.votable.tree import VOTableFile, Resource, Table, Field
import numpy.ma as ma
import matplotlib.pyplot as plt
from numpy import *
import pickle
import collections
from collections import OrderedDict
import json 
from astropy.io.votable import parse,writeto,parse_single_table,from_table
import math as m 
import os
from astroquery.gaia import Gaia
from itertools import groupby
import pandas as pd
heasarc = Heasarc()

import warnings
warnings.filterwarnings("ignore")

viz = Vizier(columns=['all', '+_r'])
viz.ROW_LIMIT = inf

Simb_v2 = Simbad()
Simb_v2.add_votable_fields('biblio',  'otypes', 'flux(V)', 'ids', 'coo_bibcode')
Simb_v2.TIMEOUT = 36000 # sets the timeout to an 10 hours


Rit_cat = 'B/cb/lmxbdata'
Rit_cat_note = 'B/cb/whoswho1'

Liu_cat = 'J/A+A/469/807/lmxb'
Liu_cat_note = 'J/A+A/469/807/lmxbnote'

with open('pulsar_list_table.txt', 'r') as text: 
	table_list = text.readlines() 

CRSF_dict={}

for x in range(9, 241): 
	output_dict={} 
	output_dict["Spin"] = table_list[x][126:137].strip() 
	output_dict["CRSF"] = table_list[x][137:143].strip() 
	CRSF_dict[table_list[x][14:38].strip()] = output_dict 

#------------------------------------------------------------------------------------------------------------------------------------------------------------

Simbadobjects = Simb_v2.query_criteria(otype = 'LXB')
#Simbadobjects = Simb_v2.query_criteria(otypes = 'LXB')
mylist = Simbadobjects['MAIN_ID'].data.data

subs = ['[OWB2005] NGC 4214 P3', 'Bol 158', '[FAB2017] NGC4697 S27', 'CXO J073753.0+653407', 'CXO J073654.3+654016', '[LI', '[Z', '2MAS', '[CHP', '[S',  '[B', '2XMM J0', '2XMM J12', '[LB', 'RX J13', 'RX J0', '[H2', 'EV', 'CXOU', 'CXOM3', '[I', '[KK', '[KG', '2E 18']
execp = ['2MASS J17', '[BS', '[SK']
compl = ['[ZGV2011] 1', '[ZGV2011] 3', '[ZGV2011] 4', '[ZGV2011] 5', '[ZGV2011] 6', '[ZGV2011] 9']

lis = mylist

for j in subs:
	lis = [i for i in lis if j not in i]

for h in execp:
	tr = [i for i in mylist if h in i]
	lis = lis + tr

for k in compl:
	kr = [i for i in mylist if k == i]
	lis = lis + kr

lis_arr = list(mylist)
ind = [lis_arr.index(i) for i in lis]
Simb_tab = Simbadobjects[ind]
names = Simb_tab['MAIN_ID'].data.data

Simbad_name_list = []

for nam in names:
	result_table = Simb_v2.query_object(nam)
	Simbad_name_list.append(result_table['IDS'][0])

#------------------------------------------------------------------------------------------------------------------------------------------------------------

Rit = viz.get_catalogs(catalog = Rit_cat)[0]
Liu = viz.get_catalogs(catalog = Liu_cat)[0]
Liu_note = viz.get_catalogs(catalog = Liu_cat_note)[0]


Name_1 = list(Rit['Name'].data.data)
Name_2_ = list(Liu['Name'].data.data)

# (AX J1745.6-2901, 1A 1742-289); (1E 1743.1-2852, SAX J1747.0-2853); (CXO J212958.1+121002, 4U 2129+12)

#print(Name_2_.index('1E 1743.1-2852'), Name_2_.index('AX J1745.6-2901'), Name_2_.index('CXO J212958.1+121002'))
#Name_2 = [e for e in Name_2_ if e not in ('1E 1743.1-2852', 'AX J1745.6-2901', 'CXO J212958.1+121002')]

print(Name_2_.index('1E 1743.1-2852'), Name_2_.index('AX J1745.6-2901'))
Name_2 = [e for e in Name_2_ if e not in ('1E 1743.1-2852', 'AX J1745.6-2901')]

#------------------------------------------------------------------------------------------------------------------------------------------------------------

AltName_R = Rit['AltName'].data.data
AltName_L = Liu['Name2'].data.data
AltName_L_ = Liu['Name3'].data.data

AltName_R1 = AltName_R
AltName_L1 = AltName_L
AltName_L_1 = AltName_L_

AltName_R1[np.where(AltName_R == '')] = 'KICK'
AltName_L1[np.where(AltName_L == '')] = 'LICK'
AltName_L_1[np.where(AltName_L_ == '')] = 'PICK'

AltName_R1 = list(AltName_R1)
AltName_L1 = list(AltName_L1)
AltName_L_1 = list(AltName_L_1)

#AltName_L1.pop(Name_2_.index('CXO J212958.1+121002'))
AltName_L1.pop(Name_2_.index('AX J1745.6-2901'))
AltName_L1.pop(Name_2_.index('1E 1743.1-2852'))

#AltName_L_1.pop(Name_2_.index('CXO J212958.1+121002'))
AltName_L_1.pop(Name_2_.index('AX J1745.6-2901'))
AltName_L_1.pop(Name_2_.index('1E 1743.1-2852'))

AltName_R = AltName_R1
AltName_L = AltName_L1
AltName_L_ = AltName_L_1

#------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------

alt_same_l = []
alt_same_r = []
for y in range(len(AltName_R)):
	for x in range(len(Name_2)):
		if AltName_R[y] == Name_2[x]: 
			alt_same_l.append(AltName_R[y])
			alt_same_r.append(Name_1[y])
		else:
			continue

#------------------------------------------------------------------------------------------------------------------------------------------------------------

alt_same_l_1 = []
alt_same_r_1 = []
for y in range(len(AltName_L)):
	for x in range(len(Name_1)):
		if AltName_L[y] == Name_1[x]:
			alt_same_r_1.append(Name_1[x])
			alt_same_l_1.append(Name_2[y])
		else:
			continue

#------------------------------------------------------------------------------------------------------------------------------------------------------------

alt_same_l_1_ = []
alt_same_r_1_ = []
for y in range(len(AltName_L_)):
	for x in range(len(Name_1)):
		if AltName_L_[y] == Name_1[x]:
			alt_same_r_1_.append(Name_1[x])
			alt_same_l_1_.append(Name_2[y])
		else:
			continue

#------------------------------------------------------------------------------------------------------------------------------------------------------------

alt_same_r_2 = []
alt_same_l_2 = []
for y in range(len(AltName_L)):
	for x in range(len(AltName_R)):
		if AltName_L[y] == AltName_R[x]:
			alt_same_l_2.append(Name_2[y])
			alt_same_r_2.append(Name_1[x])
		else:
			continue

#------------------------------------------------------------------------------------------------------------------------------------------------------------

alt_same_r_2_ = []
alt_same_l_2_ = []
for y in range(len(AltName_L_)):
	for x in range(len(AltName_R)):
		if AltName_L_[y] == AltName_R[x]:
			alt_same_l_2_.append(Name_2[y])
			alt_same_r_2_.append(Name_1[x])
		else:
			continue

#------------------------------------------------------------------------------------------------------------------------------------------------------------

dop_names_1 = []
dop_names_1_r = []
for y in range(len(Name_1)):
	for x in range(len(Simbad_name_list)):
		for z in range(len(Simbad_name_list[x].split("|"))):
			if Name_1[y] == Simbad_name_list[x].split("|")[z]:
				#print(Simbad_name_list[x])
				dop_names_1.append(Simbad_name_list[x])
				dop_names_1_r.append(Name_1[y])
			else:
				continue

dop_same_1 = []
dop_same_1_r = []
for y in range(len(Name_2)):
	for x in range(len(dop_names_1)):
		for z in range(len(dop_names_1[x].split("|"))):
			if Name_2[y] == dop_names_1[x].split("|")[z]:
				dop_same_1.append(Name_2[y])
				dop_same_1_r.append(dop_names_1_r[x])
			else:
				continue

#-----------------------------------------------------------------------------------------------------------------------------------------------------------

dop_alt_names_1 = []
dop_alt_names_r = []
for y in range(len(AltName_R)):
	for x in range(len(Simbad_name_list)):
		for z in range(len(Simbad_name_list[x].split("|"))):
			if AltName_R[y] == Simbad_name_list[x].split("|")[z]:
				dop_alt_names_1.append(Simbad_name_list[x])
				dop_alt_names_r.append(Name_1[y])
			else:
				continue

dop_alt_same_l = []
dop_alt_same_r = []
for y in range(len(Name_2)):
	for x in range(len(dop_alt_names_1)):
		for z in range(len(dop_alt_names_1[x].split("|"))):
			if Name_2[y] == dop_alt_names_1[x].split("|")[z]:
				dop_alt_same_l.append(Name_2[y])
				dop_alt_same_r.append(dop_alt_names_r[x])
			else:
				continue

#------------------------------------------------------------------------------------------------------------------------------------------------------------

dop_alt_names_2 = []
dop_alt_names_l = []
for y in range(len(AltName_L)):
	for x in range(len(Simbad_name_list)):
		for z in range(len(Simbad_name_list[x].split("|"))):
			if AltName_L[y] == Simbad_name_list[x].split("|")[z]:
				dop_alt_names_2.append(Simbad_name_list[x])
				dop_alt_names_l.append(Name_2[y])
			else:
				continue

dop_alt_same_r_ = []
dop_alt_same_l_ = []
for y in range(len(Name_1)):
	for x in range(len(dop_alt_names_2)):
		for z in range(len(dop_alt_names_2[x].split("|"))):
			if Name_1[y] == dop_alt_names_2[x].split("|")[z]:
				dop_alt_same_r_.append(Name_1[y])
				dop_alt_same_l_.append(dop_alt_names_l[x])
			else:
				continue

#------------------------------------------------------------------------------------------------------------------------------------------------------------

dop_alt_names_2_ = []
dop_alt_names_l_ = []
#Integral_name_list[x]=Integral_name_list[x].split('|')

#Names in Rit and Kolb
for y in range(len(AltName_L_)):
	for x in range(len(Simbad_name_list)):
		for z in range(len(Simbad_name_list[x].split("|"))):
			if AltName_L_[y] == Simbad_name_list[x].split("|")[z]:
				dop_alt_names_2_.append(Simbad_name_list[x])
				dop_alt_names_l_.append(Name_2[y])
			else:
				continue
		
dop_alt_same_r_1_ = []
dop_alt_same_l_1_ = []
for y in range(len(Name_1)):
	for x in range(len(dop_alt_names_2_)):
		for z in range(len(dop_alt_names_2_[x].split("|"))):
			if Name_1[y] == dop_alt_names_2_[x].split("|")[z]:
				dop_alt_same_r_1_.append(Name_1[y])
				dop_alt_same_l_1_.append(dop_alt_names_l_[x])
			else:
				continue

#------------------------------------------------------------------------------------------------------------------------------------------------------------

dop_names_2 = []
dop_names_2_l = []
for y in range(len(Name_2)):
	for x in range(len(Simbad_name_list)):
		for z in range(len(Simbad_name_list[x].split("|"))):
			if Name_2[y] == Simbad_name_list[x].split("|")[z]:  
				dop_names_2.append(Simbad_name_list[x])
				dop_names_2_l.append(Name_2[y])
			else:
				continue

dop_same_2 = []
dop_same_2_r = []
for y in range(len(Name_1)):
	for x in range(len(dop_names_2)):
		for z in range(len(dop_names_2[x].split("|"))):
			if Name_1[y] in dop_names_2[x].split("|")[z]:
				dop_same_2.append(dop_names_2_l[x])
				dop_same_2_r.append(Name_1[y])
			else:
				continue

#------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------


b = alt_same_l + alt_same_l_1 + alt_same_l_1_ + alt_same_l_2 + alt_same_l_2_ + dop_same_1 + dop_alt_same_l + dop_alt_same_l_ + dop_alt_same_l_1_ + dop_same_2
b_r = alt_same_r + alt_same_r_1 + alt_same_r_1_ + alt_same_r_2 + alt_same_r_2_ + dop_same_1_r + dop_alt_same_r + dop_alt_same_r_ + dop_alt_same_r_1_ + dop_same_2_r

def f(l):
	n = []
	for i in l:
		if i not in n:
			n.append(i)
	return n

final = f(b)
final_r = f(b_r)

L_n_R = []
R_n_L = []

for x in range(len(Name_2)):
	if Name_2[x] in final:
		continue
	else:
		L_n_R.append(Name_2[x])
		
for y in range(len(Name_1)):
	if Name_1[y] in final_r:
		continue
	else:
		R_n_L.append(Name_1[y])


def distance(nam):
	#print(nam)
	try:
		liu_tab = viz.query_constraints(catalog = 'J/A+A/469/807/lmxb', Name = nam)
		liu_coord = SkyCoord(ra = liu_tab[0]['RAJ2000'][0], dec = liu_tab[0]['DEJ2000'][0], unit = (u.hourangle, u.deg), frame='fk5')
		rit_tab = viz.query_region(liu_coord, catalog='B/cb/lmxbdata', radius=Angle(25, 'arcsec'))
		rit_coord = SkyCoord(ra = rit_tab[0]['RAJ2000'][0], dec = rit_tab[0]['DEJ2000'][0], unit = (u.hourangle, u.deg), frame='fk5') 
		d2d = liu_coord.separation(rit_coord)
		return nam, rit_tab[0]['Name'][0], d2d
	except IndexError:
		return  float('NaN'), float('NaN'), float('NaN')


nam_arr_ = [] 
nam_arr_1_ = []
d2d_arr_ = []
for name in L_n_R:
	nam, nam1, d2d = distance(name)
	nam_arr_.append(nam) 
	nam_arr_1_.append(nam1)
	d2d_arr_.append(d2d)

nam_arr = [x for x in nam_arr_ if pd.isnull(x) == False]
nam_arr_1 = [x for x in nam_arr_1_ if pd.isnull(x) == False]
d2d_arr = [x for x in d2d_arr_ if pd.isnull(x) == False]

#Not the same but close
#2E 1742.5-2858 J1745-2900 0d00m20.9963s
#EXO 1745-248 J1748-2446 0d00m05.6254s
print(nam_arr[7], nam_arr_1[7], d2d_arr[7])
print(nam_arr[10], nam_arr_1[10], d2d_arr[10])

nam_arr.pop(10)
nam_arr.pop(7)
nam_arr.append('Swift J1753.5-0127')
nam_arr.append('AX J1538.3-5541')
nam_arr.append('1A 1742-289')
#nam_arr.append('4U 2129+12')

nam_arr_1.pop(10)
nam_arr_1.pop(7)
nam_arr_1.append('J1753-0127')
nam_arr_1.append('J1538-5542')
nam_arr_1.append('J1745-2901')
#nam_arr_1.append('NGC 7078-X2')

The_Same_L = final + nam_arr
The_Same_R = final_r + nam_arr_1


LnR = []
RnL = []

for x in range(len(Name_2)):
	if Name_2[x] in The_Same_L:
		continue
	else:
		LnR.append(Name_2[x])
		
for y in range(len(Name_1)):
	if Name_1[y] in The_Same_R:
		continue
	else:
		RnL.append(Name_1[y])


print(len(LnR), len(The_Same_L), len(Name_2))
print(len(RnL), len(The_Same_R), len(Name_1))


#Extra = ['M 51 X-7', 'J0043+4112', 'J0043+4107', 'J0042+4118', 'J1242+3232', 'J0055-3738' , 'J1047+1234']
RnL_Galactic = [e for e in RnL if e not in ('M 51 X-7', 'J0043+4112', 'J0043+4107', 'J0042+4118', 'J1242+3232', 'J0055-3738' , 'J1047+1234')]
Cat_names = Name_2 + RnL_Galactic
Simb_names = list(names)

cat_not_sim = []
simb_lxb = []
simb_lxb_s = []	

No_simb_lxb = []
No_simb_lxb_s = []

for y in range(len(Cat_names)):
	if Cat_names[y] == 'CXO J212958.1+121002':
		print(Cat_names[y])
		cat_not_sim.append(Cat_names[y])
	else:
		try:
			tab = Simb_v2.query_object(Cat_names[y])
			if tab['MAIN_ID'][0] in Simb_names:
				simb_lxb.append(Cat_names[y])
				simb_lxb_s.append(tab['MAIN_ID'][0]) 
			else:
				No_simb_lxb.append(Cat_names[y])
				No_simb_lxb_s.append(tab['MAIN_ID'][0])	   
		except TypeError:
			cat_not_sim.append(Cat_names[y])


# # print([item for item, count in collections.Counter(simb_lxb_s_).items() if count > 1]) ['[SKM2002] 27', '[KRL2007b] 275', '4U 2129+12']
# # (AX J1745.6-2901, 1A 1742-289); (1E 1743.1-2852, SAX J1747.0-2853); (CXO J212958.1+121002, 4U 2129+12)
# # Первая и третья пары точно одно и то же

# simb_lxb_s = list(OrderedDict.fromkeys(simb_lxb_s_))
# simb_lxb = [e for e in simb_lxb_ if e not in ('AX J1745.6-2901', '1E 1743.1-2852', 'CXO J212958.1+121002')]

# Cat_names = [e for e in Cat_names_ if e not in ('AX J1745.6-2901', '1E 1743.1-2852', 'CXO J212958.1+121002')]
# Simb_names = list(names)

print(len(simb_lxb), len(simb_lxb_s), len(Cat_names), len(Simb_names))

def Simb_distance(nam):
	#print(nam)
	try:
		liu_tab = viz.query_constraints(catalog = 'J/A+A/469/807/lmxb', Name = nam)
		liu_coord = SkyCoord(ra = liu_tab[0]['RAJ2000'][0], dec = liu_tab[0]['DEJ2000'][0], unit = (u.hourangle, u.deg), frame='fk5')
		simb_tab = Simb_v2.query_region(liu_coord, radius = Angle(15, 'arcsec'))
		try:
			simb_coord = SkyCoord(ra = simb_tab['RA'][0], dec = simb_tab['DEC'][0], unit = (u.hourangle, u.deg), frame='icrs')
			d2d = simb_coord.separation(liu_coord)
			if simb_tab['MAIN_ID'][0] in Simb_names:
				return nam, simb_tab['MAIN_ID'][0], d2d, 'liu and simb lxb'
			else:
				return nam, simb_tab['MAIN_ID'][0], d2d, 'liu and no simb lxb'
		except TypeError:
				#print(nam,  float('NaN'), float('NaN'), '	 a')
				return nam, float('NaN'), float('NaN'), 'liu and no one around'
	except IndexError:
		rit_tab = viz.query_constraints(catalog = 'B/cb/lmxbdata', Name = nam)
		rit_coord = SkyCoord(ra = rit_tab[0]['RAJ2000'][0], dec = rit_tab[0]['DEJ2000'][0], unit = (u.hourangle, u.deg), frame='fk5')
		simb_tab = Simb_v2.query_region(rit_coord, radius = Angle(15, 'arcsec'))
		try:
			simb_coord = SkyCoord(ra = simb_tab['RA'][0], dec = simb_tab['DEC'][0], unit = (u.hourangle, u.deg), frame='icrs')
			d2d = simb_coord.separation(rit_coord)
			if simb_tab['MAIN_ID'][0] in Simb_names:
				#print(nam, simb_tab['MAIN_ID'][0], d2d, '	 3')
				return nam, simb_tab['MAIN_ID'][0], d2d, 'rit and simb lxb'
			else:
				#print(nam, simb_tab['MAIN_ID'][0], float('NaN'), '	 4')
				return nam, simb_tab['MAIN_ID'][0], d2d, 'rit and no simb lxb'
		except TypeError:
				#print(nam, float('NaN'), float('NaN'), '	 b')
				return nam, float('NaN'), float('NaN'), 'rit and no one around'

nam_s_arr_ = []
nam_s_arr_s_ = []
for name in cat_not_sim:
	n1, n2, _, _ = Simb_distance(name)
	nam_s_arr_.append(n1)
	nam_s_arr_s_.append(n2)

replacements = {'PSR J1417-4402':'1FGL J1417.7-4407', '[TDD2012] 2':'IGR J17480-2446', '[DCM2011] 8':'2XMM J174931.6-280805', '[KTR2013] V58':'[BPH2004] CX 1'}
replacer = replacements.get  # For faster gets.

nam_s_arr_s = [replacer(n, n) for n in nam_s_arr_s_]
nam_s_arr = nam_s_arr_

arr_same_obj = []
arr_n_same_obj = []

arr_same_obj_s = []
arr_n_same_obj_s = []

for y in arange(len(nam_s_arr_s)):
	if nam_s_arr_s[y] in Simb_names:
		arr_same_obj.append(nam_s_arr[y])
		arr_same_obj_s.append(nam_s_arr_s[y])
	else:
		arr_n_same_obj.append(nam_s_arr[y])
		arr_n_same_obj_s.append(nam_s_arr_s[y])

Sim_NoSimlxb = No_simb_lxb_s + arr_n_same_obj_s
Cat_NoSimlxb = No_simb_lxb + arr_n_same_obj

Sim_Simlxb = simb_lxb_s + arr_same_obj_s
Cat_Simlxb = simb_lxb + arr_same_obj

Simb_add = []
for y in range(len(Simb_names)):
	if Simb_names[y] in Sim_Simlxb:
		continue
	else:
		Simb_add.append(Simb_names[y])

pulsar = ['MAXI J0911-655', 'IGR J16597-3704', 'IGR J17379-3747', 'IGR J17591-2342']
All_names = Cat_names + Simb_add + pulsar

print(len(Sim_Simlxb), len(Cat_NoSimlxb), len(Cat_names), len(Simb_add), len(Simb_names), len(All_names))


# ------------------------------------------------------------------------------------------------------------------------------------------------------------


class LMXB:
	'''

	docstring for LMXB

	'''
	def __init__(self, name):

		 # Say my name
		self.name = name

		# SIMBAD functions
		self.SIMBAD_info = self.simb_get_info()
		self.SIMBAD_coord = self.simb_get_coord()
		self.Simbad_ra = self.simb_get_coord().ra.deg
		self.Simbad_dec = self.simb_get_coord().dec.deg
		self.SIMBAD_ids = self.simb_get_ids()
		self.SIMBAD_Vmag = self.simb_get_vmag()

		 # Сhecking GAIA and 2MASS for identifer
		self.GAIA_id_ = self.gaia_look_id() # Looking is there any id in gaia
		self.MASS_id_ = self.mass_look_id() # Looking is there any id in 2MASS

		 # Сoordinates and errors final getting 
		self.coord = self.get_coord()[0] 
		self.ra = self.coord.ra
		self.dec = self.coord.dec 
		self.coor_bib = self.get_coord()[1]
		self.coord_error = self.get_coord_error()[0]
		self.err_orig = self.get_coord_error()[1]
	
		 # Galactic coordinates
		self.gal_coord = self.gal_coord()
		self.lon = self.gal_coord.l
		self.lat = self.gal_coord.b

		 # GAIA
		self.GAIA_id = self.gaia_get_id()[0]
		self.Gmag = self.gaia_get_id()[1]
		self.distance = self.gaia_get_dist()

		self.CRSF = self.get_crsf()
		self.Spin = self.get_spin()
	

		self.Incl = self.catalog_get_incl()
		self.Lx_Lopt = self.get_lx_lopt()
		self.Mx_Mopt = self.get_mx_mopt()
		self.Mx = self.get_mx_mass()
		self.Mopt =  self.get_mopt_mass()

		self.b_v = self.get_b_v()
		self.eb_v = self.get_eb_v()

		 # Fluxes
		self.XMM_flux = self.xmm_get_flux()
		self.XRT_flux = self.xrt_get_flux()
		self.BAT_flux = self.bat_get_flux()
		self.INTEGRAL_flux = self.integral_get_flux()
		self.Chandra = self.chandra_get_flux()
		self.Chandra_min_flux = self.chandra_get_flux()[0]
		self.Chandra_max_flux = self.chandra_get_flux()[1]

		self.MASS_flux = self.mass_get_flux()
		self.GAIA_flux = self.gaia_get_flux()

		
		# Parameters from catalogs

	
		self.alt_name = self.catalog_alt_name()
		self.Opt_name = self.catalog_opt_name()
		self.Xtype = self.catalog_get_type()
		self.Porb_info = self.catalog_porb_info()
		self.Porb = self.catalog_get_porb()
		self.Ppulse = self.catalog_get_pulse()
		self.Pos = self.catalog_get_pos()


		self.gaia_all = self.gaia_get_info()
		self.BP_RP = self.gaia_all[0]
		self.BP_G = self.gaia_all[1]
		self.Teff = self.gaia_all[2]
		self.b_Teff = self.gaia_all[3]
		self.B_Teff = self.gaia_all[4]
		self.AG = self.gaia_all[5]
		self.b_AG = self.gaia_all[6]
		self.B_AG = self.gaia_all[7]
		self.EBP_RP = self.gaia_all[8]
		self.b_EBP_RP = self.gaia_all[9]
		self.B_EBP_RP = self.gaia_all[10]
		self.Lum = self.gaia_all[11]
		self.b_Lum = self.gaia_all[12]
		self.B_Lum = self.gaia_all[13]

		self.Plx = self.gaia_get_parallax()[0]
		self.e_Plx = self.gaia_get_parallax()[1]

		self.Gmag1 = self.get_gmag()
		self.comments = self.get_comments()
		self.Starhorse = self.starhouse_get_info()

		self.Chandra_flux = self.chandra_get_flux()[2]

		self.cross_check_DR2_1 = self.cross_check_dr2_1()
		self.cross_check_DR2_2 = self.cross_check_dr2_2()
		self.cross_check_DR2_3 = self.cross_check_dr2_3()
		self.cross_check_DR3_1 = self.cross_check_dr3_1()
		self.cross_check_DR3_2 = self.cross_check_dr3_2()
		self.cross_check_DR3_3 = self.cross_check_dr3_3()
		#self.cross_check = self.final_cross_check() 
		

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# SIMBAD

	def simb_get_info(self):
		'''
		Getting all info in table format
		'''
		try:
			tabl_info = Simb_v2.query_object(self.name)
			ra = tabl_info['RA'][0]
			return tabl_info
		except TypeError:
			if self.name in nam_s_arr:
				#print(name, nam_s_arr_s[nam_s_arr.index(name)])
				tabl_info = Simb_v2.query_object(nam_s_arr_s[nam_s_arr.index(self.name)])
				return tabl_info 
			elif self.name in The_Same_L:
				#print(name, nam_s_arr_s[nam_s_arr.index(name)])
				tabl_info = Simb_v2.query_object(The_Same_R[The_Same_L.index(self.name)])
				return tabl_info 
			else:
				print('Something is wrong')
				return float('NaN')

	def simb_get_coord(self):
		'''
		Getting Simbad coordinates
		'''
		Simb_RA = self.SIMBAD_info['RA'][0]
		Simb_DEC = self.SIMBAD_info['DEC'][0]
		#print(Simb_RA)
		Simb_coord = SkyCoord(ra = Simb_RA, dec = Simb_DEC, unit=(u.hourangle, u.deg))
		return Simb_coord

	def simb_get_ids(self):
		'''
		Getting diffrent identifiers for same object and making a list of it
		'''
		Sim_inf = self.SIMBAD_info
		IDS = Sim_inf['IDS'][0]
		#IDS = Sim_inf['IDS'][0].decode('utf-8')
		splited_ids = IDS.split('|')
		lis_ids = [x for x in splited_ids]
		return lis_ids 


	def simb_get_vmag(self):
		'''
		Getting diffrent identifiers for same object and making a list of it
		'''
		Sim_inf = self.SIMBAD_info
		if ma.is_masked(Sim_inf["FLUX_V"][0]) == False:
			return Sim_inf["FLUX_V"][0]
		else:
			return float('NaN')
	
# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# Gaia and 2MASS

	def gaia_look_id(self):
		'''
		Checking list of identifiers for Gaia's one
		'''
		for x in self.SIMBAD_ids:
			if 'Gaia DR2' in x:
				return int(x.split()[2])
			else:
				continue

	def mass_look_id(self):
		'''
		Checking list of identifiers for 2MASS's one
		'''
		for x in self.SIMBAD_ids:
			if '2MASS' in x:
				return x.split('J')[1]
			else:
				continue

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# Chandra

	def chandra_look_coord(self):
		'''
		Looking for Chandra observations and coordinates
		'''
		Chandra_v2 = 'IX/57/csc2master' #newest version
		Chan_cat = viz(catalog = Chandra_v2, columns = ['2CXO', 'RAICRS', 'DEICRS', 'r0', 'r1', 'Fluxw', 'Fluxb'])
		try:
			Chan_tab = Chan_cat.query_object(self.name, radius = Angle(15, 'arcsec'))
			Chan_RA = Chan_tab[0]['RAICRS'][0]
		except (TypeError, IndexError, KeyError):
			if self.name in nam_s_arr:
				Chan_tab = Chan_cat.query_object(nam_s_arr_s[nam_s_arr.index(self.name)], radius = Angle(15, 'arcsec'))
			elif self.name in The_Same_L:
				Chan_tab = Chan_cat.query_object(The_Same_R[The_Same_L.index(self.name)], radius = Angle(15, 'arcsec'))
		return Chan_tab 

	def chandra_get_coord(self):
		'''
		Getting Chandra's data if there are any
		'''
		Chandra_v2 = 'IX/57/csc2master' #newest version
		Chan_cat = viz(catalog = Chandra_v2, columns = ['2CXO', 'RAICRS', 'DEICRS', 'r0', 'r1', 'Fluxw', 'Fluxb'])
		try:
			Chan_tab = Chan_cat.query_object(self.name, radius = Angle(15, 'arcsec'))
			Chan_RA = Chan_tab[0]['RAICRS'][0] #RAICRS
			Chan_DEC = Chan_tab[0]['DEICRS'][0] #DEICRS
			Chan_coord = SkyCoord(ra = Chan_RA, dec = Chan_DEC, unit = (u.deg, u.deg))
			r0 = Chan_tab[0]['r0'][0]
			r1 = Chan_tab[0]['r1'][0]
			r_ = m.sqrt((r0)**2+(r1)**2) 
			r = Angle(r_, 'arcsec').degree
		except (TypeError, IndexError, KeyError):
			if self.name in nam_s_arr:
				Chan_tab = Chan_cat.query_object(nam_s_arr_s[nam_s_arr.index(self.name)], radius = Angle(15, 'arcsec'))
				Chan_RA = Chan_tab[0]['RAICRS'][0] #RAICRS
				Chan_DEC = Chan_tab[0]['DEICRS'][0] #DEICRS
				Chan_coord = SkyCoord(ra = Chan_RA, dec = Chan_DEC, unit = (u.deg, u.deg))
				r0 = Chan_tab[0]['r0'][0]
				r1 = Chan_tab[0]['r1'][0]
				r_ = m.sqrt((r0)**2+(r1)**2) 
				r = Angle(r_, 'arcsec').degree
			elif self.name in The_Same_L:
				Chan_tab = Chan_cat.query_object(The_Same_R[The_Same_L.index(self.name)], radius = Angle(15, 'arcsec'))
				Chan_RA = Chan_tab[0]['RAICRS'][0] #RAICRS
				Chan_DEC = Chan_tab[0]['DEICRS'][0] #DEICRS
				Chan_coord = SkyCoord(ra = Chan_RA, dec = Chan_DEC, unit = (u.deg, u.deg))
				r0 = Chan_tab[0]['r0'][0]
				r1 = Chan_tab[0]['r1'][0]
				r_ = m.sqrt((r0)**2+(r1)**2) 
				r = Angle(r_, 'arcsec').degree
		return Chan_coord, r

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# XMM

	def xmm_get_coord(self):
		'''
		Getting XMM's data if there are any
		'''
		try:
			XMM = heasarc.query_object(self.name, mission='xmmssc', radius = Angle(15,'arcsec'), fields = 'NAME, RA, DEC, ERROR_RADIUS')	
			XMM_coord = SkyCoord(ra = XMM['RA'][0], dec = XMM['DEC'][0], unit = (u.deg, u.deg))
			XMM_err = XMM['ERROR_RADIUS'][0]
			return XMM_coord, XMM_err
		except (astroquery.exceptions.InvalidQueryError, requests.exceptions.SSLError):
			if self.name in nam_s_arr:
				XMM = heasarc.query_object(nam_s_arr_s[nam_s_arr.index(self.name)], mission='xmmssc', radius = Angle(15,'arcsec'), fields = 'NAME, RA, DEC, ERROR_RADIUS')
				XMM_coord = SkyCoord(ra = XMM['RA'][0], dec = XMM['DEC'][0], unit = (u.deg, u.deg))
				XMM_err = XMM['ERROR_RADIUS'][0]
				return XMM_coord, XMM_err
			elif self.name in The_Same_L:
				XMM = heasarc.query_object(The_Same_R[The_Same_L.index(self.name)], mission='xmmssc', radius = Angle(15,'arcsec'), fields = 'NAME, RA, DEC, ERROR_RADIUS')
				XMM_coord = SkyCoord(ra = XMM['RA'][0], dec = XMM['DEC'][0], unit = (u.deg, u.deg))
				XMM_err = XMM['ERROR_RADIUS'][0]	
				return XMM_coord, XMM_err


# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# Swift BAT

	def swift_look_coord(self):
		'''
		Looking for Swift XRT observations using simbad coordinates
		'''
		XRT = 'IX/58/2sxps' # https://arxiv.org/abs/1911.11710 (An improved and expanded Swift X-ray telescope point source catalog)
		XRT_cat = viz(columns = ['RAJ2000', 'DEJ2000',  'Err90', '+_r', 'Det', '2SXPS'], column_filters={'Det': '=0'})
		try:
			XRT_tab = XRT_cat.query_object(self.name, catalog = XRT, radius = Angle(15, 'arcsec'))
			XRT_RA = XRT_tab[0]['RAJ2000'][0]
		except (TypeError, IndexError, KeyError):
			if self.name in nam_s_arr:
				XRT_tab = XRT_cat.query_object(nam_s_arr_s[nam_s_arr.index(self.name)], catalog = XRT, radius = Angle(15, 'arcsec'))
			elif self.name in The_Same_L:
				XRT_tab = XRT_cat.query_object(The_Same_R[The_Same_L.index(self.name)], catalog = XRT, radius = Angle(15, 'arcsec'))
		#print(XRT_tab)
		return XRT_tab

	def swift_get_coord(self):
		'''
		Getting Swift XRT's data if there are any
		'''
		XRT = 'IX/58/2sxps' # https://arxiv.org/abs/1911.11710 (An improved and expanded Swift X-ray telescope point source catalog)
		XRT_cat = viz(columns = ['RAJ2000', 'DEJ2000',  'Err90', '+_r', 'Det', '2SXPS'], column_filters={'Det': '=0'})
		try:
			XRT_tab = XRT_cat.query_object(self.name, catalog = XRT, radius = Angle(15, 'arcsec'))
			XRT_name = XRT_tab[0]['_2SXPS'][0]
			XRT_coord = SkyCoord(ra = XRT_tab[0]['RAJ2000'][0], dec = XRT_tab[0]['DEJ2000'][0], unit = (u.deg, u.deg))
			XRT_err = XRT_tab[0]['Err90'][0]
			XRT_r = XRT_tab[0]['_r'][0]
		except (TypeError, IndexError, KeyError):
			if self.name in nam_s_arr:
				XRT_tab = XRT_cat.query_object(nam_s_arr_s[nam_s_arr.index(self.name)], catalog = XRT, radius = Angle(15, 'arcsec'))
				XRT_name = XRT_tab[0]['_2SXPS'][0]
				XRT_coord = SkyCoord(ra = XRT_tab[0]['RAJ2000'][0], dec = XRT_tab[0]['DEJ2000'][0], unit = (u.deg, u.deg))
				XRT_err = XRT_tab[0]['Err90'][0]
				XRT_r = XRT_tab[0]['_r'][0]
			elif self.name in The_Same_L:
				XRT_tab = XRT_cat.query_object(The_Same_R[The_Same_L.index(self.name)], catalog = XRT, radius = Angle(15, 'arcsec'))
				XRT_name = XRT_tab[0]['_2SXPS'][0]
				XRT_coord = SkyCoord(ra = XRT_tab[0]['RAJ2000'][0], dec = XRT_tab[0]['DEJ2000'][0], unit = (u.deg, u.deg))
				XRT_err = XRT_tab[0]['Err90'][0]
				XRT_r = XRT_tab[0]['_r'][0]
		return XRT_coord, XRT_err, XRT_r, XRT_name

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# Gettig Gaia coordinates
	
	def gaia_get_coord(self):
		'''
		Gets gaia coordinates from gaia id

		'''
		try:
			gaia_tab = viz(catalog = 'I/350/gaiaedr3', columns = ['RA_ICRS', 'DE_ICRS']).query_constraints(Source = self.GAIA_id_)
			gaia_coord = SkyCoord(ra = gaia_tab[0]['RA_ICRS'][0], dec = gaia_tab[0]['DE_ICRS'][0], unit = (u.deg, u.deg))
		except IndexError:
			try:
				gaia_tab = viz(catalog = 'I/350/gaiaedr3', columns = ['RA_ICRS', 'DE_ICRS', 'Source']).query_object(self.name, radius = Angle(1,'arcsec'))
				gaia_coord = SkyCoord(ra = gaia_tab[0]['RA_ICRS'][0], dec = gaia_tab[0]['DE_ICRS'][0], unit = (u.deg, u.deg))
			except (TypeError, IndexError, KeyError):
				if self.name in nam_s_arr:
					gaia_tab = viz(catalog = 'I/350/gaiaedr3', columns = ['RA_ICRS', 'DE_ICRS', 'Source']).query_object(nam_s_arr_s[nam_s_arr.index(self.name)], radius = Angle(1,'arcsec'))
					gaia_coord = SkyCoord(ra = gaia_tab[0]['RA_ICRS'][0], dec = gaia_tab[0]['DE_ICRS'][0], unit = (u.deg, u.deg))
				elif self.name in The_Same_L:
					gaia_tab = viz(catalog = 'I/350/gaiaedr3', columns = ['RA_ICRS', 'DE_ICRS', 'Source']).query_object(The_Same_R[The_Same_L.index(self.name)], radius = Angle(1,'arcsec'))
					gaia_coord = SkyCoord(ra = gaia_tab[0]['RA_ICRS'][0], dec = gaia_tab[0]['DE_ICRS'][0], unit = (u.deg, u.deg))
		return gaia_coord	
			


	def gaia_get_err(self):
		'''
		Gets Position error for Gaia

		'''
		try:
			gaia_tab = viz(catalog = 'I/350/gaiaedr3', columns = ['e_RA_ICRS','e_DE_ICRS']).query_constraints(Source = self.GAIA_id_)
			err_1 = gaia_tab[0]['e_RA_ICRS'][0]
			err_2 = gaia_tab[0]['e_DE_ICRS'][0]
			sum_err = ( (err_1)**2 + (err_2)**2 )**(1/2)
		except IndexError:
			try:
				gaia_tab = viz(catalog='I/350/gaiaedr3',columns=['e_RA_ICRS','e_DE_ICRS','Source']).query_object(self.name, radius=Angle(1,'arcsec'))
				err_1 = gaia_tab[0]['e_RA_ICRS'][0]
				err_2 = gaia_tab[0]['e_DE_ICRS'][0]
				sum_err = ( (err_1)**2 + (err_2)**2 )**(1/2)
			except (TypeError, IndexError, KeyError):
				if self.name in nam_s_arr:
					gaia_tab = viz(catalog='I/350/gaiaedr3',columns=['e_RA_ICRS','e_DE_ICRS','Source']).query_object(nam_s_arr_s[nam_s_arr.index(self.name)], radius=Angle(1,'arcsec'))
					err_1 = gaia_tab[0]['e_RA_ICRS'][0]
					err_2 = gaia_tab[0]['e_DE_ICRS'][0]
					sum_err = ( (err_1)**2 + (err_2)**2 )**(1/2)
				elif self.name in The_Same_L:
					gaia_tab = viz(catalog='I/350/gaiaedr3',columns=['e_RA_ICRS','e_DE_ICRS','Source']).query_object(The_Same_R[The_Same_L.index(self.name)], radius=Angle(1,'arcsec'))
					err_1 = gaia_tab[0]['e_RA_ICRS'][0]
					err_2 = gaia_tab[0]['e_DE_ICRS'][0]
					sum_err = ( (err_1)**2 + (err_2)**2 )**(1/2)
		return sum_err
	


	def gaia_get_coord_vmag(self):
		''' 
		Gets Gaia coordinates from comparing Vmag and Gmag 

		''' 
		try:
			#all gaia observation around coordinates with search radius 10 arcsec
			gaia  = viz(columns=['Source','RA_ICRS','DE_ICRS','Gmag','e_RA_ICRS','e_DE_ICRS'])
			gaia_tab =  gaia.query_region(self.SIMBAD_coord, catalog = 'I/350/gaiaedr3', radius = Angle(10, 'arcsec'))

			#Checks if Gmag is in tolerance range of Vmag  true = in range, false = out of range 
			check_Vmag = [m.isclose(self.SIMBAD_Vmag, gaia_tab[0][x]['Gmag'], abs_tol=2) for x in range(len(gaia_tab[0]))]

			#Gets Gaia observations of which are in tolerance range
			checked_gaia_id = [gaia_tab[0][x] for x in range(len(check_Vmag)) if check_Vmag[x] == True]
			
			#subtract Vmag from Gmag of obersevations in tolerance range 
			final_gaia = [abs(checked_gaia_id[x]['Gmag'] - self.SIMBAD_Vmag) for x in range(len(checked_gaia_id))]

			#outputs Gaia observation with the smallest differance of V- and Gmag
			final_gaia_id = checked_gaia_id[final_gaia.index(min(final_gaia))]
			gaia_coord = SkyCoord(ra = final_gaia_id['RA_ICRS'],dec = final_gaia_id['DE_ICRS'], unit = (u.deg, u.deg))

			return gaia_coord
		except (TypeError, IndexError, ValueError): 
			return 'VmagError'



	def gaia_get_err_Vmag(self):
		try:
			#all gaia observation around coordinates with search radius 10 arcsec
			gaia  = viz(columns=['Source','RA_ICRS','DE_ICRS','Gmag','e_RA_ICRS','e_DE_ICRS'])
			gaia_tab =  gaia.query_region(self.SIMBAD_coord, catalog = 'I/350/gaiaedr3', radius = Angle(10, 'arcsec'))

			#Checks if Gmag is in tolerance range of Vmag  true = in range, false = out of range 
			check_Vmag = [m.isclose(self.SIMBAD_Vmag, gaia_tab[0][x]['Gmag'], abs_tol=2) for x in range(len(gaia_tab[0]))]

			#Gets Gaia observations of which are in tolerance range
			checked_gaia_id = [gaia_tab[0][x] for x in range(len(check_Vmag)) if check_Vmag[x]== True]
			
			#subtract Vmag from Gmag of obersevations in tolerance range 
			final_gaia = [abs(checked_gaia_id[x]['Gmag'] - self.SIMBAD_Vmag) for x in range(len(checked_gaia_id))]

			#outputs Gaia observation with the smallest differance of V- and Gmag
			final_gaia_id = checked_gaia_id[final_gaia.index(min(final_gaia))]
			gaia_coord = SkyCoord(ra = final_gaia_id['RA_ICRS'],dec = final_gaia_id['DE_ICRS'], unit = (u.deg, u.deg))
			err_1 = final_gaia_id['e_RA_ICRS']
			err_2 = final_gaia_id['e_DE_ICRS']
			sum_err = ( (err_1)**2 + (err_2)**2 )**(1/2)
			return sum_err
		except (IndexError,ValueError): 
			return 'VmagError'	

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# Gettig 2MASS coordinates (or Gaia)

	def mass_get_coord(self):
		'''
		Use Simbad Coordinates to get 2Mass coordinates and check if observation has the right 2Mass id

		'''
		Mass_tab = viz(catalog = 'II/246/out', column_filters = {'2MASS': self.MASS_id_}).query_region(self.SIMBAD_coord, radius = Angle(10, 'arcsec'))
		Mass_coord = SkyCoord(ra = Mass_tab[0]['RAJ2000'][0], dec = Mass_tab[0]['DEJ2000'][0], unit = (u.deg, u.deg)) #unit = (u.hourangle, u.deg)
		return Mass_coord

	def mass_get_err(self):
		'''
		Use Simbad Coordinates to get 2Mass coordinates errors

		'''
		Mass_tab =  viz(catalog = 'II/246/out', columns=['2MASS','errMaj','errMin'], column_filters = {'2MASS': self.MASS_id_}).query_region(self.SIMBAD_coord, radius = Angle(10, 'arcsec'))
		err_1 = Mass_tab[0]['errMaj'][0]
		err_2 = Mass_tab[0]['errMin'][0]
		sum_err = ( (err_1)**2 + (err_2)**2 )**(1/2)
		return sum_err


	def gaia_get_coord_mass(self):
		'''
		Looks for observation around 2Mass coordinates with 1 arcsec radius in
		Gaia Catalog, outputs gaia coordinates if only one observation is available 
		
		'''
		gaia = viz(columns=['RA_ICRS', 'DE_ICRS', 'Source', 'Gmag'])
		gaia_tab = gaia.query_region(self.mass_get_coord(),catalog='I/350/gaiaedr3',radius=Angle(1,'arcsec'))
		try:
			if len(gaia_tab[0]) == 1:
				gaia_coord = SkyCoord(ra = gaia_tab[0]['RA_ICRS'][0], dec = gaia_tab[0]['DE_ICRS'][0], unit = (u.deg, u.deg))
				return gaia_coord
			else: 
				return 'None'
		except IndexError:
			return 'None'
			

	def gaia_get_err_mass(self):
		'''
		Gets Position error for Gaia from 2Mass coordinates 
		
		'''
		gaia = viz(columns=['e_RA_ICRS', 'e_DE_ICRS'])
		gaia_tab = gaia.query_region(self.gaia_get_coord_mass(), catalog='I/350/gaiaedr3', radius=Angle(1,'arcsec'))
		err_1 = gaia_tab[0]['e_RA_ICRS'][0]
		err_2 = gaia_tab[0]['e_DE_ICRS'][0]
		sum_err = ( (err_1)**2 + (err_2)**2 )**(1/2)
		return sum_err

#------------------------------------------------------------------------------------------------------------------------------------------------------------
#Gettig catalog coordinates


	def catalog_get_coord(self):
		''' 
		info from Ritter and Kolb 2004 (Rit_cat = 'B/cb/lmxbdata')
		https://ui.adsabs.harvard.edu/abs/2003A%26A...404..301R
		or 
		from Liu, van Paradijs and van den Heuvel 2008 (Liu_cat = 'J/A+A/469/807/lmxb')
		https://ui.adsabs.harvard.edu/abs/2007A%26A...469..807L

		'''
		try:
			cat_tab = viz.query_constraints(catalog = Rit_cat, Name = self.name)
			cat_coord = SkyCoord(ra = cat_tab[0]['RAJ2000'][0], dec = cat_tab[0]['DEJ2000'][0], unit = (u.hourangle, u.deg))
		except IndexError:
			cat_tab = viz.query_constraints(catalog = Liu_cat, Name = self.name)
			cat_coord = SkyCoord(ra = cat_tab[0]['RAJ2000'][0], dec = cat_tab[0]['DEJ2000'][0], unit = (u.hourangle, u.deg))
		return cat_coord
		

#------------------------------------------------------------------------------------------------------------------------------------------------------------
#Gettig coordinates. Main function


	def get_coord(self):
		print(self.name)
		if self.GAIA_id_ != None:
			print('GAIA')
			return self.gaia_get_coord(), '2020yCat.1350....0G'
		elif self.MASS_id_ != None:
			if self.gaia_get_coord_mass() is not 'None':
				print('GAIA from 2MASS')
				return self.gaia_get_coord_mass(), '2020yCat.1350....0G'
			else:
				print('2MASS')
				return self.mass_get_coord(), '2003yCat.2246....0C'
		elif ma.is_masked(self.SIMBAD_Vmag) == False and self.gaia_get_coord_vmag() is not 'VmagError':
			print('GAIA Vmag')
			return self.gaia_get_coord_vmag(), '2020yCat.1350....0G'
		else:
			if self.chandra_look_coord() != []:
				print('Chandra')
				return self.chandra_get_coord()[0],'2010ApJS..189...37E'
			else:
				try:
					XMM_coord = self.xmm_get_coord()[0]
					print('XMM')
					return XMM_coord, '2020A&A...641A.136W'
				except (TypeError, IndexError, KeyError):
					if self.swift_look_coord() != []:
						print('SWIFT BAT')
						return self.swift_get_coord()[0],'2020yCat.9058....0E'
					else:
						try:
							cat_coord = self.catalog_get_coord()
							print('CATALOG')
							return cat_coord, '2007A&A...469..807L'
						except IndexError:
							print('SIMBAD')
							return self.SIMBAD_coord, 'Simbad Nothing'

	def get_coord_error(self):
		'''
		Get corresponding Position errors 
		
		'''
		if self.GAIA_id_ != None:
			return self.gaia_get_err(), 'Gaia'
		elif self.MASS_id_ != None:
			if self.gaia_get_coord_mass() is not 'None':
				return self.gaia_get_err_mass(), 'Gaia'
			else:
				return self.mass_get_err(), '2MASS'
		elif ma.is_masked(self.SIMBAD_Vmag) == False and self.gaia_get_coord_vmag() is not 'VmagError':
			return self.gaia_get_err_Vmag(), 'Gaia'
		else:
			if self.chandra_look_coord() != []:
				return self.chandra_get_coord()[1],'Chandra'	
			else:
				try:
					XMM_error = self.xmm_get_coord()[1]
					return XMM_error, 'XMM'
				except (TypeError, IndexError, KeyError):					
					if self.swift_look_coord() != []:
						return self.swift_get_coord()[1], 'SWIFT XRT'
					else:
						return float('NaN'),'Catalog' 


	def gaia_get_id(self):
		''' 
		Gets gaia id and Gmag for source
	
		'''
		if self.GAIA_id_ != None:
			return self.GAIA_id_, float('NaN')

		elif self.MASS_id_ != None:
			if self.gaia_get_coord_mass() is not 'None':
				
				gaia = viz(columns=['Source'])
				gaia_tab = gaia.query_region(self.gaia_get_coord_mass(), catalog='I/350/gaiaedr3', radius=Angle(1, 'arcsec'))
				gaia_id = gaia_tab[0]['Source'][0]
				return gaia_id, float('NaN')
			else:
				try:

					coord = SkyCoord(ra  = self.ra, dec = self.dec, unit= (u.deg, u.deg))
					gaia = viz(columns = ['Source','RA_ICRS','DE_ICRS','Gmag','e_RA_ICRS','e_DE_ICRS'])

					#all gaia observation around coordinates with search radius 10 arcsec
					gaia_tab = gaia.query_region(coord, catalog = 'I/350/gaiaedr3', radius = Angle(1,'arcsec'))
					
					#Checks if Gmag is in tolerance range of Vmag  true= in range, false= out of range 
					check_Vmag = [m.isclose(self.SIMBAD_Vmag, gaia_tab[0][x]['Gmag'], abs_tol = 2) for x in range(len(gaia_tab[0]))]
					
					#Gets gaia observations of which are in tolerance range
					checked_gaia_id = [gaia_tab[0][x] for x in range(len(check_Vmag)) if check_Vmag[x] == True]
					
					#subtract Vmag from Gmag of obersevations in tolerance range 
					final_gaia = [abs(checked_gaia_id[x]['Gmag'] - self.SIMBAD_Vmag) for x in range(len(checked_gaia_id))]
					
					#outputs gaia observation with the smallest differance of V- and Gmag
					final_gaia_id = checked_gaia_id[final_gaia.index(min(final_gaia))]
					

					gaia_id = final_gaia_id['Source']
					Gmag = final_gaia_id['Gmag']
	
					return gaia_id, Gmag
				except (IndexError,ValueError):
					return float('NaN'), float('NaN')
		else:
			try:
				if not m.isnan(self.SIMBAD_Vmag):  

					gaia = viz(columns=['Source','RA_ICRS','DE_ICRS','Gmag','e_RA_ICRS','e_DE_ICRS'])

					#all gaia observation around coordinates with search radius 10 arcsec
					gaia_tab = gaia.query_region(self.SIMBAD_coord, catalog='I/350/gaiaedr3', radius=Angle(10,'arcsec'))

					#Checks if Gmag is in tolerance range of Vmag  true= in range, false= out of range 
					check_Vmag = [m.isclose(self.SIMBAD_Vmag, gaia_tab[0][x]['Gmag'], abs_tol=2) for x in range(len(gaia_tab[0]))]

					#Gets Gaia observations of which are in tolerance range
					checked_gaia_id = [gaia_tab[0][x] for x in range(len(check_Vmag)) if check_Vmag[x] == True]

					#subtract Vmag from Gmag of obersevations in tolerance range 
					final_gaia = [abs(checked_gaia_id[x]['Gmag'] - self.SIMBAD_Vmag) for x in range(len(checked_gaia_id))]

					#outputs Gaia observation with the smallest differance of V- and Gmag
					final_gaia_id = checked_gaia_id[final_gaia.index(min(final_gaia))]
		
					gaia_id = final_gaia_id['Source']
					Gmag = final_gaia_id['Gmag']

					return gaia_id, Gmag
				else: 
					return float('NaN'), float('NaN')
			except (TypeError, IndexError, ValueError):
				return float('NaN'), float('NaN')		


	def gaia_get_dist(self):
		'''
		Uses Gaia identifier to get gaia distance 
		
		'''
		coord = self.coord
		if not m.isnan(self.GAIA_id):

			try:
				query = 'select top 10 s.source_id , s.r_med_geo,s.r_lo_geo, s.r_hi_geo,s.r_med_photogeo,s.r_lo_photogeo,s.r_hi_photogeo from external.gaiaedr3_distance as s where s.source_id ={}'.format(self.GAIA_id)
				job = Gaia.launch_job(query)
				r = job.get_results()
				med_geo = r['r_med_geo'][0]
				lo_geo = r['r_lo_geo'][0]
				hi_geo = r['r_hi_geo'][0]
				med_photogeo = r['r_med_photogeo'][0]
				lo_photogeo = r['r_lo_photogeo'][0]
				hi_photogeo = r['r_hi_photogeo'][0]  
	
				return med_geo, lo_geo, hi_geo, med_photogeo, lo_photogeo, hi_photogeo
			except IndexError:
				return float('Nan'), float('NaN'), float('NaN'), float('Nan'), float('NaN'), float('NaN')
		else: 
			return float('Nan'), float('NaN'), float('NaN'), float('Nan'), float('NaN'), float('NaN')


	def gaia_get_info(self):
		'''
		Use Gaia-ID to get Information from Gaia DR2
		
		'''
		if not m.isnan(self.GAIA_id):
			try:
				gaia = viz(catalog = 'I/345/gaia2', columns = ['BP-RP', 'BP-G', 'Teff', 'b_Teff', 'B_Teff', 'AG', 'b_AG', 'B_AG', 'E(BP-RP)', 'b_E(BP-RP)', 'B_E(BP-RP)', 'Lum', 'b_Lum', 'B_Lum'])
				gaia_tab = gaia.query_constraints(Source = self.GAIA_id)
				#print(self.GAIA_id, gaia_tab)
				BP_RP = gaia_tab[0]['BP-RP'][0]
				BP_G = gaia_tab[0]['BP-G'][0]
				Teff = gaia_tab[0]['Teff'][0]
				b_Teff = gaia_tab[0]['b_Teff'][0]
				B_Teff = gaia_tab[0]['B_Teff'][0]
				AG = gaia_tab[0]['AG'][0]
				b_AG = gaia_tab[0]['b_AG'][0]
				B_AG = gaia_tab[0]['B_AG'][0]
				EBP_RP = gaia_tab[0]['E_BP-RP_'][0]
				b_EBP_RP = gaia_tab[0]['b_E_BP-RP_'][0]
				B_EBP_RP = gaia_tab[0]['B_E_BP-RP_'][0]
				Lum = gaia_tab[0]['Lum'][0]
				b_Lum = gaia_tab[0]['b_Lum'][0]
				B_Lum = gaia_tab[0]['B_Lum'][0]

				return BP_RP, BP_G, Teff, b_Teff, B_Teff, AG, b_AG, B_AG, EBP_RP, b_EBP_RP, B_EBP_RP, Lum, b_Lum, B_Lum
			except IndexError:
				#print(self.name, self.GAIA_id, gaia_tab)
				return float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN')
		else:
			return float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN')
	
	def gaia_get_parallax(self):
		if not m.isnan(self.GAIA_id):
			try:
				gaia = viz(catalog = 'I/350/gaiaedr3', columns = ['Plx','e_Plx'])	
				gaia_tab = gaia.query_constraints(Source = self.GAIA_id)
				Plx = gaia_tab[0]['Plx'][0]
				e_Plx = gaia_tab[0]['e_Plx'][0]
				return Plx, e_Plx
			except IndexError:
				return float('NaN'),float('NaN')	
		else:
			return float('NaN'),float('NaN')	


	def get_gmag(self):
		if not m.isnan(self.GAIA_id):
			try:
				gaia = viz(catalog = 'I/350/gaiaedr3', columns = ['Source, Gmag, e_Gmag'])
				gaia_tab = gaia.query_constraints(Source = self.GAIA_id)
				Gmag = gaia_tab[0]['Gmag'][0]
				e_Gmag = gaia_tab[0]['e_Gmag'][0]
				return Gmag, e_Gmag
			except IndexError:
				try:
					gaia = viz(catalog='I/345/gaia2', columns = ['Source, Gmag, e_Gmag'])
					gaia_tab = gaia.query_constraints(Source = self.GAIA_id)
					Gmag = gaia_tab[0]['Gmag'][0]
					e_Gmag = gaia_tab[0]['e_Gmag'][0]
					return Gmag, e_Gmag
				except IndexError:
					return float('NaN'), float('NaN')
		else:
			return float('NaN'),float('NaN')


 
	def gal_coord(self):
		'''
		Transforms  ra and dec into galactic coordinates

		'''
		return self.coord.galactic
	
#------------------------------------------------------------------------------------------------------------------------------------------------------------
#Fluxes

	def xmm_get_flux(self):	
		'''
		Gets XMM Flux
		'''	
		try:
			#print(self.name, 'XMM flux')
			XMM_flux = heasarc.query_region(self.coord, mission='xmmssc', radius=Angle(10, 'arcsec'), fields = 'RA, DEC, EP_8_FLUX, EP_8_FLUX_ERROR, ERROR_RADIUS')
			#print(self.name, XMM_flux)
			XMM_min_flux = min(XMM_flux['EP_8_FLUX'])			
			XMM_min_Error = XMM_flux['EP_8_FLUX_ERROR'][np.argmin(XMM_flux['EP_8_FLUX'])]
			XMM_max_flux = max(XMM_flux['EP_8_FLUX'])
			XMM_max_Error = XMM_flux['EP_8_FLUX_ERROR'][np.argmax(XMM_flux['EP_8_FLUX'])]
			XMM_RA = XMM_flux['RA']
			XMM_DEC = XMM_flux['DEC']
			XMM_Error = [Angle(XMM_flux['ERROR_RADIUS'][x],'arcsec').degree for x in range(len(XMM_flux['ERROR_RADIUS']))]

			return XMM_min_flux, XMM_max_flux, XMM_RA, XMM_DEC, XMM_Error, XMM_min_Error, XMM_max_Error
		except (KeyError, TypeError, AttributeError, requests.exceptions.SSLError):
			return float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')


	def xrt_get_flux(self):
		'''
	 	Gets SWIFT XRT flux

		'''
		XRT = 'IX/58/2sxps' # https://arxiv.org/abs/1911.11710 (An improved and expanded Swift X-ray telescope point source catalog)
		XRT_cat = viz(columns = ['RAJ2000','DEJ2000','Err90','CR0','CFO0','CFU0','FPO0','E_FPO0','e_FPO0','FPU0','FpPO0','FpPU0','NH1H','E_NH1H','e_NH1H'])
		try:
			XRT_tab = XRT_cat.query_region(self.coord, catalog = XRT, radius=Angle(10,'arcsec'))
			XRT_min_flux = min(XRT_tab[0]['FPO0'])
			XRT_min_Error_pos = XRT_tab[0]['E_FPO0'][ma.argmin(XRT_tab[0]['FPO0'])]
			XRT_min_Error_neg = XRT_tab[0]['e_FPO0'][ma.argmin(XRT_tab[0]['FPO0'])]

			XRT_max_flux = max(XRT_tab[0]['FPO0'])
			XRT_max_Error_pos = XRT_tab[0]['E_FPO0'][ma.argmax(XRT_tab[0]['FPO0'])]
			XRT_max_Error_neg = XRT_tab[0]['e_FPO0'][ma.argmax(XRT_tab[0]['FPO0'])]


			XRT_min_peak_flux = min(XRT_tab[0]['FpPO0'])
			XRT_max_peak_flux = max(XRT_tab[0]['FpPO0'])

			if ma.is_masked(XRT_tab[0]['NH1H']) == True:
				XRT_nH_column_density = float('NaN')
				pos_nH_error = float('NaN')
				neg_nH_error = float('NaN')
			else:	
				XRT_nH_column_density = XRT_tab[0]['NH1H']
				pos_nH_error= XRT_tab[0]['E_NH1H']
				neg_nH_error= XRT_tab[0]['e_NH1H']

			XRT_Err90 = [Angle(XRT_tab[0]['Err90'][x],'arcsec').degree for x in range(len(XRT_tab[0]['Err90']))]
			XRT_RA = XRT_tab[0]['RAJ2000']
			XRT_DEC = XRT_tab[0]['DEJ2000']
			return XRT_min_flux, XRT_max_flux, XRT_min_peak_flux, XRT_max_peak_flux, XRT_nH_column_density, XRT_Err90, XRT_RA, XRT_DEC, XRT_min_Error_pos, XRT_min_Error_neg, XRT_max_Error_pos, XRT_max_Error_neg, pos_nH_error, neg_nH_error #5,6,7

		except IndexError:
			return float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')

	
	def bat_get_flux(self):
		'''
		Gets SWIFT BAT flux

		'''
		try:
			BAT_flux = heasarc.query_region(self.coord, mission='swbat105m', radius = Angle(1,'arcmin'), fields='RA, DEC, FLUX, FLUX_LOWER, FLUX_UPPER')
			BAT_min_flux = min(BAT_flux['FLUX'])
			BAT_max_flux = max(BAT_flux['FLUX'])
			BAT_RA = BAT_flux['RA']
			BAT_DEC = BAT_flux['DEC']
			return BAT_min_flux, BAT_max_flux, BAT_RA, BAT_DEC
		except (KeyError, TypeError,AttributeError, requests.exceptions.SSLError):
			return float('NaN'), float('NaN'), float('NaN'), float('NaN')


	def integral_get_flux(self):
		'''
		Gets INTEGRAl flux

		'''
		#print(self.name)
		INTEGRAL_Catalog = 'J/A+A/545/A27'
		try:
			INTEGRAL = viz(catalog = INTEGRAL_Catalog, column_filters = {'Name':self.name}, columns = ['RAJ2000', 'DEJ2000', 'Flux1', 'Flux2', 'Flux3', 'Name'])
			INTEGRAL_tab = INTEGRAL.query_region(self.coord, radius = Angle(10, 'deg'))
			INTEGRAL_min_flux = min(INTEGRAL_tab[0]['Flux1'])
			INTEGRAL_max_flux = max(INTEGRAL_tab[0]['Flux1'])
			INTEGRAL_min_flux = INTEGRAL_min_flux*(1.43*10**(-11))
			INTEGRAL_max_flux = INTEGRAL_max_flux*(1.43*10**(-11))
			Integral_RA = INTEGRAL_tab[0]['RAJ2000']
			Integral_DEC = INTEGRAL_tab[0]['DEJ2000']
			return INTEGRAL_min_flux, INTEGRAL_max_flux, Integral_RA, Integral_DEC
		except IndexError:
			try:
				if self.name in nam_s_arr: 			
					INTEGRAL = viz(catalog = INTEGRAL_Catalog, column_filters = {'Name':nam_s_arr_s[nam_s_arr.index(self.name)]}, columns = ['RAJ2000', 'DEJ2000', 'Flux1', 'Flux2', 'Flux3', 'Name'])
					INTEGRAL_tab = INTEGRAL.query_region(self.coord, radius = Angle(10, 'deg'))
					INTEGRAL_min_flux = min(INTEGRAL_tab[0]['Flux1'])
					INTEGRAL_max_flux = max(INTEGRAL_tab[0]['Flux1'])
					INTEGRAL_min_flux = INTEGRAL_min_flux*(1.43*10**(-11))
					INTEGRAL_max_flux = INTEGRAL_max_flux*(1.43*10**(-11))
					Integral_RA = INTEGRAL_tab[0]['RAJ2000']
					Integral_DEC = INTEGRAL_tab[0]['DEJ2000']
					return INTEGRAL_min_flux, INTEGRAL_max_flux, Integral_RA, Integral_DEC
				elif self.name in The_Same_L: 
					INTEGRAL = viz(catalog = INTEGRAL_Catalog, column_filters = {'Name':The_Same_R[The_Same_L.index(self.name)]}, columns = ['RAJ2000', 'DEJ2000', 'Flux1', 'Flux2', 'Flux3', 'Name'])
					INTEGRAL_tab = INTEGRAL.query_region(self.coord, radius = Angle(10, 'deg'))
					INTEGRAL_min_flux = min(INTEGRAL_tab[0]['Flux1'])
					INTEGRAL_max_flux = max(INTEGRAL_tab[0]['Flux1'])
					INTEGRAL_min_flux = INTEGRAL_min_flux*(1.43*10**(-11))
					INTEGRAL_max_flux = INTEGRAL_max_flux*(1.43*10**(-11))
					Integral_RA = INTEGRAL_tab[0]['RAJ2000']
					Integral_DEC = INTEGRAL_tab[0]['DEJ2000']
					return INTEGRAL_min_flux, INTEGRAL_max_flux, Integral_RA, Integral_DEC
				else:
					return float('NaN'), float('NaN'), float('NaN'), float('NaN')
			except IndexError:
				return float('NaN'), float('NaN'), float('NaN'), float('NaN')

	def chandra_get_flux(self):
		#Gets Chandra flux
		Chandra_v2 = 'IX/57/csc2master'
		Chandra = viz(catalog = Chandra_v2, columns = ['2CXO','RAICRS','DEICRS','r0','r1','Fluxw','Fluxb','b_Fluxb','B_Fluxb','b_Fluxw','B_Fluxw'])
		try:
			Chandra_flux=Chandra.query_region(self.coord, radius=Angle(2, 'arcsec'))
			if not ma.is_masked(Chandra_flux[0]['Fluxb']):

				Chandra_min_flux = min(Chandra_flux[0]['Fluxb'])
				Chandra_min_Error_pos = Chandra_flux[0]['B_Fluxb'][ma.argmin(Chandra_flux[0]['Fluxb'])]
				Chandra_min_Error_neg = Chandra_flux[0]['b_Fluxb'][ma.argmin(Chandra_flux[0]['Fluxb'])]
	
				Chandra_max_flux = max(Chandra_flux[0]['Fluxb'])
				Chandra_max_Error_pos = Chandra_flux[0]['B_Fluxb'][ma.argmax(Chandra_flux[0]['Fluxb'])]
				Chandra_max_Error_neg = Chandra_flux[0]['b_Fluxb'][ma.argmax(Chandra_flux[0]['Fluxb'])]
				Chandra_Instrument = 'ACIS'
			else:
				Chandra_min_flux = min(Chandra_flux[0]['Fluxw'])
				Chandra_min_Error_pos = Chandra_flux[0]['B_Fluxw'][ma.argmin(Chandra_flux[0]['Fluxw'])]
				Chandra_min_Error_neg = Chandra_flux[0]['b_Fluxw'][ma.argmin(Chandra_flux[0]['Fluxw'])]
	
				Chandra_max_flux = max(Chandra_flux[0]['Fluxw'])
				Chandra_max_Error_pos = Chandra_flux[0]['B_Fluxw'][ma.argmax(Chandra_flux[0]['Fluxw'])]
				Chandra_max_Error_neg = Chandra_flux[0]['b_Fluxw'][ma.argmax(Chandra_flux[0]['Fluxw'])]
				Chandra_Instrument = 'HRC'

			Chandra_RA=Chandra_flux[0]['RAICRS']
			Chandra_DEC=Chandra_flux[0]['DEICRS']

			r0 = Chandra_flux[0]['r0']
			r1 = Chandra_flux[0]['r1']
			r = [m.sqrt((r0[x])**2+(r1[x])**2) for x in range(len(r0))] 
			r_final = [Angle(r[x],'arcsec').degree for x in  range(len(r))]

			return Chandra_min_flux, Chandra_max_flux, Chandra_flux, Chandra_RA, Chandra_DEC, r_final, Chandra_min_Error_pos, Chandra_min_Error_neg, Chandra_max_Error_pos, Chandra_max_Error_neg, Chandra_Instrument
		except IndexError:
			return float('NaN'), float('NaN'), 'None', float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), 'None'
	

	def starhouse_get_info(self):
		Starhorse = 'I/349/starhorse'
		Starhorse_cat = viz(catalog = Starhorse, columns=['all'])
		if not m.isnan(self.GAIA_id):
			try:
				Starhorse_tab = Starhorse_cat.query_constraints(Source = self.GAIA_id)
				dist50 = Starhorse_tab[0]['dist50']*1000
				dist16 = Starhorse_tab[0]['dist16']*1000
				dist84 = Starhorse_tab[0]['dist84']*1000
				AV50 = Starhorse_tab[0]['AV50']
				AV16 = Starhorse_tab[0]['AV16']
				AV84 = Starhorse_tab[0]['AV84']

				AG50 = Starhorse_tab[0]['AG50']
				teff50 = Starhorse_tab[0]['teff50']
				teff16 = Starhorse_tab[0]['teff16']
				teff84 = Starhorse_tab[0]['teff84']
				logg50 = Starhorse_tab[0]['logg50']
				logg16 = Starhorse_tab[0]['logg50']
				logg84 = Starhorse_tab[0]['logg50']
				met50 = Starhorse_tab[0]['met50']
				met16 = Starhorse_tab[0]['met16']
				met84 = Starhorse_tab[0]['met84']
				mass50 = Starhorse_tab[0]['mass50']
				mass16 = Starhorse_tab[0]['mass16']
				mass84 = Starhorse_tab[0]['mass84']

				return dist50, dist16, dist84, AV50, AV16, AV84, teff50, teff16, teff84, logg50, logg16, logg84, met50, met16, met84, mass50, mass16, mass84, Starhorse_tab, AG50
			except IndexError:
				return float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), 'None', float('NaN')
		else:
			return float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), 'None', float('NaN')


	# def gaia_get_coord(self):
	# 	'''
	# 	Gets gaia coordinates from gaia id

	# 	'''
	# 	try:
	# 		gaia_tab = viz(catalog = 'I/350/gaiaedr3', columns = ['RA_ICRS', 'DE_ICRS']).query_constraints(Source = self.GAIA_id_)
	# 		gaia_coord = SkyCoord(ra = gaia_tab[0]['RA_ICRS'][0], dec = gaia_tab[0]['DE_ICRS'][0], unit = (u.deg, u.deg))
	# 	except IndexError:
	# 		try:
	# 			gaia_tab = viz(catalog = 'I/350/gaiaedr3', columns = ['RA_ICRS', 'DE_ICRS', 'Source']).query_object(self.name, radius = Angle(1,'arcsec'))
	# 			gaia_coord = SkyCoord(ra = gaia_tab[0]['RA_ICRS'][0], dec = gaia_tab[0]['DE_ICRS'][0], unit = (u.deg, u.deg))
	# 		except (TypeError, IndexError, KeyError):
	# 			if self.name in nam_s_arr:
	# 				gaia_tab = viz(catalog = 'I/350/gaiaedr3', columns = ['RA_ICRS', 'DE_ICRS', 'Source']).query_object(nam_s_arr_s[nam_s_arr.index(self.name)], radius = Angle(1,'arcsec'))
	# 				gaia_coord = SkyCoord(ra = gaia_tab[0]['RA_ICRS'][0], dec = gaia_tab[0]['DE_ICRS'][0], unit = (u.deg, u.deg))
	# 			elif self.name in The_Same_L:
	# 				gaia_tab = viz(catalog = 'I/350/gaiaedr3', columns = ['RA_ICRS', 'DE_ICRS', 'Source']).query_object(The_Same_R[The_Same_L.index(self.name)], radius = Angle(1,'arcsec'))
	# 				gaia_coord = SkyCoord(ra = gaia_tab[0]['RA_ICRS'][0], dec = gaia_tab[0]['DE_ICRS'][0], unit = (u.deg, u.deg))
	# 	return gaia_coord	
			


	# def gaia_get_err(self):
	# 	'''
	# 	Gets Position error for Gaia

	# 	'''
	# 	try:
	# 		gaia_tab = viz(catalog = 'I/350/gaiaedr3', columns = ['e_RA_ICRS','e_DE_ICRS']).query_constraints(Source = self.GAIA_id_)
	# 		err_1 = gaia_tab[0]['e_RA_ICRS'][0]
	# 		err_2 = gaia_tab[0]['e_DE_ICRS'][0]
	# 		sum_err = ( (err_1)**2 + (err_2)**2 )**(1/2)
	# 	except IndexError:
	# 		try:
	# 			gaia_tab = viz(catalog='I/350/gaiaedr3',columns=['e_RA_ICRS','e_DE_ICRS','Source']).query_object(self.name, radius=Angle(1,'arcsec'))
	# 			err_1 = gaia_tab[0]['e_RA_ICRS'][0]
	# 			err_2 = gaia_tab[0]['e_DE_ICRS'][0]
	# 			sum_err = ( (err_1)**2 + (err_2)**2 )**(1/2)
	# 		except (TypeError, IndexError, KeyError):
	# 			if self.name in nam_s_arr:
	# 				gaia_tab = viz(catalog='I/350/gaiaedr3',columns=['e_RA_ICRS','e_DE_ICRS','Source']).query_object(nam_s_arr_s[nam_s_arr.index(self.name)], radius=Angle(1,'arcsec'))
	# 				err_1 = gaia_tab[0]['e_RA_ICRS'][0]
	# 				err_2 = gaia_tab[0]['e_DE_ICRS'][0]
	# 				sum_err = ( (err_1)**2 + (err_2)**2 )**(1/2)
	# 			elif self.name in The_Same_L:
	# 				gaia_tab = viz(catalog='I/350/gaiaedr3',columns=['e_RA_ICRS','e_DE_ICRS','Source']).query_object(The_Same_R[The_Same_L.index(self.name)], radius=Angle(1,'arcsec'))
	# 				err_1 = gaia_tab[0]['e_RA_ICRS'][0]
	# 				err_2 = gaia_tab[0]['e_DE_ICRS'][0]
	# 				sum_err = ( (err_1)**2 + (err_2)**2 )**(1/2)
	# 	return sum_err


	def mass_get_flux(self):
		'''
	 	Gets 2MASS flux

		'''
		Mass = 'II/246/out'
		Mass_cat =  viz(catalog = 'II/246/out', columns=['2MASS', 'RAJ2000', 'DEJ2000', 'errMaj','errMin'])
		#Mass_cat =  viz(catalog = 'II/246/out', columns=['2MASS','errMaj','errMin'], column_filters = {'2MASS': self.MASS_id_})
		try:
			Mass_tab = Mass_cat.query_region(self.coord, radius = Angle(5, 'arcsec'))
			if len(Mass_tab [0]) == 1:
				Mass_RA = Mass_tab[0]['RAJ2000'][0]
				Mass_DEC = Mass_tab[0]['DEJ2000'][0]
				err_1 = Mass_tab[0]['errMaj'][0]
				err_2 = Mass_tab[0]['errMin'][0]
				sum_err = ( (err_1)**2 + (err_2)**2 )**(1/2)
				return Mass_RA, Mass_DEC, sum_err
			else:
				return float('NaN'), float('NaN'), float('NaN')
		except IndexError:
			return float('NaN'), float('NaN'), float('NaN')

	def gaia_get_flux(self):
		'''
	 	Gets Gaia flux

		'''
		gaia = 'I/350/gaiaedr3'
		gaia_cat = viz(catalog = gaia, columns = ['RA_ICRS', 'DE_ICRS', 'e_RA_ICRS', 'e_DE_ICRS', 'Source', 'Gmag'])
		try:
			gaia_tab = gaia_cat.query_region(self.coord, radius = Angle(5, 'arcsec'))
			if len(gaia_tab[0]) == 1:
				gaia_RA = gaia_tab[0]['RA_ICRS'][0]
				gaia_DEC = gaia_tab[0]['DE_ICRS'][0]
				err_1 = gaia_tab[0]['e_RA_ICRS'][0]
				err_2 = gaia_tab[0]['e_DE_ICRS'][0]
				sum_err = ( (err_1)**2 + (err_2)**2 )**(1/2)
				return gaia_RA, gaia_DEC, sum_err
			else:
				return float('NaN'), float('NaN'), float('NaN')
		except IndexError:
			return float('NaN'), float('NaN'), float('NaN')

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# Catalogs

	
	def catalog_get_type(self):
		'''
		Uses Rit and Kolb 2004 or Liu et al 2007 to get type of LMXB

		'''
		Rit_cat = 'B/cb/lmxbdata'
		Liu_cat = 'J/A+A/469/807/lmxb'
		
		Rit_Kolb = viz(columns = ['Type1', 'Type2', 'Type3', 'Type4' ])

		if self.name in The_Same_L: 
			Xray_type_R = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = The_Same_R[The_Same_L.index(self.name)])
		else: 
			Xray_type_R = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = self.name)

		Liu_etal = viz(columns = ['Type'])
		Xray_type_L = Liu_etal.query_constraints(catalog = Liu_cat, Name = self.name)
		
		#print(self.name, Xray_type_L[0],  Xray_type_R[0])
		
		if Xray_type_L != [] and Xray_type_R == []:
	#		 if self.name in The_Same_L:
	#			 print(The_Same_R[The_Same_L.index(self.name)])
	#		 elif self.name in Name_1:
	#			 print(self.name)
			if Xray_type_L[0]['Type'][0] != '':
				return Xray_type_L[0]['Type'][0]
			else:
				return 'None'

		elif Xray_type_L == [] and Xray_type_R != []:
			if Xray_type_R[0]['Type1'] == '' and  Xray_type_R[0]['Type2'] == '' and  Xray_type_R[0]['Type3'] == '' and  Xray_type_R[0]['Type4'] == '':
				return 'None'
			else:
				arr_list = []
				for i in range(1, 5):
					if Xray_type_R[0]['Type'+ str(i)] != '':
						arr_list.append(Xray_type_R[0]['Type'+ str(i)][0])
				types = arr_list[0]
				for j in range(1, len(arr_list)):
					types = types + ', ' + arr_list[j]
				return types

		elif Xray_type_L != [] and Xray_type_R != []:
			if Xray_type_R[0]['Type1'] == '' and  Xray_type_R[0]['Type2'] == '' and  Xray_type_R[0]['Type3'] == '' and  Xray_type_R[0]['Type4'] == '':
				if Xray_type_L[0]['Type'][0] != '':
					return Xray_type_L[0]['Type'][0]
				else:
					return 'None'
			else:
				arr_list = []
				for i in range(1, 5):
					if Xray_type_R[0]['Type'+ str(i)] != '':
						arr_list.append(Xray_type_R[0]['Type'+ str(i)][0])
				types = arr_list[0]
				for j in range(1, len(arr_list)):
					types = types + ', ' + arr_list[j]
				if Xray_type_L[0]['Type'][0] != '':
					return Xray_type_L[0]['Type'][0] + '|' + types
				else:
					return types
		else:
	#		 if self.name in The_Same_L:
	#			 print(The_Same_R[The_Same_L.index(self.name)])
	#		 elif self.name in Name_1:
	#			 print(self.name)
			return 'None'
	

	def catalog_alt_name(self):
		'''
		Uses Riiter and Kolb 2004 or Liu et al 2007 to get alternative name for source 

		'''
		Rit_Kolb = viz(columns = ['Name', 'AltName'])
		
		if self.name in The_Same_L: 
			Alt_Name_R = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = The_Same_R[The_Same_L.index(self.name)])
			Alt_Name_R_1 = Rit_Kolb.query_constraints(catalog = Rit_cat_note, Name = The_Same_R[The_Same_L.index(self.name)])
		else:
			Alt_Name_R = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = self.name)
			Alt_Name_R_1 = Rit_Kolb.query_constraints(catalog = Rit_cat_note, Name = self.name)

		Liu_etal = viz(columns = ['Name', 'Name2, Name3'])
		Alt_Name_L = Liu_etal.query_constraints(catalog = Liu_cat, Name = self.name)


		if Alt_Name_L != [] and Alt_Name_R != []:
			if Alt_Name_L[0]['Name2'][0] != '' and Alt_Name_R[0]['AltName'][0] != '':
				if Alt_Name_L[0]['Name2'][0] != Alt_Name_R[0]['AltName'][0] and self.name != Alt_Name_R[0]['AltName'][0]:
					return Alt_Name_L[0]['Name2'][0],  Alt_Name_R[0]['AltName'][0]
				else:
					if self.name in The_Same_L:
						if self.name != The_Same_R[The_Same_L.index(self.name)] and Alt_Name_L[0]['Name2'][0] != The_Same_R[The_Same_L.index(self.name)]:
							return Alt_Name_L[0]['Name2'][0], The_Same_R[The_Same_L.index(self.name)]
						else:
							if Alt_Name_R_1 != []:
								for j in range(len(Alt_Name_R_1[0])):
									if Alt_Name_R_1[0]['AltName'][j] != '' and Alt_Name_R_1[0]['AltName'][j] != Alt_Name_L[0]['Name2'][0] and self.name != Alt_Name_R_1[0]['AltName'][j]:
										alt_dop_name = Alt_Name_R_1[0]['AltName'][j]
										break
									else:
										continue
								try:
									return Alt_Name_L[0]['Name2'][0], alt_dop_name
								except NameError:
									if Alt_Name_L[0]['Name3'][0] == '':
										return Alt_Name_L[0]['Name2'][0], 'None'
									else:
										return Alt_Name_L[0]['Name2'][0], Alt_Name_L[0]['Name3'][0]
							else:
								if Alt_Name_L[0]['Name3'][0] == '':
									return Alt_Name_L[0]['Name2'][0], 'None'
								else:
									return Alt_Name_L[0]['Name2'][0], Alt_Name_L[0]['Name3'][0]
					else:
						print('smth is wrong, v1,', self.name)
						return 'None', 'None'

			elif Alt_Name_L[0]['Name2'][0] != '' and  Alt_Name_R[0]['AltName'][0] == '':
				if self.name in The_Same_L:
					if self.name != The_Same_R[The_Same_L.index(self.name)] and Alt_Name_L[0]['Name2'][0] != The_Same_R[The_Same_L.index(self.name)]:
						return Alt_Name_L[0]['Name2'][0], The_Same_R[The_Same_L.index(self.name)]
					else:
						if Alt_Name_R_1 != []:
							for j in range(len(Alt_Name_R_1[0])):
								if Alt_Name_R_1[0]['AltName'][j] != '' and Alt_Name_R_1[0]['AltName'][j] != Alt_Name_L[0]['Name2'][0] and self.name != Alt_Name_R_1[0]['AltName'][j]:
									alt_dop_name = Alt_Name_R_1[0]['AltName'][j]
									break
								else:
									continue
							try:
								return Alt_Name_L[0]['Name2'][0], alt_dop_name
							except NameError:
								if Alt_Name_L[0]['Name3'][0] == '':
									return Alt_Name_L[0]['Name2'][0], 'None'
								else:
									return Alt_Name_L[0]['Name2'][0], Alt_Name_L[0]['Name3'][0]
						else:
							if Alt_Name_L[0]['Name3'][0] == '':
								return Alt_Name_L[0]['Name2'][0], 'None'
							else:
								return Alt_Name_L[0]['Name2'][0], Alt_Name_L[0]['Name3'][0]
				else:
					print('smth is wrong, v2,', self.name)
					return 'None', 'None'


			elif Alt_Name_L[0]['Name2'][0] == '' and  Alt_Name_R[0]['AltName'][0] != '':
				if self.name in The_Same_L:
					if self.name != The_Same_R[The_Same_L.index(self.name)] and self.name != Alt_Name_R[0]['AltName'][0]:
						return Alt_Name_R[0]['AltName'][0], The_Same_R[The_Same_L.index(self.name)]
					elif self.name != The_Same_R[The_Same_L.index(self.name)] and self.name == Alt_Name_R[0]['AltName'][0]:
						if Alt_Name_R_1 != []:
							for j in range(len(Alt_Name_R_1[0])):
								if Alt_Name_R_1[0]['AltName'][j] != '' and Alt_Name_R_1[0]['AltName'][j] != The_Same_R[The_Same_L.index(self.name)] and self.name != Alt_Name_R_1[0]['AltName'][j]:
									alt_dop_name = Alt_Name_R_1[0]['AltName'][j]
									break
								else:
									continue
							try:
								return The_Same_R[The_Same_L.index(self.name)], alt_dop_name
							except NameError:
								if Alt_Name_L[0]['Name3'][0] != '' and Alt_Name_L[0]['Name3'][0] != The_Same_R[The_Same_L.index(self.name)]:
									return The_Same_R[The_Same_L.index(self.name)], Alt_Name_L[0]['Name3'][0]  
								else:
									return The_Same_R[The_Same_L.index(self.name)], 'None'							
						else:		
							if Alt_Name_L[0]['Name3'][0] != '' and Alt_Name_L[0]['Name3'][0] != The_Same_R[The_Same_L.index(self.name)]:
								return The_Same_R[The_Same_L.index(self.name)], Alt_Name_L[0]['Name3'][0]  
							else:
								return The_Same_R[The_Same_L.index(self.name)], 'None' 
					elif self.name == The_Same_R[The_Same_L.index(self.name)] and self.name != Alt_Name_R[0]['AltName'][0]:
						if Alt_Name_R_1 != []:
							for j in range(len(Alt_Name_R_1[0])):
								if Alt_Name_R_1[0]['AltName'][j] != '' and Alt_Name_R_1[0]['AltName'][j] != Alt_Name_R[0]['Name2'][0] and self.name != Alt_Name_R_1[0]['AltName'][j]:
									alt_dop_name = Alt_Name_R_1[0]['AltName'][j]
									break
								else:
									continue
							try:
								return Alt_Name_R[0]['Name2'][0], alt_dop_name
							except NameError:
								if Alt_Name_L[0]['Name3'][0] != '' and Alt_Name_L[0]['Name3'][0] != Alt_Name_R[0]['AltName'][0]:
									return Alt_Name_R[0]['AltName'][0], Alt_Name_L[0]['Name3'][0]  
								else:
									return Alt_Name_R[0]['AltName'][0], 'None'							
						else:		
							if Alt_Name_L[0]['Name3'][0] != '' and Alt_Name_L[0]['Name3'][0] != Alt_Name_R[0]['AltName'][0]:
								return Alt_Name_R[0]['AltName'][0], Alt_Name_L[0]['Name3'][0] 
							else:
								return Alt_Name_R[0]['AltName'][0], 'None'   
				else:
					print('smth is wrong, v3,', self.name)
					return 'None', 'None'


			else:
				if self.name in The_Same_L:
					if self.name != The_Same_R[The_Same_L.index(self.name)]:
						if Alt_Name_R_1 != []:
							for j in range(len(Alt_Name_R_1[0])):
								if Alt_Name_R_1[0]['AltName'][j] != '' and Alt_Name_R_1[0]['AltName'][j] != The_Same_R[The_Same_L.index(self.name)] and self.name != Alt_Name_R_1[0]['AltName'][j]:
									alt_dop_name = Alt_Name_R_1[0]['AltName'][j]
									break
								else:
									continue
							try:
								return The_Same_R[The_Same_L.index(self.name)], alt_dop_name
							except NameError:
								if Alt_Name_L[0]['Name3'][0] != '' and Alt_Name_L[0]['Name3'][0] != The_Same_R[The_Same_L.index(self.name)]:
									return The_Same_R[The_Same_L.index(self.name)], Alt_Name_L[0]['Name3'][0]  
								else:
									return The_Same_R[The_Same_L.index(self.name)], 'None'							
						else:		
							if Alt_Name_L[0]['Name3'][0] != '' and Alt_Name_L[0]['Name3'][0] != The_Same_R[The_Same_L.index(self.name)]:
								return The_Same_R[The_Same_L.index(self.name)], Alt_Name_L[0]['Name3'][0]  
							else:
								return The_Same_R[The_Same_L.index(self.name)], 'None'					 
				else:
					print('smth is wrong, v4,', self.name)
					return 'None', 'None'

		elif Alt_Name_L != [] and Alt_Name_R == []:
			if Alt_Name_L[0]['Name2'][0] == '':
				if Alt_Name_L[0]['Name3'][0] == '':
					return 'None', 'None'
				else:
					return 'None', Alt_Name_L[0]['Name3'][0]
			else:
				if Alt_Name_L[0]['Name3'][0] == '':
					return Alt_Name_L[0]['Name2'][0], 'None'
				else:
					return Alt_Name_L[0]['Name2'][0], Alt_Name_L[0]['Name3'][0]


		elif Alt_Name_L == [] and Alt_Name_R != []:
			if Alt_Name_R[0]['AltName'][0] == '':
				if Alt_Name_R_1 != []:
					for j in range(len(Alt_Name_R_1[0])):
						if Alt_Name_R_1[0]['AltName'][j] != '' and self.name != Alt_Name_R_1[0]['AltName'][j]:
							alt_dop_name = Alt_Name_R_1[0]['AltName'][j]
							break
						else:
							continue
					try:
						for j in range(len(Alt_Name_R_1[0])):
							if Alt_Name_R_1[0]['AltName'][j] != '' and self.name != Alt_Name_R_1[0]['AltName'][j] and Alt_Name_R_1[0]['AltName'][j] != alt_dop_name:
								alt_dop_name_1 = Alt_Name_R_1[0]['AltName'][j]
								break
							else:
								continue
						try:
							return alt_dop_name, alt_dop_name_1
						except NameError:
							return alt_dop_name, 'None'
					except NameError:
						return 'None', 'None'
				else:
					return 'None', 'None'
			else:
				if Alt_Name_R_1 != []:
					for j in range(len(Alt_Name_R_1[0])):
						if Alt_Name_R_1[0]['AltName'][j] != '' and self.name != Alt_Name_R_1[0]['AltName'][j] and Alt_Name_R_1[0]['AltName'][j] != Alt_Name_R[0]['AltName'][0]:
							alt_dop_name = Alt_Name_R_1[0]['AltName'][j]
							break
						else:
							continue
					try:
						return Alt_Name_R[0]['AltName'][0], alt_dop_name
					except NameError:
						return Alt_Name_R[0]['AltName'][0], 'None'
				else:
					return Alt_Name_R[0]['AltName'][0], 'None'

		else:
			return 'None', 'None'

	def catalog_opt_name(self):
		'''
		Uses Liu et al 2007 and Riiter and Kolb 2004 to get name (or spectral type) of optical counterpart 

		'''
		Rit_Kolb = viz(columns = ['SpType2'])

		if self.name in The_Same_L: 
			Opt_name_R = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = The_Same_R[The_Same_L.index(self.name)])
		else: 
			Opt_name_R = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = self.name)
			
		Liu_etal = viz(columns=['Opt'])
		Opt_name_L = Liu_etal.query_constraints(catalog = Liu_cat, Name = self.name)

		if Opt_name_L != [] and Opt_name_R == []:
			if Opt_name_L[0]['Opt'][0] == '':
				if not m.isnan(self.GAIA_id):
					return self.GAIA_id, 'None', 'Gaia'
				else:
					return 'None', 'None', 'None'
			else:
				return Opt_name_L[0]['Opt'][0], 'None', 'Liu'
		elif Opt_name_L == [] and Opt_name_R != []:
			if Opt_name_R[0]['SpType2'][0] == '':
				if not m.isnan(self.GAIA_id):
					return self.GAIA_id, 'None', 'Gaia'
				else:
					return 'None', 'None', 'None'
			else:
				if not m.isnan(self.GAIA_id):
					return self.GAIA_id, Opt_name_R[0]['SpType2'][0], 'Gaia'
				else:
					return 'None', Opt_name_R[0]['SpType2'][0], 'None'
		elif Opt_name_L == [] and Opt_name_R == []:
			if not m.isnan(self.GAIA_id):
				return self.GAIA_id, 'None', 'Gaia'
			else: 
				return 'None',  'None', 'None'
		else:
			#print('smth happend!')
			if Opt_name_L[0]['Opt'][0] == '' and Opt_name_R[0]['SpType2'][0] == '':
				if not m.isnan(self.GAIA_id):
					return self.GAIA_id, 'None', 'Gaia'
				else:
					return 'None', 'None', 'None'
			elif Opt_name_L[0]['Opt'][0] != '' and Opt_name_R[0]['SpType2'][0] == '':
				return Opt_name_L[0]['Opt'][0], 'None', 'Liu'
			elif Opt_name_L[0]['Opt'][0] == '' and Opt_name_R[0]['SpType2'][0] != '':
				if not m.isnan(self.GAIA_id):
					return self.GAIA_id, Opt_name_R[0]['SpType2'][0], 'Gaia'
				else:
					return 'None', Opt_name_R[0]['SpType2'][0], 'None'
			else:
				return Opt_name_L[0]['Opt'][0], Opt_name_R[0]['SpType2'][0], 'Liu'



	def catalog_porb_info(self):
		'''
		Uses Ritter and Kolb 2004 or Liu et al 2007 to get Porb and Ppulse info

		'''
		Rit_Kolb = viz(columns = ['Orb.Per', '3.__Per'])

		if self.name in The_Same_L: 
			Orb_Per_R = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = The_Same_R[The_Same_L.index(self.name)])
		else: 
			Orb_Per_R = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = self.name)

		Liu_etal = viz(columns = ['Porb', 'Ppulse'])
		Orb_Per_L = Liu_etal.query_constraints(catalog = Liu_cat, Name = self.name)

		if Orb_Per_L != [] and Orb_Per_R == []:
			return Orb_Per_L
		elif Orb_Per_L == [] and Orb_Per_R != []:
			return Orb_Per_R
		elif Orb_Per_L == [] and Orb_Per_R == []:
			return 0
		else:
			#print('smth happend!')
			return Orb_Per_L + Orb_Per_R


	def catalog_get_porb(self):
		'''
		Uses Ritter and Kolb 2004 or Liu et al 2007 to get Porb

		'''
		if self.Porb_info == 0:
			return float('Nan')
		else:
			try:
				if ma.is_masked(self.Porb_info[0]['Porb']) != True and ma.is_masked(self.Porb_info[1]['Orb.Per']) != True:
					lis_excep = ['GRS 1915+105 ', '1E 1740.7-2942', '3A 1728-247', '4U 1700+24']
					if self.name in lis_excep:
						return self.Porb_info[0]['Porb'][0]
						#return self.Porb_info[0]['Porb'][0], self.Porb_info[0]['Orb.Per'][0]
					else:
						return self.Porb_info[0]['Porb'][0] * 24
						#return self.Porb_info[0]['Porb'][0] * 24, self.Porb_info[0]['Orb.Per'][0]
				else:
					return float('Nan')
					#return float('Nan'), float('Nan')
			except (KeyError, IndexError):
				try:
					if ma.is_masked(self.Porb_info[0]['Porb']) != True:
						lis_excep = ['GRS 1915+105 ', '1E 1740.7-2942', '3A 1728-247', '4U 1700+24']
						if self.name  in lis_excep:
							return self.Porb_info[0]['Porb'][0]
						else:
							return self.Porb_info[0]['Porb'][0]*24
					else:
						return float('Nan')
				except (KeyError, IndexError):
					if ma.is_masked(self.Porb_info[0]['Orb.Per']) != True:
						return self.Porb_info[0]['Orb.Per'][0]
					else:
						return float('Nan')			


	def catalog_get_pulse(self):
		'''
		Uses Ritter and Kolb 2004 or Liu et al 2007 to get Ppulse

		'''
		if self.Porb_info == 0:
			return float('Nan')
		else:
			try:
				if ma.is_masked(self.Porb_info[0]['Ppulse']) != True and ma.is_masked(self.Porb_info[1]['_3.__Per']) != True:
					return self.Porb_info[0]['Ppulse'][0]
					#return self.Porb_info[0]['Ppulse'][0], self.Porb_info[0]['_3.__Per'][0]
				else:
					return float('Nan')
					#return float('Nan'), float('Nan')
			except (KeyError, IndexError):
				try:
					if ma.is_masked(self.Porb_info[0]['Ppulse']) != True:
						return self.Porb_info[0]['Ppulse'][0]
					else:
						return float('Nan')
				except (KeyError, IndexError):
					if ma.is_masked(self.Porb_info[0]['_3.__Per']) != True:
						return self.Porb_info[0]['_3.__Per'][0]
					else:
						return float('Nan')	

	def catalog_get_pos(self):
		'''
		Uses Ritter and Kolb 2004 or Liu et al 2007 to get type of observation to derive position (accuracy in this case) 

		'''
		Rit_Kolb = viz(columns=['epos'])
		
		if self.name in The_Same_L: 
			Pos_R = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = The_Same_R[The_Same_L.index(self.name)])
		else: 	
			Pos_R = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = self.name)

		Liu_etal = viz(columns=['Pos'])
		Pos_L = Liu_etal.query_constraints(catalog = Liu_cat, Name = self.name)

		if Pos_L != [] and Pos_R == []:
			return Pos_L[0]['Pos'][0]
		elif Pos_L == [] and Pos_R != []:
			return Pos_R[0]['epos'][0]
		elif Pos_L == [] and Pos_R == []:
			return 'None'
		else:
			#print('smth happend!')
			return Pos_L[0]['Pos'][0] + ';  ' + Pos_R[0]['epos'][0]




# ------------------------------------------------------------------------------------------------------------------------------------------------------------
# Add more

	def get_lx_lopt(self):
		'''
		Uses Ritter and Kolb 2004 or smth else to get Lx/Lopt

		'''
		Rit_Kolb = viz(columns=['all'])
		if self.name in The_Same_L: 
			Lx_Lopt = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = The_Same_R[The_Same_L.index(self.name)])
		else: 
			Lx_Lopt = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = self.name)

		if Lx_Lopt != []:
			if ma.is_masked(Lx_Lopt[0]['LX_Lopt'][0]) != True:
				return Lx_Lopt[0]['LX_Lopt'][0]
			else:
				return float('Nan')
		else:
		# Should be smth more!!!!
			return float('Nan')

	def get_mx_mopt(self):
		'''
		Uses Ritter and Kolb 2004 or smth else to get Mx/Mopt

		'''
		Rit_Kolb = viz(columns=['all'])
		if self.name in The_Same_L: 
			Mx_Mopt = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = The_Same_R[The_Same_L.index(self.name)])
		else: 
			Mx_Mopt = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = self.name)

		if Mx_Mopt != []:
			if ma.is_masked(Mx_Mopt[0]['M1_M2'][0]) != True:
				return Mx_Mopt[0]['M1_M2'][0]
			else:
				return float('Nan')
		else:
		# Should be smth more!!!!
			return float('Nan')
		
	def get_mx_mass(self):
		'''
		Uses Ritter and Kolb 2004 or smth else to get Mx mass

		'''
		Rit_Kolb = viz(columns=['all'])
		if self.name in The_Same_L: 
			Mx_Mopt = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = The_Same_R[The_Same_L.index(self.name)])
		else: 
			Mx_Mopt = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = self.name)

		if Mx_Mopt != []:
			if ma.is_masked(Mx_Mopt[0]['M1'][0]) != True:
				return Mx_Mopt[0]['M1'][0]
			else:
				return float('Nan')
		else:
		# Should be smth more!!!!
			return float('Nan')
		
		
	def get_mopt_mass(self):
		'''
		Uses Ritter and Kolb 2004 or smth else to get Mopt mass

		'''
		Rit_Kolb = viz(columns=['all'])
		if self.name in The_Same_L: 
			Mx_Mopt = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = The_Same_R[The_Same_L.index(self.name)])
		else: 
			Mx_Mopt = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = self.name)

		if Mx_Mopt != []:
			if ma.is_masked(Mx_Mopt[0]['M2'][0]) != True:
				return Mx_Mopt[0]['M2'][0]
			else:
				return float('Nan')
		else:
		# Should be smth more!!!!
			return float('Nan')



	def catalog_get_incl(self):
		'''
		Uses Ritter and Kolb 2004 to get inclination
		'''
		Rit_Kolb = viz(columns=['all'])
		if self.name in The_Same_L: 
			Lx_Lopt = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = The_Same_R[The_Same_L.index(self.name)])
		else: 
			Lx_Lopt = Rit_Kolb.query_constraints(catalog = Rit_cat, Name = self.name)

		if Lx_Lopt != []:
			if ma.is_masked(Lx_Lopt[0]['Incl'][0]) != True:
				return Lx_Lopt[0]['Incl'][0]
			else:
				return float('Nan')
		else:
		# Should be smth more!!!!
			return float('Nan')

	def get_b_v(self):
		'''
		Uses Ritter and Kolb 2004 to get inclination
		'''
		Liu_etal = viz(columns=['all'])
		B_V_tab = Liu_etal.query_constraints(catalog = Liu_cat, Name = self.name)

		if B_V_tab != []:
			if ma.is_masked(B_V_tab[0]['B-V'][0]) != True:
				return B_V_tab[0]['B-V'][0]
			else:
				return float('Nan')
		else:
		# Should be smth more!!!!
			return float('Nan')
		
	def get_eb_v(self):
		'''
		Uses Ritter and Kolb 2004 to get inclination
		'''
		Liu_etal = viz(columns=['all'])
		E_B_V_tab = Liu_etal.query_constraints(catalog = Liu_cat, Name = self.name)

		if E_B_V_tab != []:
			if ma.is_masked(E_B_V_tab[0]['E_B-V_'][0]) != True:
				return E_B_V_tab[0]['E_B-V_'][0]
			else:
				return float('Nan')
		else:
		# Should be smth more!!!!
			return float('Nan')

# ------------------------------------------------------------------------------------------------------------------------------------------------------------

	def get_crsf(self):
		'''
		Uses http://orma.iasfbo.inaf.it:7007/~mauro/pulsar_list.html to get CRSF

		'''
		if self.SIMBAD_ids is not None:
			for key in CRSF_dict:
				if key in self.SIMBAD_ids:
					return CRSF_dict[key]['CRSF']
				else: 
					continue
		else:
			try:
				CRF = CRSF_dict[self.name]['CRSF']
				return CRF
			except KeyError:
				try:
					if self.name in nam_s_arr:
						CRF = CRSF_dict[nam_s_arr_s[nam_s_arr.index(self.name)]]['CRSF']
						return CRF
					elif self.name in The_Same_L: 
						CRF = CRSF_dict[The_Same_R[The_Same_L.index(self.name)]]['CRSF']
						return CRF
				except KeyError:
					return 'None'

	def get_spin(self):
		'''
		Uses http://orma.iasfbo.inaf.it:7007/~mauro/pulsar_list.html to get Spin
		
		'''
		if self.SIMBAD_ids is not None:
			for key in CRSF_dict:
				if key in self.SIMBAD_ids:
					return CRSF_dict[key]['Spin']
				else: 
					continue
		else:
			try:
				Spin = CRSF_dict[self.name]['Spin']
				return Spin
			except KeyError:
				try:
					if self.name in nam_s_arr:
						Spin = CRSF_dict[nam_s_arr_s[nam_s_arr.index(self.name)]]['Spin']
						return Spin
					elif self.name in The_Same_L: 
						Spin = CRSF_dict[The_Same_R[The_Same_L.index(self.name)]]['Spin']
						return Spin
				except KeyError:
					return 'None'

	#-----------------------------------------------------------------------------------------------------

	def get_comments(self):
		return 'None'

	#-----------------------------------------------------------------------------------------------------


	#	return ChandraCatalog,Chandra_Flux,Chandra_pwl_Flux

	def cross_check_dr2_1(self):
		Gaia = 'I/345/gaia2'
		Gaia_tab = viz(catalog = Gaia, columns=['RA_ICRS','DE_ICRS','e_RA_ICRS','e_DE_ICRS','Source'])
		try:
			gaia_tab = Gaia_tab.query_object(self.name, radius = Angle(1, 'arcsec'))
			gaia_RA = gaia_tab[0]['RA_ICRS'] # gaia_tab[0]['RA_ICRS'][0]
			gaia_DE = gaia_tab[0]['DE_ICRS'] # gaia_tab[0]['DE_ICRS']]0]
			gaia_ID = gaia_tab[0]['Source']  # gaia_tab[0]['Source'][0]
			return gaia_RA, gaia_DE, gaia_ID 
		except IndexError:
			try:
				if self.name in nam_s_arr:
					gaia_tab = Gaia_tab.query_object(nam_s_arr_s[nam_s_arr.index(self.name)], radius = Angle(1, 'arcsec'))
					gaia_RA = gaia_tab[0]['RA_ICRS'] # gaia_tab[0]['RA_ICRS'][0]
					gaia_DE = gaia_tab[0]['DE_ICRS'] # gaia_tab[0]['DE_ICRS']]0]
					gaia_ID = gaia_tab[0]['Source']  # gaia_tab[0]['Source'][0]
					return gaia_RA, gaia_DE, gaia_ID 
				elif self.name in The_Same_L: 
					gaia_tab = Gaia_tab.query_object(The_Same_R[The_Same_L.index(self.name)], radius = Angle(1, 'arcsec'))
					gaia_RA = gaia_tab[0]['RA_ICRS'] # gaia_tab[0]['RA_ICRS'][0]
					gaia_DE = gaia_tab[0]['DE_ICRS'] # gaia_tab[0]['DE_ICRS']]0]
					gaia_ID = gaia_tab[0]['Source']  # gaia_tab[0]['Source'][0]
					return gaia_RA, gaia_DE, gaia_ID 
				else:
					return float('NaN'), float('NaN'), float('NaN')
			except IndexError:
				return float('NaN'), float('NaN'), float('NaN')


	def cross_check_dr2_2(self):
		Gaia = 'I/345/gaia2'
		Gaia_tab = viz(catalog = Gaia, columns=['RA_ICRS','DE_ICRS','e_RA_ICRS','e_DE_ICRS','Source'])
		try:
			gaia_tab = Gaia_tab.query_constraints(Source = self.GAIA_id)
			gaia_RA = gaia_tab[0]['RA_ICRS'] # gaia_tab[0]['RA_ICRS'][0]
			gaia_DE = gaia_tab[0]['DE_ICRS'] # gaia_tab[0]['DE_ICRS'][0]
			gaia_ID = gaia_tab[0]['Source']  # gaia_tab[0]['Source'][0]
			return gaia_RA, gaia_DE, gaia_ID 
		except IndexError:
			return float('NaN'), float('NaN'), float('NaN')


	def cross_check_dr2_3(self):
		Gaia = 'I/345/gaia2'
		Gaia_tab = viz(catalog = Gaia, columns=['RA_ICRS','DE_ICRS','e_RA_ICRS','e_DE_ICRS','Source'])
		try:
			gaia_tab = Gaia_tab.query_region(self.coord, radius = Angle(1, 'arcsec'))
			gaia_RA = gaia_tab[0]['RA_ICRS'] # gaia_tab[0]['RA_ICRS'][0]
			gaia_DE = gaia_tab[0]['DE_ICRS'] # gaia_tab[0]['DE_ICRS'][0]
			gaia_ID = gaia_tab[0]['Source']  # gaia_tab[0]['Source'][0]
			return gaia_RA, gaia_DE, gaia_ID 
		except IndexError:
			return float('NaN'), float('NaN'), float('NaN')


	def cross_check_dr3_1(self):
		Gaia3 ='I/350/gaiaedr3'
		Gaia_tab_ = viz(catalog = Gaia3, columns=['RA_ICRS','DE_ICRS','e_RA_ICRS','e_DE_ICRS','Source'])
		try:
			gaia_tab_ = Gaia_tab_.query_object(self.name, radius = Angle(1, 'arcsec'))
			gaia_RA = gaia_tab_[0]['RA_ICRS'] # gaia_tab_[0]['RA_ICRS'][0]
			gaia_DE = gaia_tab_[0]['DE_ICRS'] # gaia_tab_[0]['DE_ICRS'][0]
			gaia_ID = gaia_tab_[0]['Source']  # gaia_tab_[0]['Source'][0]
			return gaia_RA, gaia_DE, gaia_ID 
		except IndexError:
			try:
				if self.name in nam_s_arr:
					gaia_tab_ = Gaia_tab_.query_object(nam_s_arr_s[nam_s_arr.index(self.name)], radius = Angle(1, 'arcsec'))
					gaia_RA = gaia_tab_[0]['RA_ICRS'] # gaia_tab_[0]['RA_ICRS'][0]
					gaia_DE = gaia_tab_[0]['DE_ICRS'] # gaia_tab_[0]['DE_ICRS'][0]
					gaia_ID = gaia_tab_[0]['Source']  # gaia_tab_[0]['Source'][0]
					return gaia_RA, gaia_DE, gaia_ID 
				elif self.name in The_Same_L: 
					gaia_tab_ = Gaia_tab_.query_object(The_Same_R[The_Same_L.index(self.name)], radius = Angle(1, 'arcsec'))
					gaia_RA = gaia_tab_[0]['RA_ICRS'] # gaia_tab_[0]['RA_ICRS'][0]
					gaia_DE = gaia_tab_[0]['DE_ICRS'] # gaia_tab_[0]['DE_ICRS'][0]
					gaia_ID = gaia_tab_[0]['Source']  # gaia_tab_[0]['Source'][0]
					return gaia_RA, gaia_DE, gaia_ID
				else: 
					return float('NaN'), float('NaN'), float('NaN')
			except IndexError:
				return float('NaN'), float('NaN'), float('NaN')

	def cross_check_dr3_2(self):
		Gaia3 ='I/350/gaiaedr3'
		Gaia_tab_ = viz(catalog = Gaia3, columns=['RA_ICRS','DE_ICRS','e_RA_ICRS','e_DE_ICRS','Source'])
		try:
			gaia_tab_ = Gaia_tab_.query_constraints(Source = self.GAIA_id)
			gaia_RA = gaia_tab_[0]['RA_ICRS'] # gaia_tab_[0]['RA_ICRS'][0]
			gaia_DE = gaia_tab_[0]['DE_ICRS'] # gaia_tab_[0]['DE_ICRS'][0]
			gaia_ID = gaia_tab_[0]['Source']  # gaia_tab_[0]['Source'][0]
			return gaia_RA, gaia_DE, gaia_ID 
		except IndexError:
			return float('NaN'), float('NaN'), float('NaN')

	def cross_check_dr3_3(self):
		Gaia3 ='I/350/gaiaedr3'
		Gaia_tab_ = viz(catalog = Gaia3, columns=['RA_ICRS','DE_ICRS','e_RA_ICRS','e_DE_ICRS','Source'])
		try:
			gaia_tab_ = Gaia_tab_.query_region(self.coord, radius = Angle(1, 'arcsec'))
			gaia_RA = gaia_tab_[0]['RA_ICRS'] # gaia_tab_[0]['RA_ICRS'][0]
			gaia_DE = gaia_tab_[0]['DE_ICRS'] # gaia_tab_[0]['DE_ICRS'][0]
			gaia_ID = gaia_tab_[0]['Source']  # gaia_tab_[0]['Source'][0]
			return gaia_RA, gaia_DE, gaia_ID 
		except IndexError:
			return float('NaN'), float('NaN'), float('NaN')

	# def final_cross_check(self):

	# 	cross1 = self.cross_check_DR2_1[2] - self.cross_check_DR2_2[2]
	# 	cross2 = self.cross_check_DR2_1[2] - self.cross_check_DR2_3[2]
	# 	cross3 = self.cross_check_DR2_2[2] - self.cross_check_DR2_3[2]

	# 	cross4 = self.cross_check_DR2_1[2] - self.cross_check_DR3_1[2]
	# 	cross5 = self.cross_check_DR2_1[2] - self.cross_check_DR3_2[2]
	# 	cross6 = self.cross_check_DR2_1[2] - self.cross_check_DR3_3[2]

	# 	cross7 = self.cross_check_DR2_2[2] - self.cross_check_DR3_1[2]
	# 	cross8 = self.cross_check_DR2_2[2] - self.cross_check_DR3_2[2]
	# 	cross9 = self.cross_check_DR2_2[2] - self.cross_check_DR3_3[2]

	# 	cross10 = self.cross_check_DR2_3[2] - self.cross_check_DR3_1[2]
	# 	cross11 = self.cross_check_DR2_3[2] - self.cross_check_DR3_2[2]
	# 	cross12 = self.cross_check_DR2_3[2] - self.cross_check_DR3_3[2]

	# 	return cross1, cross2, cross3, cross4, cross5, cross6, cross7, cross8, cross9, cross10, cross11, cross12



#Writing all informations to VOTable

#Makes List of all sources in LMXB catalogs

#Apply LMXB calls on every source in Ritter and KOlb
All_sources = [LMXB(All_names[x]) for x in range(len(All_names)) ]  

votable = VOTableFile()

resource = Resource()

votable.resources.append(resource)

table = Table(votable)
resource.tables.append(table)

table.fields.extend([

		# Source Name 
		Field (votable, name = 'Name', datatype = 'unicodeChar', arraysize = '*', links = 'https://www.ivoa.net/documents/VOTable/20130920/REC-VOTable-1.3-20130920.html'),
		
		# Coordinates
		Field(votable, name = 'RA' , datatype = 'double', arraysize = '1'),
		Field(votable, name = 'DEC' , datatype = 'double', arraysize = '1'),

		# Coordinates
		Field(votable, name = 'Coord Error' , datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Coord Error Source' , datatype = 'unicodeChar', arraysize = '*'),

		Field(votable, name = 'G LON', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'G LAT', datatype = 'double', arraysize = '1'),

		Field (votable, name = 'X-ray Type', datatype = 'unicodeChar', arraysize= '*' ),
		Field (votable, name = 'PosType', datatype = 'unicodeChar', arraysize = '*' ),
		Field(votable, name = 'Porb', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Ppulse', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'AltName1', datatype = 'unicodeChar', arraysize = '*'),
		Field(votable, name = 'AltName2 ', datatype = 'unicodeChar', arraysize = '*'),
		Field(votable, name = 'OptName', datatype = 'unicodeChar', arraysize = '*'),
		Field(votable, name = 'OptRef', datatype = 'unicodeChar', arraysize = '*'),
		Field(votable, name = 'SpType', datatype = 'unicodeChar', arraysize = '*'),
		
		Field(votable, name = 'Spin', datatype = 'unicodeChar', arraysize = '*'),
		Field(votable, name = 'CRSF', datatype = 'unicodeChar', arraysize = '*'),
		Field(votable, name = 'Type', datatype = 'unicodeChar', arraysize = '*'),

		Field(votable, name = 'Gmag', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'GmagError', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Vmag', datatype = 'double', arraysize = '1'),

		Field(votable, name = 'Jmag', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Hmag', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Kmag', datatype = 'double', arraysize = '1'),

		Field(votable, name = 'JmagError', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'HmagError', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'KmagError', datatype = 'double', arraysize = '1'),

		Field(votable, name = 'XRT nH column density', datatype = 'double', arraysize = '*'),
		
		Field(votable, name = 'Flux_range', datatype = 'double', arraysize = '*'),

		Field(votable, name = 'Incl', datatype = 'double', arraysize = '*'),
		
		Field(votable, name = 'Lx_Lopt', datatype = 'double', arraysize = '*'),
		
		Field(votable, name = 'Mx_Mopt', datatype = 'double', arraysize = '*'),
		Field(votable, name = 'Mx', datatype = 'double', arraysize = '*'),
		Field(votable, name = 'Mopt', datatype = 'double', arraysize = '*'),

		Field(votable, name = 'IDS', datatype = 'unicodeChar', arraysize = '*'),


		Field (votable, name = 'Coordinates', datatype = 'unicodeChar', arraysize = '*'),
		Field(votable, name = 'COO_BIBCODE', datatype = 'unicodeChar', arraysize = '*'),

		# Distance
		Field(votable, name = 'Geometric Distance', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Geometric Distance lower', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Geometric Distance higher', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'PhotoGeometric Distance', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'PhotoGeometric Distance lower', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'PhotoGeometric Distance higher', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Gaia Parallax', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Gaia Parallax Error', datatype = 'double', arraysize = '1'),

		# Gaia ID
		Field(votable, name = 'Gaia DR2 Name', datatype = 'unicodeChar', arraysize = '*'),

		# Gaia Info 
		Field(votable, name = 'BP-RP', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'BP-G', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Teff', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Teff lower', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Teff upper', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'AG', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'AG lower', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'AG upper', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'E(BP-RP)', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'E(BP-RP) lower', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'E(BP-RP) upper', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Lum', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Lum lower', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Lum upper', datatype = 'double', arraysize = '1'),

	

		# XMM Flux
		Field(votable, name = 'XMM min Flux', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'XMM max Flux', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'XMM min Flux Error', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'XMM max Flux Error', datatype = 'double', arraysize = '1'),
		# Field(votable, name = 'XMM Flux Error', datatype = 'double', arraysize = '*'),


		
	
		# Integral Flux
		Field(votable, name = 'INTEGRAL min Flux', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'INTEGRAL max Flux', datatype = 'double', arraysize = '1'),

		# BAT Flux
		Field(votable, name = 'BAT min Flux', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'BAT max Flux', datatype = 'double', arraysize = '1'),
		# Field(votable, name = 'BAT Flux Lower', datatype = 'double', arraysize = '*'),
		# Field(votable, name = 'BAT Flux Upper', datatype = 'double', arraysize = '*'),


		# Chandra Flux
		Field(votable, name = 'Chandra min Flux', datatype = 'double', arraysize = '1'),
	 	Field(votable, name = 'Chandra min Flux Error pos', datatype = 'double', arraysize = '1'),
	 	Field(votable, name = 'Chandra min Flux Error neg', datatype = 'double', arraysize = '1'),

		Field(votable, name = 'Chandra max Flux', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Chandra max Flux Error pos', datatype = 'double', arraysize = '1'),
	 	Field(votable, name = 'Chandra max Flux Error neg', datatype = 'double', arraysize = '1'),
	 	Field(votable, name = 'Chandra Instrument', datatype = 'unicodeChar', arraysize = '*'),
	 	# Field(votable, name = 'Chandra powerlaw Flux', datatype = 'double', arraysize = '*'),
		
		# SWIFT XRT Flux
	 	Field(votable, name = 'XRT min Flux observed', datatype = 'double', arraysize = '1'),
	 	Field(votable, name = 'XRT min Flux Error pos', datatype = 'double', arraysize = '1'),
	 	Field(votable, name = 'XRT min Flux Error neg', datatype = 'double', arraysize = '1'),

	 	Field(votable, name = 'XRT max Flux observed', datatype = 'double', arraysize = '1'),
	 	Field(votable, name = 'XRT max Flux Error pos', datatype = 'double', arraysize = '1'),
	 	Field(votable, name = 'XRT max Flux Error neg', datatype = 'double', arraysize = '1'),

	 	# Field(votable, name = 'XRT Flux unabsorbed', datatype = 'unsignedByte', arraysize = '*'),

		Field(votable, name = 'XRT min Peak Flux observed', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'XRT max Peak Flux observed', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'positive Error XRT nH column density', datatype = 'double', arraysize = '*'),
		Field(votable, name = 'negative XRT nH column density', datatype = 'double', arraysize = '*'),

		Field(votable, name = 'Comments', datatype = 'unicodeChar', arraysize = '*'),

		Field(votable, name = 'Dist 50', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Dist 50 lower', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Dist 50 higher', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'AV 50', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'AV 50 lower', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'AV 50 higher', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'AG 50', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'teff 50', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'teff 50 lower', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'teff 50 higher', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'logg 50 ', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'logg 50 lower ', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'logg 50 higher', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'met 50', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'met 50 lower', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'met 50 higher', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'mass 50 ', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'mass 50 lower', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'mass 50 higher', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Starhorse', datatype = 'unicodeChar', arraysize = '*'),

		Field(votable, name = 'Chandra', datatype = 'unicodeChar', arraysize = '*'),

		Field(votable, name = 'XMM_RA',datatype = 'double', arraysize = '*'),
		Field(votable, name = 'XMM_DEC',datatype = 'double', arraysize = '*'),
		Field(votable, name = 'XMM_Error',datatype = 'double', arraysize = '*'),

		Field(votable, name = 'Integral_RA',datatype = 'double', arraysize = '*'),
		Field(votable, name = 'Integral_DEC',datatype = 'double', arraysize = '*'),

		Field(votable, name = 'BAT_RA', datatype = 'double', arraysize = '*'),
		Field(votable, name = 'BAT_DEC', datatype = 'double', arraysize = '*'),

		Field(votable, name = 'Chandra_RA', datatype = 'double', arraysize = '*'),
		Field(votable, name = 'Chandra_DEC', datatype = 'double', arraysize = '*'),
		Field(votable, name = 'Chandra_Error', datatype = 'double', arraysize = '*'),

		Field(votable, name = 'XRT_RA', datatype = 'double', arraysize = '*'),
		Field(votable, name = 'XRT_DEC', datatype = 'double', arraysize = '*'),
		Field(votable, name = 'XRT_Error', datatype = 'double', arraysize = '*'),

		Field(votable, name = 'SimbadRA', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'SimbadDEC', datatype = 'double', arraysize = '1'),

		Field(votable, name = '2MASS_RA', datatype = 'double', arraysize = '*'),
		Field(votable, name = '2MASS_DEC', datatype = 'double', arraysize = '*'),
		Field(votable, name = '2MASS_Error', datatype = 'double', arraysize = '*'),

		Field(votable, name = 'GAIA_RA', datatype = 'double', arraysize = '*'),
		Field(votable, name = 'GAIA_DEC', datatype = 'double', arraysize = '*'),
		Field(votable, name = 'GAIA_Error', datatype = 'double', arraysize = '*'),



		#Field(votable, name = 'XRT Peak Flux unabsorbed', datatype = 'bit', arraysize = '*')


		])


table.create_arrays(len(All_sources))

for x in range(len(All_sources)):
#print(All_sources[x].name, All_sources[x].Xtype, All_sources[x].Pos, All_sources[x].Porb, All_sources[x].Ppulse, All_sources[x].Opt_name, All_sources[x].alt_name, All_sources[x].Spin, All_sources[x].CRSF, All_sources[x].SIMBAD_info['OTYPES'][0])
	table.array[x]=(

		# Source Name
		All_sources[x].name,

		#Coordinates
		All_sources[x].ra,
		All_sources[x].dec,
		All_sources[x].coord_error,
		All_sources[x].err_orig,
		All_sources[x].lon,
		All_sources[x].lat,


		All_sources[x].Xtype,
		All_sources[x].Pos,
		All_sources[x].Porb,
		All_sources[x].Ppulse,
		All_sources[x].alt_name[0],
		All_sources[x].alt_name[1],
		All_sources[x].Opt_name[0],
		All_sources[x].Opt_name[2],
		All_sources[x].Opt_name[1],
		All_sources[x].Spin,
		All_sources[x].CRSF,
		All_sources[x].SIMBAD_info['OTYPES'][0],

		All_sources[x].Gmag1[0],
		All_sources[x].Gmag1[1],
		All_sources[x].SIMBAD_Vmag,

		0,
		0,
		0,
		0,
		0,
		0,

		All_sources[x].XRT_flux[4],
		0,
		All_sources[x].Incl,
		All_sources[x].Lx_Lopt,
		All_sources[x].Mx_Mopt,
		All_sources[x].Mx,
		All_sources[x].Mopt,

		All_sources[x].SIMBAD_ids,


		




		#Coordinates
		All_sources[x].coord,
		All_sources[x].coor_bib,
		#Distance
		All_sources[x].distance[0],
		All_sources[x].distance[1],
		All_sources[x].distance[2],
		All_sources[x].distance[3],
		All_sources[x].distance[4],
		All_sources[x].distance[5],
		All_sources[x].Plx,
		All_sources[x].e_Plx,

		# Gaia ID
		All_sources[x].GAIA_id,
		
		# Gaia Info
		All_sources[x].BP_RP,
		All_sources[x].BP_G,
		All_sources[x].Teff,
		All_sources[x].b_Teff,
		All_sources[x].B_Teff,
		All_sources[x].AG,
		All_sources[x].b_AG,
		All_sources[x].B_AG,
		All_sources[x].EBP_RP,
		All_sources[x].b_EBP_RP,
		All_sources[x].B_EBP_RP,
		All_sources[x].Lum,
		All_sources[x].b_Lum,
		All_sources[x].B_Lum,


		# XMM Flux
		All_sources[x].XMM_flux[0],
		All_sources[x].XMM_flux[1],
		All_sources[x].XMM_flux[5],
		All_sources[x].XMM_flux[6],
		# All_sources[x].XMM_flux['EP_8_FLUX_ERROR'],

		# Integral Flux
		All_sources[x].INTEGRAL_flux[0],
		All_sources[x].INTEGRAL_flux[1],
		
		# BAT Flux
		All_sources[x].BAT_flux[0],
		All_sources[x].BAT_flux[1], 
		
		# All_sources[x].BAT_flux[0]['FLUX_LOWER'],
		# All_sources[x].BAT_flux[0]['FLUX_UPPER'],
		
		# Chandra Flux
		All_sources[x].Chandra_min_flux,
		All_sources[x].Chandra[6],
		All_sources[x].Chandra[7],
		All_sources[x].Chandra_max_flux,
		All_sources[x].Chandra[8],
		All_sources[x].Chandra[9],
		All_sources[x].Chandra[10],
		# All_sources[x].CHANDRA_flux[1],
		
		# Chandra Flux Powerlaw
		# All_sources[x].CHANDRA_flux[2],
		

		# XRT Flux
		# Flux observed
		All_sources[x].XRT_flux[0],
		All_sources[x].XRT_flux[8],
		All_sources[x].XRT_flux[9],
		All_sources[x].XRT_flux[1],
		All_sources[x].XRT_flux[10],
		All_sources[x].XRT_flux[11],
		All_sources[x].XRT_flux[2],
		All_sources[x].XRT_flux[3],
		All_sources[x].XRT_flux[12],
		All_sources[x].XRT_flux[13],

		# Flux unabsorbed
		# All_sources[x].XRT_flux['FPU0'],
		All_sources[x].comments,
		All_sources[x].Starhorse[0],
		All_sources[x].Starhorse[1],
		All_sources[x].Starhorse[2],
		All_sources[x].Starhorse[3],
		All_sources[x].Starhorse[4],
		All_sources[x].Starhorse[5],
		All_sources[x].Starhorse[19],
		All_sources[x].Starhorse[6],
		All_sources[x].Starhorse[7],
		All_sources[x].Starhorse[8],
		All_sources[x].Starhorse[9],
		All_sources[x].Starhorse[10],
		All_sources[x].Starhorse[11],
		All_sources[x].Starhorse[12],
		All_sources[x].Starhorse[13],
		All_sources[x].Starhorse[14],
		All_sources[x].Starhorse[15],
		All_sources[x].Starhorse[16],
		All_sources[x].Starhorse[17],
		All_sources[x].Starhorse[18],
		All_sources[x].Chandra_flux,
		All_sources[x].XMM_flux[2],
		All_sources[x].XMM_flux[3],
		All_sources[x].XMM_flux[4],
		All_sources[x].INTEGRAL_flux[2],
		All_sources[x].INTEGRAL_flux[3],
		All_sources[x].BAT_flux[2],
		All_sources[x].BAT_flux[3],
		All_sources[x].Chandra[3],
		All_sources[x].Chandra[4],
		All_sources[x].Chandra[5],
		All_sources[x].XRT_flux[6],
		All_sources[x].XRT_flux[7],
		All_sources[x].XRT_flux[5],
		All_sources[x].Simbad_ra,
		All_sources[x].Simbad_dec,

		All_sources[x].MASS_flux[0],
		All_sources[x].MASS_flux[1],
		All_sources[x].MASS_flux[2],

		All_sources[x].GAIA_flux[0],
		All_sources[x].GAIA_flux[1],
		All_sources[x].GAIA_flux[2],

		# Peak Flux unabsorbed
		# All_sources[x].XRT_flux['FpPU0']
		)


votable.to_xml('Test_Source3.xml')


votable = VOTableFile()

resource = Resource()

votable.resources.append(resource)

table = Table(votable)
resource.tables.append(table)

table.fields.extend([

		# Source Name 
		Field (votable, name = 'Name', datatype = 'unicodeChar', arraysize = '*', links = 'https://www.ivoa.net/documents/VOTable/20130920/REC-VOTable-1.3-20130920.html'),
		
		# Coordinates
		Field(votable, name = 'RA' , datatype = 'double', arraysize = '1'),
		Field(votable, name = 'DEC' , datatype = 'double', arraysize = '1'),

		# Coordinates
		Field(votable, name = 'CoordError' , datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Origin' , datatype = 'unicodeChar', arraysize = '*'),

		Field(votable, name = 'GLON', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'GLAT', datatype = 'double', arraysize = '1'),

		Field (votable, name = 'XrayType', datatype = 'unicodeChar', arraysize= '*' ),
		Field (votable, name = 'PosType', datatype = 'unicodeChar', arraysize = '*' ),
		Field(votable, name = 'Porb', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Ppulse', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'AltName1', datatype = 'unicodeChar', arraysize = '*'),
		Field(votable, name = 'AltName2 ', datatype = 'unicodeChar', arraysize = '*'),
		Field(votable, name = 'OptName', datatype = 'unicodeChar', arraysize = '*'),
		Field(votable, name = 'OptRef', datatype = 'unicodeChar', arraysize = '*'),
		Field(votable, name = 'SpType', datatype = 'unicodeChar', arraysize = '*'),
		
		Field(votable, name = 'Spin', datatype = 'unicodeChar', arraysize = '*'),
		Field(votable, name = 'CRSF', datatype = 'unicodeChar', arraysize = '*'),
		Field(votable, name = 'Type', datatype = 'unicodeChar', arraysize = '*'),

		Field(votable, name = 'Gmag', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'GmagError', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Vmag', datatype = 'double', arraysize = '1'),

		Field(votable, name = 'Jmag', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Hmag', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Kmag', datatype = 'double', arraysize = '1'),

		Field(votable, name = 'JmagError', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'HmagError', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'KmagError', datatype = 'double', arraysize = '1'),

		Field(votable, name = 'XRTnH', datatype = 'double', arraysize = '1'),
		
		Field(votable, name = 'Flux_range', datatype = 'double', arraysize = '1'),

		Field(votable, name = 'Incl', datatype = 'double', arraysize = '1'),
		
		Field(votable, name = 'Lx_Lopt', datatype = 'double', arraysize = '1'),
		
		Field(votable, name = 'Mx_Mopt', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Mx', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'Mopt', datatype = 'double', arraysize = '1'),

		Field(votable, name = 'B_V', datatype = 'double', arraysize = '1'),
		Field(votable, name = 'E_B_V', datatype = 'double', arraysize = '1'),

		Field(votable, name = 'IDS', datatype = 'unicodeChar', arraysize = '1')
		])

table.create_arrays(len(All_sources))
for x in range(len(All_sources)):
	#print(All_sources[x].name)
	table.array[x]=(
		
		# Source Name
		All_sources[x].name,

		#Coordinates
		All_sources[x].ra,
		All_sources[x].dec,
		All_sources[x].coord_error,
		All_sources[x].err_orig,
		All_sources[x].lon,
		All_sources[x].lat,


		All_sources[x].Xtype,
		All_sources[x].Pos,
		All_sources[x].Porb,
		All_sources[x].Ppulse,
		All_sources[x].alt_name[0],
		All_sources[x].alt_name[1],
		All_sources[x].Opt_name[0],
		All_sources[x].Opt_name[2],
		All_sources[x].Opt_name[1],
		All_sources[x].Spin,
		All_sources[x].CRSF,
		All_sources[x].SIMBAD_info['OTYPES'][0],

		All_sources[x].Gmag1[0],
		All_sources[x].Gmag1[1],
		All_sources[x].SIMBAD_Vmag,

		0,
		0,
		0,
		0,
		0,
		0,

		All_sources[x].XRT_flux[4],
		0,
		All_sources[x].Incl,
		All_sources[x].Lx_Lopt,
		All_sources[x].Mx_Mopt,
		All_sources[x].Mx,
		All_sources[x].Mopt,

		All_sources[x].b_v,
		All_sources[x].eb_v,


		All_sources[x].SIMBAD_ids

		)

votable.to_xml('Cat_table1.xml')


votable = VOTableFile()

resource = Resource()

votable.resources.append(resource)

table = Table(votable)
resource.tables.append(table)

table.fields.extend([

	Field (votable, name = 'Name', datatype = 'unicodeChar', arraysize = '*'),

	Field (votable, name = 'Gaia DR2 Name RA', datatype = 'double', arraysize = '*'),
	Field (votable, name = 'Gaia DR2 Name DEC', datatype = 'double', arraysize = '*'),
	Field (votable, name = 'Gaia DR2 Name ID', datatype = 'double', arraysize = '*'),
	Field (votable, name = 'Gaia DR2 ID RA', datatype = 'double', arraysize = '*'),
	Field (votable, name = 'Gaia DR2 ID DEC', datatype = 'double', arraysize = '*'),
	Field (votable, name = 'Gaia DR2 ID ID', datatype = 'double', arraysize = '*'),
	Field (votable, name = 'Gaia DR2 Coord RA', datatype = 'double', arraysize = '*'),
	Field (votable, name = 'Gaia DR2 Coord DEC', datatype = 'double', arraysize = '*'),
	Field (votable, name = 'Gaia DR2 Coord ID', datatype = 'double', arraysize = '*'),

	Field (votable, name = 'Gaia DR3 Name RA', datatype = 'double', arraysize = '*'),
	Field (votable, name = 'Gaia DR3 Name DEC', datatype = 'double', arraysize = '*'),
	Field (votable, name = 'Gaia DR3 Name ID', datatype = 'double', arraysize = '*'),
	Field (votable, name = 'Gaia DR3 ID RA', datatype = 'double', arraysize = '*'),
	Field (votable, name = 'Gaia DR3 ID DEC', datatype = 'double', arraysize = '*'),
	Field (votable, name = 'Gaia DR3 ID ID', datatype = 'double', arraysize = '*'),
	Field (votable, name = 'Gaia DR3 Coord RA', datatype = 'double', arraysize = '*'),
	Field (votable, name = 'Gaia DR3 Coord DEC', datatype = 'double', arraysize = '*'),
	Field (votable, name = 'Gaia DR3 Coord ID', datatype = 'double', arraysize = '*'),


	# Field (votable, name = 'Cross1', datatype = 'double', arraysize = '*'),
	# Field (votable, name = 'Cross2', datatype = 'double', arraysize = '*'),
	# Field (votable, name = 'Cross3', datatype = 'double', arraysize = '*'),
	# Field (votable, name = 'Cross4', datatype = 'double', arraysize = '*'),
	# Field (votable, name = 'Cross5', datatype = 'double', arraysize = '*'),
	# Field (votable, name = 'Cross6', datatype = 'double', arraysize = '*'),
	# Field (votable, name = 'Cross7', datatype = 'double', arraysize = '*'),
	# Field (votable, name = 'Cross8', datatype = 'double', arraysize = '*'),
	# Field (votable, name = 'Cross9', datatype = 'double', arraysize = '*'),
	# Field (votable, name = 'Cross10', datatype = 'double', arraysize = '*'),
	# Field (votable, name = 'Cross11', datatype = 'double', arraysize = '*'),
	# Field (votable, name = 'Cross12', datatype = 'double', arraysize = '*'),

	])

table.create_arrays(len(All_sources))
for x in range(len(All_sources)):
	#print(All_sources[x].name, All_sources[x].cross_check_DR2_1)
	table.array[x]=(
		All_sources[x].name,
		All_sources[x].cross_check_DR2_1[0],
		All_sources[x].cross_check_DR2_1[1],
		All_sources[x].cross_check_DR2_1[2],
		All_sources[x].cross_check_DR2_2[0],
		All_sources[x].cross_check_DR2_2[1],
		All_sources[x].cross_check_DR2_2[2],
		All_sources[x].cross_check_DR2_3[0],
		All_sources[x].cross_check_DR2_3[1],
		All_sources[x].cross_check_DR2_3[2],

		All_sources[x].cross_check_DR3_1[0],
		All_sources[x].cross_check_DR3_1[1],
		All_sources[x].cross_check_DR3_1[2],
		All_sources[x].cross_check_DR3_2[0],
		All_sources[x].cross_check_DR3_2[1],
		All_sources[x].cross_check_DR3_2[2],
		All_sources[x].cross_check_DR3_3[0],
		All_sources[x].cross_check_DR3_3[1],
		All_sources[x].cross_check_DR3_3[2],

		# All_sources[x].cross_check[0],
		# All_sources[x].cross_check[1],
		# All_sources[x].cross_check[2],
		# All_sources[x].cross_check[3],
		# All_sources[x].cross_check[4],
		# All_sources[x].cross_check[5],
		# All_sources[x].cross_check[6],
		# All_sources[x].cross_check[7],
		# All_sources[x].cross_check[8],
		# All_sources[x].cross_check[9],
		# All_sources[x].cross_check[10],
		# All_sources[x].cross_check[11],


		)
votable.to_xml('CrossCheck.xml')