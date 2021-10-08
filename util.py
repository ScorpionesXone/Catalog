import astropy.units as u
import astropy.io.fits as pyfits
import astropy.coordinates
import glob
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from numpy import *
import astropy.wcs as wcs
import numpy
import pyvo as vo
from astroquery.skyview import SkyView
import aplpy, os
from astropy.coordinates import SkyCoord
from hips import HipsSurveyProperties
from hips import make_sky_image
from hips import WCSGeometry
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
import math as m
#import matplotlib
#matplotlib.use('pgf') 
#matplotlib.rcParams['pgf.preamble'] = [r'\usepackage{hyperref}', ]
 
lsun = 3.839e33
 
Simbad.add_votable_fields('biblio')
Simbad.add_votable_fields('otypes')
import matplotlib.pyplot as plt
# matplotlib.use('pgf')
# from matplotlib import rc, rcParams
# rc('text', usetex=True)
# rcParams['pgf.preamble'] = [r'\usepackage{hyperref}', ]

vvvjhips = 'http://alasky.u-strasbg.fr/VISTA/VVV_DR4/VISTA-VVV-DR4-J/'
sh = vo.dal.TAPService("https://gaia.aip.de/tap")
def drop_dupes(cat,tol=5*u.arcsec):
    """Returns only unique objects from a catalog (SkyCoord instance, here we assume that objecs within tol are the same object)"""
    while any(astropy.coordinates.match_coordinates_sky(cat,cat,nthneighbor=2)[1]<=tol):
        match = astropy.coordinates.match_coordinates_sky(cat,cat,nthneighbor=2)
        dupes = cat[match[0][match[1]<=tol]]
        pairs = [set(x) for x in transpose((astropy.coordinates.match_coordinates_sky(dupes,cat,nthneighbor=1)[0],astropy.coordinates.match_coordinates_sky(dupes,cat,nthneighbor=2)[0]))]
        upairs = []
        while len(pairs)>0:
            x = pairs.pop()
            upairs.append(x)
            try:
                i = where([x==k for k in pairs])[0][0]
                del(pairs[i])
            except:
                pass
        mask = ones(len(cat),dtype=int)
        mask[[list(x)[0] for x in upairs]]=0
        cat = cat[where(mask)]
    # cat = astropy.coordinates.SkyCoord(cat.ra.deg,cat.dec.deg,frame='fk5')
    return cat
    
def twomass2vvv(j,h,k):
    """Converts 2MASS JHK to VVV system using Eq from section 2.2. of https://www.aanda.org/articles/aa/pdf/2013/04/aa20046-12.pdf"""
    return j-0.077*(j-h), h+0.032*(j-h), k+0.01*(j-k)


def get_simbad_summary(pos, radius=25):
    """Get simbad summary for given pos"""
    sr = Simbad.query_region(pos,radius=radius*u.arcsec)
    try:
        refs = len(unique(concatenate([x.decode('utf-8').split('|') for x in sr['BIBLIO']])))
        otypes = unique(concatenate([x.decode('utf-8').split('|') for x in sr['OTYPES']]))
    except:
        refs = 0
        otypes = ['unknown']
    return refs, otypes


def resolve_name(name):
    """docstring for resolve_name"""
    query_res = astropy.coordinates.name_resolve.get_icrs_coordinates(name)
    return query_res.transform_to(astropy.coordinates.FK5)


def query_vvv(pos,radius):
    """That's to query VVVDR4 position and JHK via cone-search"""
    url = 'http://wfaudata.roe.ac.uk/vvvDR4-dsa/DirectCone?DSACAT=VVV_DR4&DSATAB=vvvSource& '
    rr = vo.dal.conesearch(url,pos,radius)
    # sid = rr['sourceID'].data
    j = rr['jAperMag3'].data
    h = rr['hAperMag3'].data
    k = rr['ksAperMag3'].data
    j[j<-10]=nan
    h[h<-10]=nan
    k[k<-10]=nan
    ra = rr['ra'].data
    dec = rr['dec'].data
    if len(ra)==0:
        raise BaseException
    return ra, dec, ones_like(ra)*u.arcsec, j,h,k


def query_allwise(pos,radius):
    """docstring for query_allwise"""
    rr = Vizier.query_region(pos,radius,catalog='II/328/allwise')[0]
    w1 = rr['W1mag'].data.data
    w2 = rr['W2mag'].data.data
    w3 = rr['W3mag'].data.data
    w4 = rr['W4mag'].data.data
    ra = rr['RAJ2000'].data.data
    dec = rr['DEJ2000'].data.data
    return ra, dec, w1,w2,w3,w4
    

def query_2mass(pos,radius):
    """That's to query 2MASS position and JHK magnitudes via cone-search"""
    # url = 'https://irsa.ipac.caltech.edu/SCS?table=fp_psc&format=votable'
    # rr = vo.dal.conesearch(url,pos,radius)
    try:
        rr = Vizier.query_region(pos,radius,catalog='II/246/out')[0]
        # sid = rr['sourceID'].data
        j = rr['Jmag'].data.data
        h = rr['Hmag'].data.data
        k = rr['Kmag'].data.data
        ra = rr['RAJ2000'].data.data
        dec = rr['DEJ2000'].data.data
    except:
        rr = Vizier.query_region(pos,radius,catalog='II/281/2mass6x')[0]
        # sid = rr['sourceID'].data
        j = rr['Jmag'].data.data
        h = rr['Hmag'].data.data
        k = rr['Kmag'].data.data
        ra = rr['RAJ2000'].data
        dec = rr['DEJ2000'].data
    j,h,k = twomass2vvv(j,h,k)
    if len(ra)==0:
        raise BaseException
    return ra, dec, 0.3*ones_like(ra)*u.arcsec, j,h,k

def query_nir(pos,radius):
    """Combine VVV and 2MASS by searching first VVV and then 2MASS"""
    try:
        return query_vvv(pos,radius)
    except:
        try:
            return query_2mass(pos,radius)
        except:
            raise BaseException


def query_starhorse(pos,radius):
    """Query starhorse DB. Seems to be unavailable for most sources"""
    query = "select top 10 g.source_id, g.ra, g.dec,g.phot_g_mean_mag,s.dist50,s.av50,s.teff50,s.ruwe,s.mg0,s.ag50 from gdr2.gaia_source AS g, gdr2_contrib.starhorse as s where g.source_id = s.source_id and 1=CONTAINS(POINT('ICRS',g.ra,g.dec), CIRCLE(’ICRS’,%f,%f,%f))"%(pos.ra.deg,pos.dec.deg,radius.to(u.deg).value)
    #d.r_est from gdr2_contrib.geometric_distance as d where d.source_id = s.source_id and
    resultset = sh.search(query)
    # return resultset
    sid = resultset['source_id'][0]
    mg = resultset['mg0'][0]
    d50 = resultset['dist50'][0]
    ag = resultset['ag50'][0]
    lum = 3.013e35*10**(-0.4*mg)
    return sid, mg, ag, d50, lum
    # return resultset


def query_gaia(pos,radius):
    """That's to query Gaia DR2"""
    rr = Gaia.cone_search(pos,radius).get_results()
    # sid = rr['sourceID'].data
    sid = rr['source_id'].data.data
    parallax = rr['parallax'].data.data
    parallax_over_error = rr['parallax_over_error'].data.data
    # print(parallax_over_error)
    parallax[parallax_over_error<=1]=nan
        # paralax=nan
    ra = rr['ra'].data.data
    dec = rr['dec'].data.data
    en = rr['astrometric_excess_noise'].data.data
    g = rr['phot_g_mean_mag'].data.data
    bprp = rr['bp_rp'].data.data
    lum = rr['lum_val'].data.data
    vari = rr['phot_variable_flag'].data.data
    ag = rr['a_g_val'].data.data
    #Now lets update distance from BJ paper
    shflag = 0
    try:        
        bjres = Vizier.query_region(pos,radius,catalog='I/347/gaia2dis')
        bjid = bjres[0]['Source'].data.data
        bjrest = bjres[0]['rest'].data.data
        for x in zip(bjid,bjrest):
            d[sid==x[0]]=x[1]
        print("Updated distance from BJ")
    except:
        pass

    #To have it in kpc
    try:
        shres = query_starhorse(pos,radius)
        if isnan(lum) and not(isnan(shres[-1])):
            lum = shres[-1]/slum
            shflag=1
        if isnan(d) and not(isnan(shres[3])):
            d = shres[3]
            shflag=1
        if isnan(ag) and not(isnan(shres[2])):
            ag = shres[2]
            shflag=1
        print("Updated distance from StarHorse")
    except:
        d = 1000./parallax
    d/=1000
    return ra, dec, g, bprp, d, sid, en, lum, parallax_over_error, ag, shflag
    #external.external.gaiadr2_geometric_distance, gaiadr2.gaiadr2.ruwe

def query_skymapper(pos,radius):
    """Query SkyMapper"""
    url = "http://skymapper.anu.edu.au/sm-cone/public/query?"
    rr = vo.dal.conesearch(url,pos,radius)
    fields = ['raj2000','dej2000','u_psf','v_psf','g_psf','r_psf','i_psf','z_psf']
    ra,de,u,v,g,r,i,z = [rr[x].data for x in fields]
    return ra, de, g,r,i,z





def query_csc2(pos,radius):
    """That's to query SCS v2 via cone-search"""
    url = 'http://cda.cfa.harvard.edu/cscvo/coneSearch'
    rr = vo.dal.conesearch(url,pos,radius)
    sid = rr['name'].data
    j = rr['flux_aper_b'].data
    ra = rr['ra'].data
    dec = rr['dec'].data
    err = rr['err_ellipse_r0'].data
    err = sqrt(err**2+1.2**2)
    pp = astropy.coordinates.SkyCoord(ra,dec,unit=u.deg,frame='fk5')
    # ra = pp[argmin(pos.separation(pp))].ra.deg
    # dec = pp[argmin(pos.separation(pp))].dec.deg
    # err = err[argmin(pos.separation(pp))]
    ra = pp[argmax(j)].ra.deg
    dec = pp[argmax(j)].dec.deg
    err = err[argmax(j)]    
    return ra, dec, sqrt(err**2+2**2)*u.arcsec

def query_ssc(pos,radius):
    """That's to query SSC DR9 via cone-search"""
    v = Vizier(columns=['RA_ICRS','DE_ICRS','ePos'])
    rr = v.query_region(pos,radius,catalog='IX/59/xmm4dr9s')
    # j = rr['ep_8_flux'].data
    ra = rr[0]['RA_ICRS'].data.data
    dec = rr[0]['DE_ICRS'].data.data
    pe = rr[0]['ePos'].data.data
    ra = median(ra)
    dec = median(dec)
    return ra,dec, sqrt(pe**2+2**2)[0]*u.arcsec

def query_sxps(pos,radius):
    """That's to query VVVDR4 via cone-search"""
    # v = Vizier(columns=['RAJ2000','DEJ2000','NH1H','E_NH1H','e_NH1H','NH2A','E_NH2A','e_NH2A','NH2C','E_NH2C','e_NH2C'])
    rr = Vizier.query_region(pos,radius,catalog='IX/58/2sxps')
    sid = rr[0]['_2SXPS'].data.data
    # j = rr['ep_8_flux'].data
    ra = rr[0]['RAJ2000'].data.data
    dec = rr[0]['DEJ2000'].data.data
   # print(ra)
    #print(dec)
    return ra, dec,25*u.arcsec

def query_2rsx(pos,radius):
    """That's to query VVVDR4 via cone-search"""
    rr = Vizier.query_region(pos,radius,catalog='J/A+A/588/A103/cat2rxs')
    sid = rr[0]['_2RXS'].data.data
    # j = rr['ep_8_flux'].data
    ra = rr[0]['RAJ2000'].data.data
    dec = rr[0]['DEJ2000'].data.data
    print(ra)
    print(dec)
    return ra, dec, 17*u.arcsec

gar = vo.dal.TAPService("https://gea.esac.esa.int/tap-server/tap")

def mk_finding_chart(SimbadRA,SimbadDEC,pos,poserr,error_name,eroid,XMM_RA,XMM_DEC,XMM_Error,XRT_Error,XRT_RA,XRT_DEC,BAT_RA,BAT_DEC,INTEGRAL_RA,INTEGRAL_DEC,CHANDRA_RA,CHANDRA_DEC,CHANDRA_Error,twoMass_RA,twoMass_DEC,twoMass_Error,GAIA_RA,GAIA_DEC,GAIA_Error,cut_size=1):
    print(eroid)
   # print(CHANDRA_RA)
   # print(XMM_RA)
   # print(type(XMM_RA))

    """Make an 2MASS J image and overplot with region corresponding to eposerr"""
    # clf()
    surl = 'http://simbad.u-strasbg.fr/simbad/sim-coo?Coord=RACOORD%2C+DECOORD&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=25&Radius.unit=arcsec&submit=submit+query'
    surl = surl.replace("RACOORD",str(pos.ra.deg)).replace('DECOORD',str(pos.dec.deg))
    if type(pos)==str:
        pos = resolve_name(pos)
    #sep = pos.separation(epos)
    #ii = argmin(sep)
    # rr = sqrt(e_poserr[ii]**2+ero_syserr**2)
    rr = poserr
    rrr=Angle(poserr,"arcsec").degree

    #Let's try to get VVV J image first. 
    try:
        geometry = WCSGeometry.create(skydir=SkyCoord(pos.ra.deg, pos.dec.deg, unit='deg', frame='icrs'),width=512, height=512, fov=cut_size*u.arcmin, coordsys='icrs', projection='TAN')
        result = make_sky_image(geometry, vvvjhips,'fits',precise=True,progress_bar=False)
        pyfits.ImageHDU(data = result.image, header=result.geometry.wcs.to_header()).writeto('tmp.fits',overwrite=True)
        ff = aplpy.FITSFigure('tmp.fits')
    except:
        paths = SkyView.get_images(position=pos,survey=['2MASS-J'],width=cut_size*u.arcmin,height=cut_size*u.arcmin)
        print(pos)
        ff = aplpy.FITSFigure(paths[0])
    # pp = epos[ii]
    pp = pos
    if eroid=="2S 0053+604":
        ff.show_grayscale(vmin=0.1,vmax=1.160e+04, invert=False,stretch='log')
    else:

        ff.show_grayscale(invert=True,stretch='log')
    # ff.show_colorscale(stretch='log',cmap=cm.seismic)
    #open('tt.reg','w').write('global color=green\nfk5\ncircle(%f,%f,%f")\ncircle(%f,%f,%f")\npoint(%f,%f) #point=cross'%(pp.ra.deg,pp.dec.deg,rr,pp.ra.deg,pp.dec.deg,e_poserr[ii],pos.ra.deg,pos.dec.deg))
    #ff.show_regions('tt.reg')
    #Now lets get gaia and 2mass results
    query = "select * from gaiadr2.gaia_source as g where 1=CONTAINS(POINT('ICRS',g.ra,g.dec), CIRCLE('ICRS',%f,%f,%f))"%(pos.ra.deg,pos.dec.deg,cut_size/60.)
    resultset = gar.search(query)
    ra, dec = resultset['ra'],resultset['dec']
    x,y = ff.world2pixel(ra,dec)
    plt.scatter(x,y,log10(resultset['phot_g_mean_flux']),c='g',label='Gaia DR2', alpha=0.5)
    query = "select * from gaiadr1.tmass_original_valid as g where 1=CONTAINS(POINT('ICRS',g.ra,g.dec), CIRCLE('ICRS',%f,%f,%f))"%(pos.ra.deg,pos.dec.deg,cut_size/60.)
    resultset = gar.search(query)
    ra, dec = resultset['ra'],resultset['dec']
    x,y = ff.world2pixel(ra,dec)
    plt.plot(x,y,'rx',ms=5,alpha=0.5,label='2MASS')
    x,y=pos.ra.deg,pos.dec.deg
    ff.show_markers(x,y,coords_frame="world",c="orange",marker="*",label="Pos: {}".format(error_name))
    ff.show_circles(x,y,rrr,coords_frame="world",alpha=0.5,edgecolor="orange",linewidth=1)    
    ff.show_markers(SimbadRA,SimbadDEC,coords_frame="world",c="red",marker="+",alpha=0.5,s=100,label="Simbad")
 
    #ra, dec = query_sxps(pos,cut_size*u.arcmin)
    #x,y = ff.world2pixel(ra,dec)
    #plt.plot(x,y,'g*',ms=5,alpha=0.5,label='2SXPS')
    try:
        ra, dec,arc= query_2rsx(pos,cut_size*u.arcmin)
        x,y = ff.world2pixel(ra,dec)
        print("VVV")
        plt.plot(x,y,'r+',ms=5,alpha=0.5,label='VVV DR4')
    except:
        pass
    #try:
        #ra, dec,arc = query_csc2(pos,cut_size*u.arcmin)
        #x,y = ff.world2pixel(ra,dec)
        #print("CSC")
        #plt.plot(x,y,'go',ms=5,alpha=0.5,label='CSC2')
    if m.isnan(CHANDRA_RA[0])==False:
        #print(m.isnan(CHANDRA_RA))
        ff.show_markers(CHANDRA_RA,CHANDRA_DEC,coords_frame="world",c="gold",marker="3",alpha=0.5,label="Chandra")
        ff.show_circles(CHANDRA_RA,CHANDRA_DEC,CHANDRA_Error,coords_frame="world",alpha=0.5,edgecolor="gold",linewidth=1)
        print("Chandra")
    else:
        pass
    #except:
    #    pass
    if m.isnan(XMM_RA[0])==False:
    #ra, dec,arc = query_ssc(pos,cut_size*u.arcmin)
        x,y = ff.world2pixel(XMM_RA,XMM_DEC)
        print("SSC")
        ff.show_markers(XMM_RA,XMM_DEC,coords_frame="world",c="blue",marker="d",alpha=0.5,label="XMM")
        ff.show_circles(XMM_RA,XMM_DEC,XMM_Error,coords_frame="world",alpha=0.5,edgecolor="blue",linewidth=1)
    #print(ff.list_layer())    #x,y = ff.world2pixel(ra,dec)    print("SSC")
    #plt.plot(x,y,'gd',ms=5,alpha=0.5,label='SSC')  
    else:
        pass
    if m.isnan(twoMass_RA[0])==False:
    #ra, dec,arc = query_ssc(pos,cut_size*u.arcmin)
        x,y = ff.world2pixel(twoMass_RA,twoMass_DEC)
        print("2MASS")
        ff.show_markers(twoMass_RA,twoMass_DEC,coords_frame="world",c="cyan",marker="2",alpha=0.5,label="2Mass")
        ff.show_circles(twoMass_RA,twoMass_DEC,twoMass_Error,coords_frame="world",alpha=0.5,edgecolor="cyan",linewidth=1)
    #print(ff.list_layer())    #x,y = ff.world2pixel(ra,dec)    print("SSC")
    #plt.plot(x,y,'gd',ms=5,alpha=0.5,label='SSC')  
    else:
        pass
    if m.isnan(GAIA_RA[0])==False:
    #ra, dec,arc = query_ssc(pos,cut_size*u.arcmin)
        x,y = ff.world2pixel(GAIA_RA,GAIA_DEC)
        print("GAIA")
        ff.show_markers(GAIA_RA,GAIA_DEC,coords_frame="world",c="hotpink",marker="3",alpha=0.5,label="GAIA")
        ff.show_circles(GAIA_RA,GAIA_DEC,GAIA_Error,coords_frame="world",alpha=0.5,edgecolor="hotpink",linewidth=1)
    #print(ff.list_layer())    #x,y = ff.world2pixel(ra,dec)    print("SSC")
    #plt.plot(x,y,'gd',ms=5,alpha=0.5,label='SSC')  
    else:
        pass

    if m.isnan(XRT_RA[0])==False:
        #ra, dec,arc = query_sxps(pos,cut_size*u.arcmin)
        x,y = ff.world2pixel(XRT_RA,XRT_DEC) 
        ff.show_markers(XRT_RA,XRT_DEC,coords_frame="world",c="green",marker="*",alpha=0.5,label="Swift XRT")
        ff.show_circles(XRT_RA,XRT_DEC,XRT_Error,coords_frame="world",alpha=0.5,edgecolor="green",linewidth=1)
        print("2SXPS")

        #plt.plot(x,y,'g*',ms=5,alpha=0.5,label='2SXPS')       
    else:
        pass
    if m.isnan(BAT_RA[0])==False:
        ff.show_markers(BAT_RA,BAT_DEC,coords_frame="world",c="yellow",marker="p",alpha=0.5,label="Swift BAT")
        ff.show_circles(BAT_RA,BAT_DEC,0.23,coords_frame="world",alpha=0.5,edgecolor="yellow",linewidth=1) #Bat uncertainty  13.8 arcmin
        print("BAT")
        #ff.show_circles(CHANDRA_RA,CHANDRA_DEC,CHANDRA_Error,coords_frame="world",alpha=1,edgecolor="purple",linewidth=1)       
    else:
        pass
    if m.isnan(INTEGRAL_RA[0])==False:
        ff.show_markers(INTEGRAL_RA,INTEGRAL_DEC,coords_frame="world",c="lime",marker="v",alpha=0.5,label="Integral")
        ff.show_circles(INTEGRAL_RA,INTEGRAL_DEC,0.05,coords_frame="world",alpha=0.5,edgecolor="yellow",linewidth=1)#Integral uncertainity 3 arcmin
        print("Integral")
    else:
        pass


    try:
        refs, types = get_simbad_summary(pos)
    except:
        refs, types = 0,''
    #eroid = closest_erosource(pos)[-1]
    #plt.title(r'eID:\href{%s}{%s}, refs=%d, types=%s'%(surl,eroid,refs,",".join(types)))
    #plt.title(r"\href {}{}".format(surl,eroid))
    #ax.set_label()
    ax=plt.gca()
    print(ax)
    plt.legend(frameon=False)
    #ax.get_legend_handler_map()


    #ax.get_legend_handler_map()
    #print(test.get_legend_handler_map())
    #plt.get_legend_handler_map()
    ff.save('charts2/%s.pdf'%eroid)
    plt.close()
    #plt.savefig('charts/%s.pdf'%eroid)
    #plt.show()
    