from sherpa.astro.ui import *
from astropy.io import fits
from astropy.table import Table
from astropy.table import Column
from IPython.core.display import Image
from gammapy.image import SkyMapCollection, SkyMap
from gammapy.utils.energy import EnergyBounds
import pylab as pt
from gammapy.background import fill_acceptance_image
from astropy.coordinates import Angle
from astropy.units import Quantity
import numpy as np
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.coordinates import SkyCoord
from sherpa.models import Gauss2D
from sherpa_models import normgauss2dint
from cat_ana2 import *

"""
energy_bins = EnergyBounds.equal_log_spacing(0.5, 100, 5, 'TeV')
E1=energy_bins[0].value
E2=energy_bins[1].value
excess = SkyMapCollection.read("fov_bg_maps"+str(E1)+"_"+str(E2)+"_TeV.fits")["excess"]
load_psf("psf_SgrA", "psf_image_SgrA_" + str(E1) + "_" + str(E2) + ".fits")
data = fits.open("residual_0.5_1.44269990591_TeV.fits")
load_image(1, data)
large_gaus = Gauss2D("g2")
source_center_SgrA = SkyCoord.from_name("SgrA*")
large_gaus.xpos,large_gaus.ypos=skycoord_to_pixel(source_center_SgrA, excess.wcs)
CS_map=SkyMap.read("CStot.fits")
cs_reproj=CS_map.reproject(excess)
cs_reproj.data[np.where(np.isnan(cs_reproj.data))]=0
#cs_reproj.data[np.where(cs_reproj.data<50)]=0
cs_reproj.write("cs_map_reproj.fits", clobber=True)
load_table_model("CS","cs_map_reproj.fits")
set_full_model(psf_SgrA(large_gaus*CS))
#large_gaus.fwhm=150
#freeze(large_gaus.fwhm)
fit()
"""
pt.ion()
energy_bins = EnergyBounds.equal_log_spacing(0.5, 100, 5, 'TeV')

#E1=energy_bins[1].value
#E2=energy_bins[2].value
for i,E in enumerate(energy_bins[1:-4]):
    E1=energy_bins[i].value
    E2=energy_bins[i+1].value
    print "Energy band= E1=" +str(E1)+ "et E2="+str(E2)
    excess = SkyMapCollection.read("fov_bg_maps"+str(E1)+"_"+str(E2)+"_TeV.fits")["excess"]
    #Ecrire la map en fits pour pouvoir load the image ca ca pue aussi faudrait pourvoir direct cree l objet data sans ecrire...
    excess.write("excess_maps"+str(E1)+"_"+str(E2)+"_TeV.fits", clobber=True)
    data = fits.open("excess_maps"+str(E1)+"_"+str(E2)+"_TeV.fits")
    load_image(1, data)
    #image_data()

    psf_file_SgrA = Table.read("psf_table_SgrA_" + str(E1) + "_" + str(E2) + ".fits")
    psf_file_G0p9 = Table.read("psf_table_G0p9"+str(E1)+"_"+str(E2)+".fits")

    header = excess.to_image_hdu().header
    #Todo:voir pour la normalisation normalement il le fait tout seul mais pas sur...
    psf_image_SgrA = fill_acceptance_image(header, excess.center(), psf_file_SgrA["theta"].to("deg"),
                                      psf_file_SgrA["psf_value"].data, psf_file_SgrA["theta"].to("deg")[-1])
    source_center_SgrA = SkyCoord(359.9442,-0.0462, unit='deg',frame="galactic")
    #source_center_SgrA = SkyCoord.from_name("SgrA*")
    source_center_G0p9 = SkyCoord(0.868, 0.075, unit='deg', frame="galactic")
    psf_image_G0p9 = fill_acceptance_image(header, source_center_G0p9, psf_file_G0p9["theta"].to("deg"),
                                      psf_file_G0p9["psf_value"].data, psf_file_G0p9["theta"].to("deg")[-1])
    psf_image_SgrA.writeto("psf_image_SgrA_" + str(E1) + "_" + str(E2) + ".fits", clobber=True)
    psf_image_G0p9.writeto("psf_image_G0p9_" + str(E1) + "_" + str(E2) + ".fits", clobber=True)
    load_psf("psf_SgrA", "psf_image_SgrA_" + str(E1) + "_" + str(E2) + ".fits")
    load_psf("psf_G0p9", "psf_image_G0p9_" + str(E1) + "_" + str(E2) + ".fits")
    #set_psf("psf2")
    #mygaus_SgrA = Gauss2D("g2")
    mygaus_SgrA =normgauss2dint("g2")
    large_gaus = Gauss2D("g2")
    mygaus_SgrA.xpos,mygaus_SgrA.ypos=skycoord_to_pixel(source_center_SgrA, excess.wcs)
    large_gaus.xpos,large_gaus.ypos=skycoord_to_pixel(source_center_SgrA, excess.wcs)
    large_gaus.fwhm = 20
    #mygaus_G0p9 = Gauss2D("g2")
    mygaus_G0p9 =normgauss2dint("g2")
    mygaus_G0p9.xpos,mygaus_G0p9.ypos=skycoord_to_pixel(source_center_G0p9, excess.wcs)
    #On voit que c est pas bien centre ca donne de la merde mais par contre image_psf ok donc c'est un problem avec le model
    #set_source(mygaus2)
    CS_map=SkyMap.read("CStot.fits")
    cs_reproj=CS_map.reproject(excess)
    cs_reproj.data[np.where(np.isnan(cs_reproj.data))]=0
    cs_reproj.write("cs_map_reproj.fits", clobber=True)
    load_table_model("CS","cs_map_reproj.fits")
    #cs_reproj.write("cs_map_reproj.fits", clobber=True)
    #load_table_model("CS","cs_map_reproj.fits")
    #set_full_model(psf_G0p9(mygaus_G0p9)+psf_SgrA(mygaus_SgrA+CS))
    #define_asymgaussiansource("GCext",359.9442,-0.04622, 0.1,0.1,1.57,1.0)
    #freeze(GCext.xpos)
    #freeze(GCext.ypos)
    #set_full_model(psf_G0p9(mygaus_G0p9)+psf_SgrA(mygaus_SgrA)+GCext)
    set_full_model(psf_G0p9(mygaus_G0p9)+psf_SgrA(mygaus_SgrA))
    #set_full_model(psf_SgrA(mygaus_SgrA))
    #set_full_model(psf_G0p9(mygaus_G0p9)+psf_SgrA(mygaus_SgrA+large_gaus*CS))
    #set_full_model(psf_G0p9(mygaus_G0p9)+psf_SgrA(mygaus_SgrA)+large_gaus*CS)
    set_stat("cstat")
    #set_method("neldermead")

    """
    freeze(mygaus_G0p9.xpos)
    freeze(mygaus_G0p9.ypos)
    freeze(mygaus_SgrA.xpos)
    freeze(mygaus_SgrA.ypos)
    #freeze(mygaus_G0p9.fwhm)
    #freeze(mygaus_SgrA.fwhm)
    freeze(large_gaus.xpos)
    freeze(large_gaus.ypos)
    """
    large_gaus.ampl =1
    freeze(large_gaus.ampl)
    freeze(large_gaus.xpos)
    freeze(large_gaus.ypos)
    fit()
    save_resid("residual_"+str(E1)+"_"+str(E2)+"_TeV.fits", clobber=True)
    #image_psf()
    #image_model()
    #image_fit()
    #image_resid()
    binsize_carte=np.fabs(data[1].header["CDELT1"])
    x_SgrA,y_SgrA=mygaus_SgrA.xpos.val, mygaus_SgrA.ypos.val
    fwhm_SgrA=mygaus_SgrA.fwhm.val*binsize_carte
    l_SgrA=pixel_to_skycoord(x_SgrA,y_SgrA,excess.wcs).l
    b_SgrA=pixel_to_skycoord(x_SgrA,y_SgrA,excess.wcs).b


    x_G0p9,y_G0p9=mygaus_G0p9.xpos.val, mygaus_G0p9.ypos.val
    fwhm_G0p9=mygaus_G0p9.fwhm.val*binsize_carte
    l_G0p9=pixel_to_skycoord(x_G0p9,y_G0p9,excess.wcs).l
    b_G0p9=pixel_to_skycoord(x_G0p9,y_G0p9,excess.wcs).b

    print("Pour SgrA*, les coord sont (l,b)= ("+str(l_SgrA.value)+","+str(b_SgrA.value)+") deg et la largeur de la source: "+str(fwhm_SgrA)+" deg")

    print("Pour G0p9, les coord sont (l,b)= ("+str(l_G0p9.value)+","+str(b_G0p9.value)+") deg et la largeur de la source: "+str(fwhm_G0p9)+" deg")

"""
txt_psf = open("psf.txt", "w")
for x, y in zip(psf_file["theta"].data, psf_file["psf_value"].data):
    txt_psf.write(str(x) + "\t" + str(y) + "\n")
txt_psf.close()
# load_table_model("psf",psf_file)
load_table_model("psf", "psf.txt")
set_model(psf)
# set_coord("physical")
set_stat("cash")
image_fit()

from astropy.coordinates import Angle

header = excess.to_image_hdu().header
norm1 = excess.data.max() / psf_file["psf_value"].data.max()
psf_image = fill_acceptance_image(header, excess.center(), psf_file["theta"].to("deg"),
                                  psf_file["psf_value"].data * norm1, psf_file["theta"].to("deg")[-1])
psf_image.writeto("psf_image.fits", clobber=True)
load_psf("psf", "psf_image.fits")
set_psf("psf")
from sherpa.models import Gauss2D
mygaus = Gauss2D("g2")
set_source(mypsf(mygaus))
fit()


from sherpa.instrument import PSFModel
mypsf=PSFModel("psf2",psf_image.data)
set_psf(mypsf)
from sherpa.models import Gauss2D
mygaus = Gauss2D("g2")
#Permet de directement convolue model et psf et d affecter du coup une psf difference a differente source
set_source(mypsf(mygaus))
fit()


set_model(psf)

image_fit()

"""
"""
#Rebin histo psfvalue
rebin=Angle(excess.meta["CDELT2"],excess.meta["CUNIT2"])
theta=psf_file["theta"].to(rebin.unit)
psf_value=psf_file["psf_value"]
binsizepsf=theta[1]-theta[0]
theta_max=theta.max()
bintab=Angle(np.arange(0, theta_max.value+rebin.value, rebin.value), rebin.unit)
theta_tab=(bintab[0:-1]+bintab[1:])/2
Nbin_group=rebin/binsizepsf
index=np.arange(0,len(theta)+Nbin_group,Nbin_group).astype(int)
psfrebin=np.array([])
for i in range(len(index[:-1])):
    value=np.sum(psf_value[index[i]:index[i+1]]*2*np.pi*theta[index[i]:index[i+1]]*binsizepsf).decompose()
    #diviser par l angle solide du nouveau bin pour que integrale de dP/dOme*dOme soit normalise
    value/=2*np.pi*theta_tab[i]*rebin
    psfrebin = np.append(psfrebin,value.to("1/sr").value)
psfrebin=Quantity(psfrebin, psf_value.unit)
#TODO: j ai un pb sur le premier bin il est trop bas peut venir du fait que premier bin en theta=0
norm=np.sum(psf_value*2*np.pi*theta*binsizepsf).decompose()
psfrebin/=norm
print np.sum(psfrebin*rebin*2*np.pi*theta_tab).decompose()
txt_psf = open("psf.txt", "w")
for y,x in zip(psfrebin.value*2*np.pi*theta_tab.value*rebin.value, theta_tab.value):
    txt_psf.write(str(x) + "\t" + str(y) + "\n")
txt_psf.close()


load_psf("psf", "psf.txt")
psf.size = 75
psf.center = 0
psf.radial = 1
psf.norm = 0
set_psf(psf)
image_psf()

table=Table()
aa = Column(psfrebin.value*2*np.pi*theta_tab.value*rebin.value, name='value')
table.add_column(aa)
table.meta["kernel"]="psf.fits"
table.meta[""]= 75
table.meta[""]= 0
table.meta[""]= 1
table.meta[""]= 0
table.write("psf.fits")

#set_stat("cash")
fit()
image_fit()
"""
