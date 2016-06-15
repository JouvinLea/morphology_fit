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
import astropy.units as u
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
# E2=energy_bins[2].value
#for i_E, E in enumerate(energy_bins[0:-3]):
for i_E, E in enumerate(energy_bins[0:-5]):
    #E1 = energy_bins[i_E].value
    #E2 = energy_bins[i_E+1].value
    E1 = energy_bins[0].value
    E2 = energy_bins[3].value
    print "Energy band= E1=" + str(E1) + "et E2=" + str(E2)
    # on = SkyMapCollection.read("fov_bg_maps"+str(E1)+"_"+str(E2)+"_TeV.fits")["excess"]
    on = SkyMapCollection.read("fov_bg_maps" + str(E1) + "_" + str(E2) + "_TeV.fits")["counts"]
    bkgmap = SkyMapCollection.read("fov_bg_maps" + str(E1) + "_" + str(E2) + "_TeV.fits")["bkg"]
    exp = SkyMapCollection.read("fov_bg_maps" + str(E1) + "_" + str(E2) + "_TeV.fits")["exposure"]
    # source_J1745_303 = SkyCoord(358.76,-0.51, unit='deg',frame="galactic")
    # On descend un peu en lattitude la region d'exclusion

    on.write("on_maps" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    # bkgmap.data=bkgmap.data/bkgmap.data.sum()
    bkgmap.write("off_maps" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    exp.write("exp_maps" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)

    data = fits.open("on_maps" + str(E1) + "_" + str(E2) + "_TeV.fits")
    load_image(1, data)
    # image_data()

    psf_file_SgrA = Table.read("psf_table_SgrA_" + str(E1) + "_" + str(E2) + ".fits")
    psf_file_G0p9 = Table.read("psf_table_G0p9" + str(E1) + "_" + str(E2) + ".fits")

    header = on.to_image_hdu().header
    # Todo:voir pour la normalisation normalement il le fait tout seul mais pas sur...
    psf_image_SgrA = fill_acceptance_image(header, on.center(), psf_file_SgrA["theta"].to("deg"),
                                           psf_file_SgrA["psf_value"].data, psf_file_SgrA["theta"].to("deg")[-1])
    source_center_SgrA = SkyCoord(359.9442, -0.0462, unit='deg', frame="galactic")
    # source_center_SgrA = SkyCoord.from_name("SgrA*")
    source_center_G0p9 = SkyCoord(0.868, 0.075, unit='deg', frame="galactic")
    psf_image_G0p9 = fill_acceptance_image(header, on.center(), psf_file_G0p9["theta"].to("deg"),
                                           psf_file_G0p9["psf_value"].data, psf_file_G0p9["theta"].to("deg")[-1])
    psf_image_SgrA.writeto("psf_image_SgrA_" + str(E1) + "_" + str(E2) + ".fits", clobber=True)
    psf_image_G0p9.writeto("psf_image_G0p9_" + str(E1) + "_" + str(E2) + ".fits", clobber=True)
    load_psf("psf_SgrA", "psf_image_SgrA_" + str(E1) + "_" + str(E2) + ".fits")
    load_psf("psf_G0p9", "psf_image_G0p9_" + str(E1) + "_" + str(E2) + ".fits")

    # modele gauss pour sgrA centre sur SgrA
    mygaus_SgrA = normgauss2dint("SgrA")
    mygaus_SgrA.xpos, mygaus_SgrA.ypos = skycoord_to_pixel(source_center_SgrA, on.wcs)
    mygaus_SgrA.xpos.val += 0.5
    mygaus_SgrA.ypos.val += 0.5
    # Modele marge gaussienne a multiplie avec CS centre sur SgrA
    large_gaus = Gauss2D("Gauss_to_CS")
    large_gaus.xpos, large_gaus.ypos = skycoord_to_pixel(source_center_SgrA, on.wcs)
    #large_gaus.fwhm = 100
    #modele asymetric large scale gauss centre sur SgrA
    asym_gaus = Gauss2D("Large_scale")
    asym_gaus.xpos, asym_gaus.ypos = skycoord_to_pixel(source_center_SgrA, on.wcs)
    #modele symetric central gauss centre sur SgrA
    central_gauss = Gauss2D("central_gauss")
    central_gauss.xpos, central_gauss.ypos = skycoord_to_pixel(source_center_SgrA, on.wcs)
    # modele gauss pour G0p9 centre sur G0p9
    mygaus_G0p9 = normgauss2dint("G0p9")
    mygaus_G0p9.xpos, mygaus_G0p9.ypos = skycoord_to_pixel(source_center_G0p9, on.wcs)
    mygaus_G0p9.xpos.val += 0.5
    mygaus_G0p9.ypos.val += 0.5

    #Arc_source
    arc_source=SkyCoord(0.130, -0.139, unit='deg', frame="galactic")
    mygaus_arcsource = Gauss2D("Arc Source")
    mygaus_arcsource.xpos, mygaus_arcsource.ypos = skycoord_to_pixel(arc_source, on.wcs)

    # Modele pour le bkg base sur la carte de fond
    load_table_model("bkg", "off_maps" + str(E1) + "_" + str(E2) + "_TeV.fits")
    bkg.ampl = 1
    load_table_model("exposure", "exp_maps" + str(E1) + "_" + str(E2) + "_TeV.fits")
    exposure.ampl = 1
    # Modele de CS
    CS_map = SkyMap.read("CStot.fits")
    cs_reproj = CS_map.reproject(on)
    cs_reproj.data[np.where(np.isnan(cs_reproj.data))] = 0
    cs_reproj.data[np.where(cs_reproj.data < 30)] = 0
    cs_reproj.write("cs_map_reproj.fits", clobber=True)
    load_table_model("CS", "cs_map_reproj.fits")

    set_full_model(bkg + psf_SgrA(mygaus_SgrA)+psf_G0p9(mygaus_G0p9))
    #set_full_model(bkg + psf_SgrA(exposure*mygaus_SgrA)+psf_G0p9(exposure*mygaus_G0p9))
    #set_full_model(bkg + psf_SgrA(exposure*(mygaus_SgrA+large_gaus * CS))+psf_G0p9(exposure*mygaus_G0p9))
    #set_full_model(bkg + psf_SgrA(mygaus_SgrA + large_gaus * CS) + psf_G0p9(mygaus_G0p9))
    #set_full_model(bkg + psf_SgrA(mygaus_SgrA + large_gaus * CS + central_gauss) + psf_G0p9(mygaus_G0p9))
    #set_full_model(bkg + psf_SgrA(mygaus_SgrA + large_gaus * CS + mygaus_arcsource) + psf_G0p9(mygaus_G0p9))
    mygaus_SgrA.fwhm = 1
    freeze(mygaus_SgrA.fwhm)
    mygaus_G0p9.fwhm = 1
    freeze(mygaus_G0p9.fwhm)
    mygaus_arcsource.fwhm =1
    freeze(mygaus_arcsource.fwhm)


    freeze(mygaus_G0p9.xpos)
    freeze(mygaus_G0p9.ypos)
    freeze(mygaus_SgrA.xpos)
    freeze(mygaus_SgrA.ypos)
    freeze(mygaus_arcsource.xpos)
    freeze(mygaus_arcsource.ypos)
    #thaw(mygaus_arcsource.ellip)
    freeze(asym_gaus.xpos)
    freeze(asym_gaus.ypos)
    thaw(asym_gaus.ellip)

    #large_gaus.fwhm = 95
    #freeze(large_gaus.fwhm)
    large_gaus.ampl = 1
    freeze(large_gaus.ampl)
    freeze(large_gaus.xpos)
    freeze(large_gaus.ypos)
    #central_gauss.fwhm =10
    #freeze(central_gauss.fwhm)
    freeze(central_gauss.xpos)
    freeze(central_gauss.ypos)


    #freeze(bkg.ampl)
    freeze(exposure.ampl)


    source_J1745_303 = SkyCoord(358.76, -0.6, unit='deg', frame="galactic")
    source_J1745_303_xpix, source_J1745_303_ypix = skycoord_to_pixel(source_J1745_303, on.wcs)
    pix_deg = on.to_image_hdu().header["CDELT2"]
    #radius = 0.4 / pix_deg
    width=100
    height=80
    name_region = "box(" + str(source_J1745_303_xpix+20) + "," + str(source_J1745_303_ypix-20) + "," + str(width) + "," + str(height) +")"
    lat=1.6/ pix_deg#Pour aller a plus et -0.8 as did Anne
    lon=4 / pix_deg#Pour aller a plus ou moins 2deg as did Anne
    x_pix_SgrA=skycoord_to_pixel(source_center_SgrA, on.wcs)[0]
    y_pix_SgrA=skycoord_to_pixel(source_center_SgrA, on.wcs)[1]
    name_interest = "box(" + str(x_pix_SgrA) + "," + str(y_pix_SgrA) + "," + str(lon) + "," + str(lat) +")"
    notice2d(name_interest)
    ignore2d(name_region)

    set_stat("cstat")
    set_method("neldermead")

    list_src = [psf_SgrA(mygaus_SgrA),psf_G0p9(mygaus_G0p9),psf_SgrA(asym_gaus),psf_SgrA(large_gaus * CS),psf_SgrA(central_gauss), psf_SgrA(mygaus_arcsource)]
    print mygaus_SgrA

    model = bkg
    for i_src,src in enumerate(list_src):
        model += src
        set_full_model(model)
        fit()
        result=get_fit_results()
        import IPython; IPython.embed()
        save_resid("residual_step_"+str(i_src)+"_" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)



