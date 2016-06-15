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
    large_gaus.fwhm = 100
    central_gauss = Gauss2D("central_gauss")
    central_gauss.xpos, central_gauss.ypos = skycoord_to_pixel(source_center_SgrA, on.wcs)
    central_gauss.fwhm = 100
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

    #set_full_model(bkg + psf_SgrA(mygaus_SgrA)+psf_G0p9(mygaus_G0p9))
    #set_full_model(bkg + psf_SgrA(exposure*mygaus_SgrA)+psf_G0p9(exposure*mygaus_G0p9))
    #set_full_model(bkg + psf_SgrA(exposure*(mygaus_SgrA+large_gaus * CS))+psf_G0p9(exposure*mygaus_G0p9))
    set_full_model(bkg + psf_SgrA(mygaus_SgrA + large_gaus * CS) + psf_G0p9(mygaus_G0p9))
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



    fit()
    import IPython; IPython.embed()
    #save_resid("residual_avec_CS_" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    #save_resid("residual_avec_CS_bkg_ampl_free_" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    #save_resid("residual_sans_CS_bkg_ampl_free_" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    #save_resid("residual_sans_CS_" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    #save_resid("residual_sans_CS_bkg_ampl_free_exposure_" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    #save_resid("residual_avec_CS_bkg_ampl_free_plus_central_gauss_" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    #save_resid("residual_avec_CS_bkg_ampl_free_plus_central_gauss_fwhmfix_" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    #save_resid("residual_avec_CS_bkg_ampl_free_plus_central_gauss_CSgauss_fwhmfix_" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    #save_resid("residual_avec_CS_bkg_ampl_free_plus_arc_source_" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    #save_resid("residual_avec_CS_bkg_ampl_free_plus_arc_source_CSgauss_fwhmfix_" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    #save_resid("residual_avec_CS_bkg_ampl_free_plus_arc_source_fwhmfix_" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    #save_resid("residual_avec_CS_bkg_ampl_free_exposure_" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    # save_model("model_"+str(E1)+"_"+str(E2)+"_TeV.fits", clobber=True)
    # save_data("data_"+str(E1)+"_"+str(E2)+"_TeV.fits", clobber=True)
    shape = np.shape(on.data)
    mask = get_data().mask.reshape(shape)
    #map_data = np.flipud(get_data().y.reshape(shape) * mask)
    #model = np.flipud(get_model()(get_data().x0, get_data().x1).reshape(shape) * mask)
    map_data = get_data().y.reshape(shape) * mask
    model = get_model()(get_data().x0, get_data().x1).reshape(shape) * mask
    resid = map_data - model
    coord = on.coordinates()
    # data_on=SkyMap.read("data_"+str(E1)+"_"+str(E2)+"_TeV.fits", clobber=True)
    # resid=SkyMap.read("residual_"+str(E1)+"_"+str(E2)+"_TeV.fits", clobber=True)
    # model=SkyMap.read("model_"+str(E1)+"_"+str(E2)+"_TeV.fits", clobber=True)
    i_b = np.where((coord.b[:, 0] < on.center().b + 0.25 * u.deg) & (coord.b[:, 0] > on.center().b - 0.25 * u.deg))[0]
    # resid.data[np.isnan(resid.data)]=0
    # model.data[np.isnan(model.data)]=0
    # data_on.data[np.isnan(data_on.data)]=0
    # on.data[np.isnan(model.data)]=0

    # TODO Faire un tableau qui pour chaque longitude compte le nombre de pixel cumule en b pour cette longitude
    npix_l=np.sum(np.flipud(mask[i_b,:]), axis=0)
    profile_l_model = np.sum(model[i_b, :], axis=0)/npix_l
    # profile_l_on=np.sum(on.data[i_b,:],axis=0)
    profile_l_on = np.sum(map_data[i_b, :], axis=0)/npix_l
    #Ca donne des coups par arcmin2 car on prend en compte qu on ne cumula pas le meme nombre de pixel pour chaque
    # longitude vu qu il y a des regions d exclusions
    profile_l_resid = np.sum(resid[i_b, :], axis=0)/npix_l
    err_l = np.sqrt(profile_l_on/npix_l)
    l = coord.l[0, :]
    l.value[np.where(l > 180 * u.deg)] = l.value[np.where(l > 180 * u.deg)] - 360

    rebin=3
    i_rebin = np.arange(0, len(profile_l_on), rebin)
    resid_l_rebin = np.array([])
    l_rebin = np.array([])
    err_l_rebin = np.array([])
    for i in range(len(i_rebin[:-1])):
        resid_l_rebin = np.append(resid_l_rebin, np.mean(profile_l_resid[i_rebin[i]:i_rebin[i + 1]]))
        l_rebin = np.append(l_rebin, np.mean(l[i_rebin[i]:i_rebin[i + 1]]))
        err_l_rebin = np.append(err_l_rebin, np.mean(err_l[i_rebin[i]:i_rebin[i + 1]]))
    resid_l_rebin = np.append(resid_l_rebin, np.mean(profile_l_resid[i_rebin[i + 1]:]))
    l_rebin = np.append(l_rebin, np.mean(l[i_rebin[i + 1]:]))
    err_l_rebin = np.append(err_l_rebin, np.mean(err_l[i_rebin[i + 1]:]))

    #fig = pt.figure(i_E)
    fig = pt.figure(1)
    ax = fig.add_subplot(2, 1, 1)
    pt.plot(l.value, profile_l_model, label="model")
    pt.plot(l.value, profile_l_on, label="on data")
    pt.xlim(-1.5, 1.5)
    pt.gca().invert_xaxis()
    pt.legend()
    ax = fig.add_subplot(2, 1, 2)
    #pt.errorbar(l.value, profile_l_resid, yerr=err_l, linestyle='None', marker="o")
    pt.errorbar(l_rebin.value, resid_l_rebin, yerr=err_l_rebin, linestyle='None', marker="o", label= str("%.2f"%E1) + "_" + str("%.2f"%E2) + "_TeV")
    pt.axhline(y=0, color='red', linewidth=2)
    pt.legend()
    pt.ylabel("residual")
    pt.xlabel("longitude (degrees)")
    pt.xlim(-1.5, 1.5)
    pt.gca().invert_xaxis()
    #pt.savefig("plot/profile_longitude_avec_CS_" + str("%.2f"%E1) + "_" + str("%.2f"%E2) + "_TeV.jpg")
    #pt.savefig("plot/profile_longitude_avec_CS_bkg_ampl_free_" + str(E1) + "_" + str(E2) + "_TeV.jpg")
    #pt.savefig("plot/profile_longitude_sans_CS_bkg_ampl_free_"+str(E1)+"_"+str(E2)+"_TeV.jpg")
    #pt.savefig("plot/profile_longitude_sans_CS_"+str(E1)+"_"+str(E2)+"_TeV.jpg")
    #pt.savefig("plot/profile_longitude_sans_CS_bkg_ampl_free_exposure_"+str(E1)+"_"+str(E2)+"_TeV.jpg")
    #pt.savefig("plot/profile_longitude_avec_CS_bkg_ampl_free_plus_central_gauss_" + str(E1) + "_" + str(E2) + "_TeV.jpg")
    #pt.savefig("plot/profile_longitude_avec_CS_bkg_ampl_free_plus_central_gauss_fwhmfix_" + str(E1) + "_" + str(E2) + "_TeV.jpg")
    #pt.savefig("plot/profile_longitude_avec_CS_bkg_ampl_free_plus_central_gauss_CSgauss_fwhmfix_" + str(E1) + "_" + str(E2) + "_TeV.jpg")
    #pt.savefig("plot/profile_longitude_avec_CS_bkg_ampl_free_plus_arc_source_" + str(E1) + "_" + str(E2) + "_TeV.jpg")
    #pt.savefig("plot/profile_longitude_avec_CS_bkg_ampl_free_plus_arc_source_CSgauss_fwhmfix_" + str(E1) + "_" + str(E2) + "_TeV.jpg")
    #pt.savefig("plot/profile_longitude_avec_CS_bkg_ampl_free_plus_arc_source_fwhmfix_" + str(E1) + "_" + str(E2) + "_TeV.jpg")
    #pt.savefig("plot/profile_longitude_avec_CS_bkg_ampl_free_exposure_"+str(E1)+"_"+str(E2)+"_TeV.jpg")
    # Ici coord l deja modif cad que les valeur sup a 360 sont mis negatives avec - 360
    l_center = on.center().l
    if l_center > 180 * u.deg:
        l_center = l_center - 360 * u.deg
    i_l = np.where((l < l_center + 1.5 * u.deg) & (l > l_center - 1.5 * u.deg))[0]
    npix_b=np.sum(np.flipud(mask[:,i_l]), axis=1)
    profile_b_model = np.sum(model[:, i_l], axis=1)/npix_b
    # profile_b_on = np.sum(on.data[:,i_l],axis=1)
    profile_b_on = np.sum(map_data[:, i_l], axis=1)/npix_b
    profile_b_resid = np.sum(resid[:, i_l], axis=1)/npix_b
    err_b = np.sqrt(profile_b_on/npix_b)

    i_rebin = np.arange(0, len(profile_b_on), rebin)
    resid_b_rebin = np.array([])
    b_rebin = np.array([])
    err_b_rebin = np.array([])
    for i in range(len(i_rebin[:-1])):
        resid_b_rebin = np.append(resid_b_rebin, np.mean(profile_b_resid[i_rebin[i]:i_rebin[i + 1]]))
        b_rebin = np.append(b_rebin, np.mean(coord.b[:, 0][i_rebin[i]:i_rebin[i + 1]]))
        err_b_rebin = np.append(err_b_rebin, np.mean(err_b[i_rebin[i]:i_rebin[i + 1]]))
    resid_b_rebin = np.append(resid_b_rebin, np.mean(profile_b_resid[i_rebin[i + 1]:]))
    b_rebin = np.append(b_rebin, np.mean(coord.b[:, 0][i_rebin[i + 1]:]))
    err_b_rebin = np.append(err_b_rebin, np.mean(err_b[i_rebin[i + 1]:]))

    #fig = pt.figure(i_E + 100)
    fig =pt.figure(100)
    ax = fig.add_subplot(2, 1, 1)
    pt.plot(coord.b[:, 0].value, profile_b_model, label="model")
    pt.plot(coord.b[:, 0].value, profile_b_on, label="on data")
    pt.xlim(-1, 1)
    pt.legend()
    ax = fig.add_subplot(2, 1, 2)
    #pt.errorbar(coord.b[:, 0].value, profile_b_resid, yerr=err_b, linestyle='None', marker="o")
    pt.errorbar(b_rebin.value, resid_b_rebin, yerr=err_b_rebin, linestyle='None', marker="o", label= str("%.2f"%E1) + "_" + str("%.2f"%E2) + "_TeV")
    pt.axhline(y=0, color='red', linewidth=2)
    pt.legend()
    pt.ylabel("residual")
    pt.xlabel("latitude (degrees)")
    pt.xlim(-1, 1)
    #pt.savefig("plot/profile_lattitude_avec_CS_" + str("%.2f"%E1) + "_" + str("%.2f"%E2) + "_TeV.jpg")
    #pt.savefig("plot/profile_lattitude_avec_CS_bkg_ampl_free_" + str(E1) + "_" + str(E2) + "_TeV.jpg")
    #pt.savefig("plot/profile_lattitude_sans_CS_bkg_ampl_free_"+str(E1)+"_"+str(E2)+"_TeV.jpg")
    #pt.savefig("plot/profile_lattitude_sans_CS_"+str(E1)+"_"+str(E2)+"_TeV.jpg")
    #pt.savefig("plot/profile_lattitude_sans_CS_bkg_ampl_free_exposure_"+str(E1)+"_"+str(E2)+"_TeV.jpg")
    #pt.savefig("plot/profile_lattitude_avec_CS_bkg_ampl_free_plus_arc_source_" + str(E1) + "_" + str(E2) + "_TeV.jpg")
    #pt.savefig("plot/profile_lattitude_avec_CS_bkg_ampl_free_plus_arc_source_CSgauss_fwhmfix_" + str(E1) + "_" + str(E2) + "_TeV.jpg")
    #pt.savefig("plot/profile_lattitude_avec_CS_bkg_ampl_free_plus_arc_source_fwhmfix_" + str(E1) + "_" + str(E2) + "_TeV.jpg")
    #pt.savefig("plot/profile_lattitude_avec_CS_bkg_ampl_free_plus_central_gauss_" + str(E1) + "_" + str(E2) + "_TeV.jpg")
    #pt.savefig("plot/profile_lattitude_avec_CS_bkg_ampl_free_plus_central_gauss_fwhmfix_" + str(E1) + "_" + str(E2) + "_TeV.jpg")
    #pt.savefig("plot/profile_lattitude_avec_CS_bkg_ampl_free_plus_central_gauss_CSgauss_fwhmfix_" + str(E1) + "_" + str(E2) + "_TeV.jpg")
    #pt.savefig("plot/profile_lattitude_avec_CS_bkg_ampl_free_exposure_"+str(E1)+"_"+str(E2)+"_TeV.jpg")
    # save_resid("residual_avec_CS_method_cstat_neldermead_fwmh_Sgra_etG0p9_fix_"+str(E1)+"_"+str(E2)+"_TeV.fits", clobber=True)
    # save_resid("residual_avec_CS_method_cstat_neldermead_"+str(E1)+"_"+str(E2)+"_TeV.fits", clobber=True)
    # save_resid("residual_sans_CS_method_cstat_neldermead_"+str(E1)+"_"+str(E2)+"_TeV.fits", clobber=True)
    # save_resid("residual_"+str(E1)+"_"+str(E2)+"_TeV.fits", clobber=True)
    # set_full_model(psf_SgrA(mygaus_SgrA+large_gaus*CS)+psf_G0p9(mygaus_G0p9))
    # image_psf()
    # image_model()
    # image_fit()
    # image_resid()
    binsize_carte = np.fabs(data[1].header["CDELT1"])
    x_SgrA, y_SgrA = mygaus_SgrA.xpos.val, mygaus_SgrA.ypos.val
    fwhm_SgrA = mygaus_SgrA.fwhm.val * binsize_carte
    l_SgrA = pixel_to_skycoord(x_SgrA, y_SgrA, on.wcs).l
    b_SgrA = pixel_to_skycoord(x_SgrA, y_SgrA, on.wcs).b

    x_G0p9, y_G0p9 = mygaus_G0p9.xpos.val, mygaus_G0p9.ypos.val
    fwhm_G0p9 = mygaus_G0p9.fwhm.val * binsize_carte
    l_G0p9 = pixel_to_skycoord(x_G0p9, y_G0p9, on.wcs).l
    b_G0p9 = pixel_to_skycoord(x_G0p9, y_G0p9, on.wcs).b

    print("Pour SgrA*, les coord sont (l,b)= (" + str(l_SgrA.value) + "," + str(
        b_SgrA.value) + ") deg et la largeur de la source: " + str(fwhm_SgrA) + " deg")

    print("Pour G0p9, les coord sont (l,b)= (" + str(l_G0p9.value) + "," + str(
        b_G0p9.value) + ") deg et la largeur de la source: " + str(fwhm_G0p9) + " deg")
