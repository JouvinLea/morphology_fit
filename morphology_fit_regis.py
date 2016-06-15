from sherpa.astro.ui import *
from astropy.io import fits
from astropy.table import Table, join
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
from utilities import *
from matplotlib.backends.backend_pdf import PdfPages

pt.ion()
energy_bins = EnergyBounds.equal_log_spacing(0.5, 100, 5, 'TeV')


#for i_E, E in enumerate(energy_bins[0:-3]):
for i_E, E in enumerate(energy_bins[0:-5]):
    #E1 = energy_bins[i_E].value
    #E2 = energy_bins[i_E+1].value
    E1 = energy_bins[0].value
    E2 = energy_bins[3].value
    print "Energy band= E1=" + str(E1) + "et E2=" + str(E2)

    #load Data
    on = SkyMapCollection.read("fov_bg_maps" + str(E1) + "_" + str(E2) + "_TeV.fits")["counts"]
    on.write("on_maps" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    data = fits.open("on_maps" + str(E1) + "_" + str(E2) + "_TeV.fits")
    load_image(1, data)

    #load bkg model
    bkgmap = SkyMapCollection.read("fov_bg_maps" + str(E1) + "_" + str(E2) + "_TeV.fits")["bkg"]
    bkgmap.write("off_maps" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    load_table_model("bkg", "off_maps" + str(E1) + "_" + str(E2) + "_TeV.fits")
    bkg.ampl = 1
    #freeze(bkg.ampl)

    #load exposure model
    exp = SkyMapCollection.read("fov_bg_maps" + str(E1) + "_" + str(E2) + "_TeV.fits")["exposure"]
    exp.write("exp_maps" + str(E1) + "_" + str(E2) + "_TeV.fits", clobber=True)
    load_table_model("exposure", "exp_maps" + str(E1) + "_" + str(E2) + "_TeV.fits")
    exposure.ampl = 1
    freeze(exposure.ampl)

    #load psf model
    psf_file_SgrA = Table.read("psf_table_SgrA_" + str(E1) + "_" + str(E2) + ".fits")
    psf_file_G0p9 = Table.read("psf_table_G0p9" + str(E1) + "_" + str(E2) + ".fits")
    header = on.to_image_hdu().header
    psf_image_SgrA = fill_acceptance_image(header, on.center(), psf_file_SgrA["theta"].to("deg"),
                                           psf_file_SgrA["psf_value"].data, psf_file_SgrA["theta"].to("deg")[-1])
    psf_image_G0p9 = fill_acceptance_image(header, on.center(), psf_file_G0p9["theta"].to("deg"),
                                           psf_file_G0p9["psf_value"].data, psf_file_G0p9["theta"].to("deg")[-1])
    psf_image_SgrA.writeto("psf_image_SgrA_" + str(E1) + "_" + str(E2) + ".fits", clobber=True)
    psf_image_G0p9.writeto("psf_image_G0p9_" + str(E1) + "_" + str(E2) + ".fits", clobber=True)
    load_psf("psf_SgrA", "psf_image_SgrA_" + str(E1) + "_" + str(E2) + ".fits")
    load_psf("psf_G0p9", "psf_image_G0p9_" + str(E1) + "_" + str(E2) + ".fits")

    # Modele de CS
    CS_map = SkyMap.read("CStot.fits")
    cs_reproj = CS_map.reproject(on)
    cs_reproj.data[np.where(np.isnan(cs_reproj.data))] = 0
    cs_reproj.data[np.where(cs_reproj.data < 30)] = 0
    cs_reproj.write("cs_map_reproj.fits", clobber=True)
    load_table_model("CS", "cs_map_reproj.fits")

    # modele gauss pour sgrA centre sur SgrA
    source_center_SgrA = SkyCoord(359.9442, -0.0462, unit='deg', frame="galactic")
    mygaus_SgrA = normgauss2dint("SgrA")
    mygaus_SgrA.xpos, mygaus_SgrA.ypos = skycoord_to_pixel(source_center_SgrA, on.wcs)
    mygaus_SgrA.xpos.val += 0.5
    mygaus_SgrA.ypos.val += 0.5
    freeze(mygaus_SgrA.xpos)
    freeze(mygaus_SgrA.ypos)
    mygaus_SgrA.fwhm = 1
    freeze(mygaus_SgrA.fwhm)

    # modele gauss pour G0p9 centre sur G0p9
    source_center_G0p9 = SkyCoord(0.868, 0.075, unit='deg', frame="galactic")
    mygaus_G0p9 = normgauss2dint("G0p9")
    mygaus_G0p9.xpos, mygaus_G0p9.ypos = skycoord_to_pixel(source_center_G0p9, on.wcs)
    mygaus_G0p9.xpos.val += 0.5
    mygaus_G0p9.ypos.val += 0.5
    freeze(mygaus_G0p9.xpos)
    freeze(mygaus_G0p9.ypos)
    mygaus_G0p9.fwhm = 1
    freeze(mygaus_G0p9.fwhm)

    #modele asymetric large scale gauss centre sur SgrA
    asym_gaus = Gauss2D("Large_scale")
    asym_gaus.xpos, asym_gaus.ypos = skycoord_to_pixel(source_center_SgrA, on.wcs)
    freeze(asym_gaus.xpos)
    freeze(asym_gaus.ypos)
    thaw(asym_gaus.ellip)

    # Modele large gaussienne  multiplie avec CS centre sur SgrA
    large_gaus = Gauss2D("Gauss_to_CS")
    large_gaus.xpos, large_gaus.ypos = skycoord_to_pixel(source_center_SgrA, on.wcs)
    large_gaus.ampl = 1
    freeze(large_gaus.ampl)
    freeze(large_gaus.xpos)
    freeze(large_gaus.ypos)

    #Modele symetric central gauss centre sur SgrA
    central_gauss = Gauss2D("central_gauss")
    central_gauss.xpos, central_gauss.ypos = skycoord_to_pixel(source_center_SgrA, on.wcs)
    freeze(central_gauss.xpos)
    freeze(central_gauss.ypos)

    #Arc_source
    arc_source=SkyCoord(0.130, -0.139, unit='deg', frame="galactic")
    mygaus_arcsource = Gauss2D("Arc Source")
    mygaus_arcsource.xpos, mygaus_arcsource.ypos = skycoord_to_pixel(arc_source, on.wcs)
    freeze(mygaus_arcsource.xpos)
    freeze(mygaus_arcsource.ypos)
    mygaus_arcsource.fwhm =1
    freeze(mygaus_arcsource.fwhm)



    #Define the region of interest and the region to ignore for the fit
    #ignore region in a box that mask J1734-303
    source_J1745_303 = SkyCoord(358.76, -0.6, unit='deg', frame="galactic")
    source_J1745_303_xpix, source_J1745_303_ypix = skycoord_to_pixel(source_J1745_303, on.wcs)
    width=100
    height=80
    name_region = "box(" + str(source_J1745_303_xpix+20) + "," + str(source_J1745_303_ypix-20) + "," + str(width) + "," + str(height) +")"
    ignore2d(name_region)

    #notice a region of interest for latitude inferior to 0.8 and longitude inferior to 2
    pix_deg = on.to_image_hdu().header["CDELT2"]
    lat=1.6/ pix_deg#Pour aller a plus et -0.8 as did Anne
    lon=4 / pix_deg#Pour aller a plus ou moins 2deg as did Anne
    x_pix_SgrA=skycoord_to_pixel(source_center_SgrA, on.wcs)[0]
    y_pix_SgrA=skycoord_to_pixel(source_center_SgrA, on.wcs)[1]
    name_interest = "box(" + str(x_pix_SgrA) + "," + str(y_pix_SgrA) + "," + str(lon) + "," + str(lat) +")"
    notice2d(name_interest)


    #set stat and method for the fit
    set_stat("cstat")
    set_method("neldermead")

    #list_src = [psf_SgrA(mygaus_SgrA),psf_G0p9(mygaus_G0p9),psf_SgrA(asym_gaus),psf_SgrA(large_gaus * CS),
     #           psf_SgrA(central_gauss), psf_SgrA(mygaus_arcsource)]
    list_src = [psf_SgrA(mygaus_SgrA), psf_G0p9(mygaus_G0p9)]
    model = bkg

    with PdfPages("profiles_" + str("%.2f"%E1) + "_" + str("%.2f"%E2) + "_TeV.pdf") as pdf:
        for i_src,src in enumerate(list_src):
            model += src
            set_full_model(model)
            fit()
            result=get_fit_results()
            try:
                table_models = join(table_models.filled(-1000), result_table(result, int(i_src)), join_type='outer')

            except NameError:
                table_models= result_table(result, int(i_src))

            covar()
            covar_res= get_covar_results()

            try:
                table_covar = join(table_covar.filled(-1000), covar_table(covar_res), join_type='outer')

            except NameError:
                table_covar = covar_table(covar_res)
            save_resid("residual_step_"+str(i_src)+"_" + str("%.2f"%E1) + "_" + str("%.2f"%E2) + "_TeV.fits", clobber=True)


            #Profil lattitude et longitude
            shape = np.shape(on.data)
            mask = get_data().mask.reshape(shape)
            map_data = get_data().y.reshape(shape) * mask
            model_map = get_model()(get_data().x0, get_data().x1).reshape(shape) * mask
            resid = map_data - model_map
            coord = on.coordinates()

            #Longitude profile
            i_b = np.where((coord.b[:, 0] < on.center().b + 0.25 * u.deg) & (coord.b[:, 0] > on.center().b - 0.25 * u.deg))[0]
            npix_l=np.sum(np.flipud(mask[i_b,:]), axis=0)
            profile_l_model = np.sum(model_map[i_b, :], axis=0)/npix_l
            profile_l_on = np.sum(map_data[i_b, :], axis=0)/npix_l
            #Ca donne des coups par arcmin2 car on prend en compte qu on ne cumula pas le meme nombre de pixel pour chaque
            # longitude vu qu il y a des regions d exclusions
            profile_l_resid = np.sum(resid[i_b, :], axis=0)/npix_l
            err_l = np.sqrt(profile_l_on/npix_l)
            l = coord.l[0, :]
            l.value[np.where(l > 180 * u.deg)] = l.value[np.where(l > 180 * u.deg)] - 360
            resid_l_rebin, l_rebin, err_l_rebin = rebin_profile(profile_l_resid, l, err_l, nrebin=3)

            #Latitude profile
            l_center = on.center().l
            if l_center > 180 * u.deg:
                l_center = l_center - 360 * u.deg
            i_l = np.where((l < l_center + 1.5 * u.deg) & (l > l_center - 1.5 * u.deg))[0]
            npix_b=np.sum(np.flipud(mask[:,i_l]), axis=1)
            profile_b_model = np.sum(model_map[:, i_l], axis=1)/npix_b
            profile_b_on = np.sum(map_data[:, i_l], axis=1)/npix_b
            profile_b_resid = np.sum(resid[:, i_l], axis=1)/npix_b
            err_b = np.sqrt(profile_b_on/npix_b)
            resid_b_rebin, b_rebin, err_b_rebin = rebin_profile(profile_b_resid, coord.b[:, 0], err_b, nrebin=3)


            fig = pt.figure()
            ax = fig.add_subplot(2, 1, 1)
            pt.plot(l.value, profile_l_model, label="model")
            pt.plot(l.value, profile_l_on, label="on data")
            pt.xlim(-1.5, 1.5)
            pt.gca().invert_xaxis()
            pt.legend()
            ax = fig.add_subplot(2, 1, 2)
            pt.errorbar(l_rebin.value, resid_l_rebin, yerr=err_l_rebin, linestyle='None', marker="o", label= "Step= "+str(i_src))
            pt.axhline(y=0, color='red', linewidth=2)
            pt.legend()
            pt.ylabel("residual")
            pt.xlabel("longitude (degrees)")
            pt.title("longitude profile")
            pt.xlim(-1.5, 1.5)
            pt.gca().invert_xaxis()
            pdf.savefig()

            fig =pt.figure()
            ax = fig.add_subplot(2, 1, 1)
            pt.plot(coord.b[:, 0].value, profile_b_model, label="model")
            pt.plot(coord.b[:, 0].value, profile_b_on, label="on data")
            pt.xlim(-1, 1)
            pt.legend()
            ax = fig.add_subplot(2, 1, 2)
            pt.errorbar(b_rebin.value, resid_b_rebin, yerr=err_b_rebin, linestyle='None', marker="o", label= "Step= "+str(i_src))
            pt.axhline(y=0, color='red', linewidth=2)
            pt.legend()
            pt.ylabel("residual")
            pt.xlabel("latitude (degrees)")
            pt.title("latitude profile")
            pt.xlim(-1, 1)
            pdf.savefig()


filename_table_result="fit_result_"+ str("%.2f"%E1) + "_" + str("%.2f"%E2) + "_TeV.txt"
table_models.write(filename_table_result, format="ascii")
filename_covar_result="covar_result_"+ str("%.2f"%E1) + "_" + str("%.2f"%E2) + "_TeV.txt"
table_covar.write(filename_covar_result, format="ascii")
