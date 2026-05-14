import numpy as np
from astropy.io import fits
from astropy.convolution import convolve_fft
import astropy.units as u
import astropy.constants as const
from astropy.table import QTable
from PIL import Image
from scipy.interpolate import interp1d
import os
from skimage.transform import rescale
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib_scalebar.dimension import _Dimension
from photutils.segmentation import make_2dgaussian_kernel
from typing import Union
import numba
import time
from .utils import read_properties, get_wavelength_unit, extend, convert_to_rgb, setup_logging, read_config, galaxygenius_data_dir

try:
    import pyfftw
    import scipy.fft as sfft
    sfft.set_backend(pyfftw.interfaces.scipy_fft)
except:
    pass

class ParsecDimension(_Dimension):
    def __init__(self):
        super().__init__('pc')
        self.add_units('kpc', 1000)
    
class AngleDimension(_Dimension):
    def __init__(self):
        super().__init__(r'${^{\prime\prime}}$')
        self.add_units(r'${^{\prime}}$', 60)

class PostProcess:

    def __init__(self, subhaloID: int):

        self.subhaloID = subhaloID
        self.dataCubeDir = f'dataCubes/Subhalo_{self.subhaloID}'
        self.properties = read_properties(self.dataCubeDir)
        self.logger = setup_logging(os.path.join(os.getcwd(), 'galaxyGenius.log'))
        self.logger.info(f'Initializing PostProcess class.')
        self.config = read_config(self.dataCubeDir)
        self.dataDir = galaxygenius_data_dir()

    def __survey_enlarge_ratio(self, survey: str) -> float:
        """Angular zoom for display-only figures; missing or invalid → 1 (legacy behavior)."""
        try:
            r = float(self.config.get(f'enlargeRatio_{survey}', 1.0))
        except (TypeError, ValueError):
            return 1.0
        if not np.isfinite(r) or r <= 0:
            return 1.0
        return max(1.0, r)

    def __center_enlarge_for_display(self, arr: np.ndarray, ratio: float) -> np.ndarray:
        """Center-crop to ``1/ratio`` of the linear extent (``ratio``>1); native pixels only, no resampling."""
        if ratio <= 1.0 or not np.isfinite(ratio):
            return arr
        a = np.asarray(arr)
        if a.ndim == 2:
            h, w = a.shape
        elif a.ndim == 3 and a.shape[2] == 3:
            h, w = a.shape[0], a.shape[1]
        else:
            return arr
        ch = max(1, int(round(h / ratio)))
        cw = max(1, int(round(w / ratio)))
        ch = min(ch, h)
        cw = min(cw, w)
        y0 = (h - ch) // 2
        x0 = (w - cw) // 2
        if a.ndim == 2:
            return np.asarray(a[y0 : y0 + ch, x0 : x0 + cw])
        return np.asarray(a[y0 : y0 + ch, x0 : x0 + cw, :])

    @staticmethod
    def __scalebar_fixed_display(value: float) -> float:
        """Round scale bar numeric label to two decimals (compact, no float noise)."""
        if not np.isfinite(value):
            return 0.0
        return float(f'{float(value):.2f}')

    def __precompile_numba(self):
        dummy_img = np.ones((10, 10, 10), dtype=np.float32)
        dummy_trans = np.ones(10, dtype=np.float32)
        dummy_wave = np.ones(10, dtype=np.float32)
        self.integrate_bandpass(dummy_img, dummy_trans, dummy_wave)
        
    # def __read_configs(self):
        
    #     with open(self.dataCubeDir + '/config.json', 'r') as file:
    #         config = json.load(file)
            
    #     return config
        
    # def __get_properties(self) -> dict:
    #     path = os.path.join(self.dataCubeDir, 'properties.json')
    #     with open(path, 'r') as file:
    #         properties = json.load(file)
        
    #     return properties
            
    def __load_method(self, format: str) -> callable:
        if format == '.npy':
            method = np.load
        else:
            method = np.loadtxt

        return method

    def __get_throughputs(self, survey: str) -> list:

        filters = self.config[f'filters_{survey}']
        filterDir = os.path.join(self.dataDir, f'filters/{survey}')
        filterLs = os.listdir(filterDir)
        filterNames = [name.split('.')[0] for name in filterLs]

        throughputs = []
        for fil in filters:
            filename = filterLs[filterNames.index(fil)]
            throughput_file = os.path.join(filterDir, filename)
            wavelength_unit = get_wavelength_unit(throughput_file)
            throughput = np.loadtxt(throughput_file).astype(np.float32)
            thr_table = QTable([throughput[:, 0] * wavelength_unit, 
                           throughput[:, 1] * u.dimensionless_unscaled],
                               names=['wavelength', 'throughput'])
            # throughput[:, 0] = throughput[:, 0] * wavelength_unit
            throughputs.append(thr_table)

        return throughputs

    def __get_PSFs(self, survey: str) -> list:

        if self.config[f'PSFFromFile_{survey}']:
            filters = self.config[f'filters_{survey}']
            psfDir = os.path.join(self.dataDir, f'PSFs/{survey}')
            psfLs = os.listdir(psfDir)

            psfNames = [name.split('.')[0] for name in psfLs]

            PSFs = []
            for fil in filters:
                filename = psfLs[psfNames.index(fil)]
                psf_file = os.path.join(psfDir, filename)
                
                psf = self.__load_method(os.path.splitext(psf_file)[1])(psf_file)
                if psf.shape[0] % 2 == 0: # make sure the size is odd
                    psf = np.pad(psf, ((0, 1), (0, 1)), 'constant')
                psf = psf / np.sum(psf)
                PSFs.append(psf)
        
        else:
            numfilters = len(self.config[f'filters_{survey}'])
            fwhm_in_arcsec = self.config[f'PSFFWHM_{survey}']
            fwhm_in_arcsec = extend(fwhm_in_arcsec, numfilters) # with arcsec
            
            pixelScales = self.properties[f'angleRes_{survey}'] # with arcsec

            stds = [(fwhm / ps).value for fwhm, ps in zip(fwhm_in_arcsec, pixelScales)]
            
            PSFs = []
            for std in stds:
                size = std * 11
                size = int(size)
                if size % 2 == 0: # make sure the size is odd
                    size += 1
                kernel = make_2dgaussian_kernel(std, size)
                PSFs.append(kernel)
            
        return PSFs

    # def __calculate_bandpass(self, img: np.ndarray, tran: np.ndarray, 
    #                           wave: np.ndarray, factor: float) -> np.ndarray:
    #     bandpass_img = factor * trapezoid(img * tran.reshape(-1, 1, 1) \
    #                                       * wave.reshape(-1, 1, 1), wave, axis=0)
    #     return bandpass_img
    
    # @staticmethod
    # @numba.njit(
    # "float32[:, :](float32[:, :, :], float32[:], float32[:], float32)", 
    # parallel=True, cache=True, fastmath=True)
    # def integrate_bandpass(img, tran, wave, factor):
    #     n = len(wave)
    #     h, w = img.shape[1], img.shape[2]
    #     out = np.zeros((h, w), dtype=np.float32)
    #     for i in numba.prange(h):
    #         for j in range(w):
    #             integral = 0.0
    #             for k in range(1, n):
    #                 y1 = img[k-1, i, j] * tran[k-1] * wave[k-1]
    #                 y2 = img[k, i, j] * tran[k] * wave[k]
    #                 dx = wave[k] - wave[k-1]
    #                 integral += (y1 + y2) / 2.0 * dx
    #             out[i, j] = factor * integral
    #     return out
    
    @staticmethod
    @numba.njit(parallel=True, cache=True, fastmath=True)
    def integrate_bandpass(img, tran, wave):
        n = len(wave)
        h, w = img.shape[1], img.shape[2]
        out = np.zeros((h, w), dtype=np.float32)
        for i in numba.prange(h):
            for j in range(w):
                integral = 0.0
                for k in range(1, n):
                    y1 = img[k-1, i, j] * tran[k-1] * wave[k-1]
                    y2 = img[k, i, j] * tran[k] * wave[k]
                    dx = wave[k] - wave[k-1]
                    integral += (y1 + y2) / 2.0 * dx
                out[i, j] = integral
        return out
    
    def __bandpass_images(self, dataCube: str, survey: str, throughputs: list,
                          PSFs: Union[list, None]=None, bkgNoise: Union[dict, None]=None) -> list:
        
        file = fits.open(dataCube)
        wave_sed = (file[1].data['grid_points']).astype(np.float32) * u.um
        numExp = self.properties[f'numExposure_{survey}']
        exposureTime = self.properties[f'exposureTime_{survey}'] # with sec
        areaMirror = np.pi * (self.properties[f'aperture_{survey}'] / 2)**2 # with m**2

        ps_in_sr = self.properties[f'ps_in_sr_{survey}'] # with sr
        interpolate_ratios = self.properties[f'ratios_{survey}'] # list
        pixelScales = self.properties[f'angleRes_{survey}'] # with arcsec
        
        filters = self.properties[f'filters_{survey}']

        data = file[0].data.astype(np.float32) * u.MJy / u.sr
        
        conversion_to_Jy = []
        gains = []
        pivots = []
        bandpass_images = []
                
        for i, thr in enumerate(throughputs):
            
            self.logger.info(f'Processing {survey} {filters[i]}.')
            
            wave_min, wave_max = np.min(thr['wavelength']), np.max(thr['wavelength'])
            idx = np.where((wave_sed > wave_min) & (wave_sed < wave_max))[0]
            
            # throughput -- wave (angstrom)
            interp = interp1d(thr['wavelength'].to(u.angstrom).value, thr['throughput'])
            wave_in = wave_sed[idx].to(u.angstrom).value.astype(np.float32)
            trans_in = interp(wave_in).astype(np.float32)
            
            # data is in f_nu (MJy/sr), converted to erg cm**-2 s**-1 angstrom**-1 sr**-1
            f_lam = (data[idx] * const.c / (wave_in.reshape(-1, 1, 1) * u.angstrom)**2)
            f_lam = f_lam.to(u.erg / u.cm**2 / u.s / u.angstrom / u.sr)
            f_lam_unit = f_lam.unit
            
            const_factor = exposureTime[i] * numExp[i] * areaMirror / (const.h * const.c)
            integral = self.integrate_bandpass(
                f_lam.value, trans_in, wave_in
            ) * f_lam_unit * u.angstrom**2
            image_in_count = const_factor * integral
            image_in_count = image_in_count.to(u.sr**-1)
            
            gain = const.h / (areaMirror * trapezoid(trans_in / wave_in, wave_in))
            
            pivot_numerator = trapezoid(trans_in, wave_in) * u.angstrom
            pivot_denominator = trapezoid(trans_in * wave_in ** -2, wave_in) * u.angstrom**-1
            pivot = np.sqrt(pivot_numerator / pivot_denominator)
            
            Jy_converter = trapezoid(trans_in * wave_in, wave_in) * u.angstrom**2
            Jy_converter = Jy_converter / (const.h * const.c) * areaMirror\
                            * numExp[i] * exposureTime[i]
            Jy_converter = pivot**2 / const.c / Jy_converter
            Jy_converter = Jy_converter.to(u.Jy)
            
            conversion_to_Jy.append(Jy_converter.astype(np.float32))
            gains.append(gain.astype(np.float32))
            pivots.append(pivot.astype(np.float32))
            bandpass_images.append(image_in_count)
            
        # for i, thr in enumerate(throughputs):
            
        #     wave_min, wave_max = np.min(thr['wavelength']), np.max(thr['wavelength'])
        #     idx = np.where((wave_sed > wave_min) & (wave_sed < wave_max))[0]
            
        #     interp = interp1d(thr['wavelength'].to(u.angstrom).value, thr['throughput'])
        #     wave_in = wave_sed[idx].to(u.angstrom).value
        #     trans_in = interp(wave_in)
            
        #     if data_norm is None or data_norm.shape != (len(idx), data.shape[1], data.shape[2]):
        #         if data_norm is not None:
        #             del data_norm
        #             gc.collect()
        #         data_norm = np.empty((len(idx), data.shape[1], data.shape[2]), dtype=np.float32)
            
        #     np.divide(data[idx], wave_in.reshape(-1, 1, 1)**2, out=data_norm)
        #     factor = (u.MJy/u.sr) * const.c / u.angstrom**2 \
        #                 * numExp[i] * exposureTime[i] \
        #                     * areaMirror / (const.c * const.h) * u.angstrom**2
        #     factor = factor.to(u.sr**-1)
        #     factor = factor.value.astype(np.float32)
            
        #     gain = const.h / (areaMirror * trapezoid(trans_in / wave_in, wave_in))
        #     numerator = trapezoid(trans_in, wave_in)
        #     denominator = trapezoid(trans_in * wave_in ** -2, wave_in)
        #     pivot = np.sqrt(numerator / denominator)
            
        #     Jy_converter = trapezoid(trans_in * wave_in, wave_in) * u.angstrom**2
        #     Jy_converter = Jy_converter / (const.h * const.c) * areaMirror\
        #                     * numExp[i] * exposureTime[i]
        #     Jy_converter = (pivot * u.angstrom)**2 / const.c / Jy_converter
        #     Jy_converter = Jy_converter.to(u.Jy).value
            
        #     conversion_to_Jy.append(Jy_converter.astype(np.float32))
        #     gains.append(gain.astype(np.float32))
        #     pivots.append(pivot.astype(np.float32))
            
        #     bandpass_images.append(self.integrate_bandpass(data_norm, trans_in, wave_in, factor))
        
        resized_imgs = []
        for i, ratio in enumerate(interpolate_ratios):
            resized_imgs.append(rescale(bandpass_images[i], ratio))

        bandpass_images = resized_imgs

        converted_imgs = []
        for i, ps_sr in enumerate(ps_in_sr):
            # remove sr, and remove unit
            converted_imgs.append((bandpass_images[i] * ps_sr).value)
        
        bandpass_images = converted_imgs

        # apply PSF effects
        if PSFs is not None: 
            images_with_psf = []
            for i, img in enumerate(bandpass_images):
                # print(img.shape, PSFs[i].shape)
                images_with_psf.append(convolve_fft(img, PSFs[i], normalize_kernel=True))
            
            bandpass_images = images_with_psf

        # add instrumental noise
        if bkgNoise is not None:
            
            if bkgNoise['noiseType'] == 'instrument':
                
                # skyBkg = np.array(bkgNoise['skyBkg'])
                # darkCurrent = np.array(bkgNoise['darkCurrent'])
                # readOut = np.array(bkgNoise['readOut'])
                # exposureTime = np.array(bkgNoise['exposureTime'])
                # numExposure = np.array(bkgNoise['numExposure'])

                skyBkg = bkgNoise['skyBkg'] # with s**-1
                darkCurrent = bkgNoise['darkCurrent'] # with s**-1
                readOut = bkgNoise['readOut'] # dimless
                exposureTime = bkgNoise['exposureTime'] # with s
                numExposure = bkgNoise['numExposure'] # dimless
                
                images_with_bkg_in_electron = []
                for i, img in enumerate(bandpass_images):
                    mean = ((skyBkg[i] + darkCurrent[i]) * exposureTime[i] * numExposure[i]).value
                    img = img + mean
                    img = np.random.poisson(img).astype(np.int32)
                    for _ in range(int(numExposure[i])):
                        readnoise = np.random.normal(loc=0, scale=readOut[i].value, size=img.shape)
                        readnoise = np.around(readnoise, 0).astype(np.int32)
                        img = img + readnoise
                    img = img.astype(np.float32) - mean.astype(np.float32)
                    images_with_bkg_in_electron.append(img)
                    
                if self.config['imageUnit'] == 'electron':
                    images_with_bkg = images_with_bkg_in_electron
                elif self.config['imageUnit'] == 'flux':
                    images_with_bkg = []
                    for i, img in enumerate(images_with_bkg_in_electron):
                        img = (img * conversion_to_Jy[i]) # with Jy
                        images_with_bkg.append(img.value)
            
            elif bkgNoise['noiseType'] == 'limitingMagnitude':
                
                limitMag = bkgNoise['limitMag'] # dimless
                limitSNR = bkgNoise['limitSNR'] # dimless
                limitAperture = bkgNoise['limitAperture'] # arcsec
                zeroPoint = bkgNoise['zeroPoint'] # dimless
                exposureTime = bkgNoise['exposureTime'] # with s
                numExposure = bkgNoise['numExposure'] # dimless
                
                npix = (np.pi * (limitAperture / 2)**2) / (pixelScales**2)
                
                images_with_bkg_in_Jy = []
                for i, img in enumerate(bandpass_images):
                    img_in_Jy = img * conversion_to_Jy[i]
                    first_term = (gains[i] * img_in_Jy) / (numExposure[i] * exposureTime[i])
                    second_term = ((1 / limitSNR[i]) * 10**((zeroPoint[i] - limitMag[i]) / 2.5) * u.Jy)**2 * (1 / npix[i])
                    third_term = gains[i] / (numExposure[i] * exposureTime[i] * npix[i]) * 10**((zeroPoint[i] - limitMag[i]) / 2.5) * u.Jy
                    
                    first_term = first_term.to(u.Jy**2)
                    second_term = second_term.to(u.Jy**2)
                    third_term = third_term.to(u.Jy**2)
                    
                    total = (first_term + second_term - third_term).value
                    
                    noise = np.sqrt(total) # in Jy, without unit
                    noises = np.random.normal(loc=0, scale=noise, size=img.shape)
                    images_with_bkg_in_Jy.append(img + noises)
                    
                if self.config['imageUnit'] == 'electron':
                    images_with_bkg = []
                    for i, img in enumerate(images_with_bkg_in_Jy):
                        img = img / conversion_to_Jy[i]
                        images_with_bkg.append(img.value)
                        
                elif self.config['imageUnit'] == 'flux':
                    images_with_bkg = images_with_bkg_in_Jy
                    
            bandpass_images = images_with_bkg
        
        return bandpass_images
            
    def __saveBandpassImages(self, images: list, survey: str):

        # images are in [numViews, [numfilters, [numPixels, numPixels]]]

        hdulist = fits.HDUList()
        primary_hdu = fits.PrimaryHDU()
        hdulist.append(primary_hdu)
        
        numfilters = self.properties[f'numfilters_{survey}']
        
        unit_dict = {'electron': 'e', 'flux': 'Jy'}
        unit_comment_dict = {
            'electron': 'Unit of image, in electron counts',
            'flux': 'Unit of image, in Jy',
        }
        unit_type = self.config['imageUnit']
        imageUnit = unit_dict[unit_type]
        imageUnitComment = unit_comment_dict[unit_type]
        
        for i in range(numfilters):
            header = fits.Header()
            header['SNAPNUM'] = (self.config['snapNum'], 'Snapshot ID')
            header['ID'] = (self.properties['subhaloID'], 'Subhalo ID')
            header['MASS'] = (np.log10(self.properties['stellarMass'].value), 'Subhalo stellar mass, in log10 scale (Msun)')
            header['SURVEY'] = (survey, 'Survey')
            header['NFILTERS'] = (numfilters, 'Number of filters')
            header['NUMVIEWS'] = (self.properties['numViews'], 'Number of views')
            for count in range(self.properties['numViews']):
                header[f'INCLI_{count:02d}'] = (self.properties['inclinations'][count], 
                                                f'Inclination angle, in deg for view {count:02d}')
            for count in range(self.properties['numViews']):
                header[f'AZIMU_{count:02d}'] = (self.properties['azimuths'][count], 
                                                  f'Azimuth angle, in deg for view {count:02d}')
            header['FILTER'] = (self.properties[f'filters_{survey}'][i], 'Filter')
            header['EXPTIME'] = (self.properties[f'exposureTime_{survey}'][i].value, 'Exposure time, in s')
            header['EXPNUM'] = (int(self.properties[f'numExposure_{survey}'][i]), 'Number of exposures')
            header['APERTURE'] = (self.properties[f'aperture_{survey}'].value, 'Aperture size, in meter')
            header['UNIT'] = (imageUnit, imageUnitComment)
            if self.properties['viewRedshift'] is not None:
                header['REDSHIFT'] = (self.properties['viewRedshift'], 'Viewing redshift')
            header['COSMO'] = (self.properties['cosmology'], 'Cosmology')
            header['BOXSIZE'] = (self.properties['boxLength'].to(u.pc).value, 'Box size, in pc')
            header['FOV'] = (self.properties['fieldOfView'].value, 'Field of view, in arcsec')
            header['PHYFOV'] = (self.properties['fovSize'].value, 'Physical field of view, in pc')
            header['LUMIDIS'] = (self.properties['lumiDis'].value, 'Luminosity distance, in Mpc')
            header['RESOL'] = (self.properties[f'resolution_{survey}'][i].value, 'Pixel scale, in pc')
            header['PS'] = (self.properties[f'angleRes_{survey}'][i].value, 'Pixel scale, in arcsec')
            
            if imageUnit == 'electron':
                header['E2JY'] = (self.conversion_to_Jy[i], 'Conversion factor from electron to Jy')
            
            if imageUnit == 'flux':
                header['JY2E'] = (1 / self.conversion_to_Jy[i], 'Conversion factor from Jy to electron')

            imgs_in_view = [images[nv][i] for nv in range(self.properties['numViews'])]
            imgs_in_view = np.stack(imgs_in_view, axis=0)

            hdu = fits.ImageHDU(data=imgs_in_view, header=header)
            hdulist.append(hdu)
            
        savedImageName = f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_images.fits'
        hdulist.writeto(savedImageName, overwrite=True)

    def __saveSEDs(self, sedFilenames: list, survey: str):
        hdulist = fits.HDUList()
        primary_hdu = fits.PrimaryHDU()
        hdulist.append(primary_hdu)
        
        for i in range(self.properties['numViews']):
            sed = np.loadtxt(sedFilenames[i])
            
            # Create column definitions for the table
            col1 = fits.Column(name='WAVELENGTH', format='E', unit='micron', array=sed[:, 0])
            col2 = fits.Column(name='FLUX', format='E', unit='Jy', array=sed[:, 1])
            
            # Create the table HDU
            cols = fits.ColDefs([col1, col2])
            tbhdu = fits.BinTableHDU.from_columns(cols)
            
            # Add header information
            header = tbhdu.header
            header['WUNIT'] = ('micron', 'Units of wavelength')
            header['FUNIT'] = ('Jy', 'Units of flux in F_nu')
            header['INCLI'] = (self.properties['inclinations'][i], 'Inclination angle, in deg')
            header['AZIMUTH'] = (self.properties['azimuths'][i], 'Azimuth angle, in deg')
            if self.properties['viewRedshift'] is not None:
                header['REDSHIFT'] = (self.properties['viewRedshift'], 'Viewing redshift')
            header['lumiDis'] = (self.properties['lumiDis'].value, 'luminosity distance, in Mpc')
            header['FOV'] = (self.properties['fieldOfView'].value, 'Field of view, in arcsec')
            header['PHYFOV'] = (self.properties['fovSize'].value, 'Physical field of view, in pc')
            header['VIEW'] = (i, 'View index')
            
            hdulist.append(tbhdu)
            
        savedSEDName = f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_SEDs.fits'
        hdulist.writeto(savedSEDName, overwrite=True)

    def __plot_image(self, image: np.ndarray, pixelscale: float, resolution: float,
                     savedFilename: Union[str, None]=None, enlarge_ratio: float = 1.0):

        arr = self.__center_enlarge_for_display(np.asarray(image), enlarge_ratio)
        # Cropped pixels keep the same on-sky / physical scale as the simulation grid.
        eff_pixelscale = pixelscale
        eff_resolution = resolution
        numPixels = arr.shape[0]
        is_rgb = arr.ndim == 3 and arr.shape[2] == 3

        if self.properties['viewRedshift'] is not None:
            z = str(np.around(self.properties['viewRedshift'], 2))
        else:
            z = '0.00'

        logM = str(np.around(np.log10(self.properties['stellarMass'].value), 1))

        if is_rgb:
            fig, ax = plt.subplots()
            ax.axis('off')
            im = ax.imshow(arr, origin='lower')
        else:
            show_img = arr.astype(float).copy()
            show_img[show_img < 0] = 0.0
            show_img = show_img + 1e-5
            log_data = np.log10(show_img)
            vmin = float(np.min(log_data))
            vmax = float(np.max(log_data))
            if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
                vmax = vmin + 1e-6
            image_unit = self.config['imageUnit']
            if image_unit == 'electron':
                cbar_label = r'$\log_{10}\,e$'
            else:
                cbar_label = r'$\log_{10}\,F_\nu\ [{\rm Jy}]$'
            fig = plt.figure(figsize=(7.0, 5.5))
            gs = GridSpec(1, 2, figure=fig, width_ratios=[1.0, 0.06], wspace=0.2)
            ax = fig.add_subplot(gs[0, 0])
            cax = fig.add_subplot(gs[0, 1])
            ax.set_axis_off()
            im = ax.imshow(
                log_data,
                origin='lower',
                cmap='gray',
                vmin=vmin,
                vmax=vmax,
            )
            cbar_text_color = 'black'
            cb = fig.colorbar(im, cax=cax)
            cb.set_label(cbar_label, fontsize=12, color=cbar_text_color)
            cb.ax.tick_params(labelsize=10, colors=cbar_text_color)
            cb.ax.yaxis.label.set_color(cbar_text_color)

        scalebarSize = 0.25 * numPixels * eff_pixelscale
        if scalebarSize > 60:
            scalebarUnit = '${^{\prime}}$'
            scalebarSize = self.__scalebar_fixed_display(scalebarSize / 60)
            ps = eff_pixelscale / 60
        else:
            scalebarUnit = '${^{\prime\prime}}$'
            scalebarSize = self.__scalebar_fixed_display(scalebarSize)
            ps = eff_pixelscale
            
        angle_dim = AngleDimension()
        
        scalebar = ScaleBar(ps, scalebarUnit, dimension=angle_dim, 
                            fixed_value=scalebarSize, fixed_units=scalebarUnit, frameon=False,
                            location='lower right', scale_loc='top',
                            color='white', font_properties={'size': 16})
        
        
        if scalebarSize < 1: # 1 in arcsec
            # fall back to show physical scale when pixel scale is too small
            scalebarSize = 0.25 * numPixels * eff_resolution
            
            if scalebarSize > 10**3:
                scalebarUnit = 'kpc'
                scalebarSize = self.__scalebar_fixed_display(scalebarSize / 10**3)
                res = eff_resolution / 10**3
            else:
                scalebarUnit = 'pc'
                scalebarSize = self.__scalebar_fixed_display(scalebarSize)
                res = eff_resolution
            
            parsec_dim = ParsecDimension()
            scalebar = ScaleBar(res, scalebarUnit, dimension=parsec_dim,
                                fixed_value=scalebarSize, fixed_units=scalebarUnit,
                                frameon=False, location='lower right', scale_loc='top',
                                color='white', font_properties={'size': 16})
        
        ax.add_artist(scalebar)
        ax.text(x=0.05, y=0.15, s=fr'${{\rm log}}M_{{\star}} = {logM}$',
                         fontsize=12, transform=ax.transAxes, color='white')
        ax.text(x=0.05, y=0.1, s=fr'$z$ = {z}', fontsize=12,
                transform=ax.transAxes, color='white')
        ax.text(x=0.05, y=0.05, s=f'ID:{self.subhaloID}', fontsize=12,
                transform=ax.transAxes, color='white')
        plt.tight_layout()
        if savedFilename is not None:
            plt.savefig(savedFilename)
            plt.close()
        else:
            plt.show()

    def __save_all_bands_montage(self, images: list, survey: str, filters: list) -> None:
        """
        One row per view: all filters side by side with log-scaled flux/e⁻ counts,
        shared color scale (gray), a colorbar, per-panel angular scalebars,
        ``{survey}.{band}`` labels (top-right, inset in points), and survey metadata
        at the lower-left of the first panel only. Figure and panels use a white
        background; band labels, metadata, and the scalebar on the first panel are
        white. The colorbar label and ticks use black. Angular scalebar is drawn
        only on the first band panel.
        """
        save_dir = f'mock_{survey}/Subhalo_{self.subhaloID}'
        os.makedirs(save_dir, exist_ok=True)
        num_views = self.properties['numViews']
        nband = len(filters)
        if nband == 0:
            return

        if self.properties['viewRedshift'] is not None:
            z_str = str(np.around(self.properties['viewRedshift'], 2))
        else:
            z_str = '0.00'
        logM = str(np.around(np.log10(self.properties['stellarMass'].value), 1))

        image_unit = self.config['imageUnit']
        if image_unit == 'electron':
            cbar_label = r'$\log_{10}\,e$'
        else:
            cbar_label = r'$\log_{10}\,F_\nu\ [{\rm Jy}]$'

        enlarge_ratio = self.__survey_enlarge_ratio(survey)

        def add_scalebar(ax, num_pixels: int, pixelscale: float, resolution: float,
                         bar_color: str):
            scalebar_size = 0.25 * num_pixels * pixelscale
            if scalebar_size > 60:
                scalebar_unit = r'${^{\prime}}$'
                scalebar_size = self.__scalebar_fixed_display(scalebar_size / 60)
                ps = pixelscale / 60
            else:
                scalebar_unit = r'${^{\prime\prime}}$'
                scalebar_size = self.__scalebar_fixed_display(scalebar_size)
                ps = pixelscale
            angle_dim = AngleDimension()
            scalebar = ScaleBar(
                ps,
                scalebar_unit,
                dimension=angle_dim,
                fixed_value=scalebar_size,
                fixed_units=scalebar_unit,
                frameon=False,
                location='lower right',
                scale_loc='top',
                color=bar_color,
                font_properties={'size': 16},
            )
            if scalebar_size < 1:
                scalebar_size = 0.25 * num_pixels * resolution
                if scalebar_size > 10**3:
                    scalebar_unit = 'kpc'
                    scalebar_size = self.__scalebar_fixed_display(scalebar_size / 10**3)
                    res = resolution / 10**3
                else:
                    scalebar_unit = 'pc'
                    scalebar_size = self.__scalebar_fixed_display(scalebar_size)
                    res = resolution
                parsec_dim = ParsecDimension()
                scalebar = ScaleBar(
                    res,
                    scalebar_unit,
                    dimension=parsec_dim,
                    fixed_value=scalebar_size,
                    fixed_units=scalebar_unit,
                    frameon=False,
                    location='lower right',
                    scale_loc='top',
                    color=bar_color,
                    font_properties={'size': 16},
                )
            ax.add_artist(scalebar)

        overlay_color = 'white'
        band_labels = [f'{survey}.{f}' for f in filters]
        band_font = 20

        for view_i in range(num_views):
            band_stack = images[view_i]
            show_imgs = []
            mins = []
            maxs = []
            for j in range(nband):
                show_img = np.asarray(band_stack[j], dtype=float).copy()
                show_img[show_img < 0] = 0.0
                show_img = show_img + 1e-5
                show_img = self.__center_enlarge_for_display(show_img, enlarge_ratio)
                show_imgs.append(show_img)
                mins.append(np.min(show_img))
                maxs.append(np.max(show_img))
            vmin = np.log10(np.min(mins))
            vmax = np.log10(np.max(maxs))

            total_width = nband * 4 + 1
            height = 4
            fig = plt.figure(figsize=(total_width, height))
            fig.patch.set_facecolor('white')
            gs = GridSpec(height, total_width, figure=fig)
            im_last = None
            for j in range(nband):
                ax = fig.add_subplot(gs[:, 4 * j : 4 * (j + 1)])
                ax.set_axis_off()
                ax.set_facecolor('white')
                im_last = ax.imshow(
                    np.log10(show_imgs[j]),
                    vmin=vmin,
                    vmax=vmax,
                    cmap='gray',
                    origin='lower',
                    aspect='equal',
                )
                if j == 0:
                    pixelscale = self.properties[f'angleRes_{survey}'][j].value
                    resolution = self.properties[f'resolution_{survey}'][j].value
                    add_scalebar(ax, show_imgs[j].shape[0], pixelscale, resolution,
                                 bar_color=overlay_color)
                ax.annotate(
                    band_labels[j],
                    xy=(1.0, 1.0),
                    xycoords='axes fraction',
                    xytext=(-8, -8),
                    textcoords='offset points',
                    ha='right',
                    va='top',
                    fontsize=band_font,
                    fontstyle='italic',
                    color=overlay_color,
                )
                if j == 0:
                    ax.annotate(
                        fr'${{\rm log}}M_{{\star}} = {logM}$',
                        xy=(0.0, 0.0),
                        xycoords='axes fraction',
                        xytext=(8, 40),
                        textcoords='offset points',
                        ha='left',
                        va='bottom',
                        fontsize=16,
                        color=overlay_color,
                    )
                    ax.annotate(
                        fr'$z$ = {z_str}',
                        xy=(0.0, 0.0),
                        xycoords='axes fraction',
                        xytext=(8, 24),
                        textcoords='offset points',
                        ha='left',
                        va='bottom',
                        fontsize=16,
                        color=overlay_color,
                    )
                    ax.annotate(
                        f'ID:{self.subhaloID}',
                        xy=(0.0, 0.0),
                        xycoords='axes fraction',
                        xytext=(8, 8),
                        textcoords='offset points',
                        ha='left',
                        va='bottom',
                        fontsize=16,
                        color=overlay_color,
                    )
            cbar_ax = fig.add_subplot(gs[:, 4 * nband :])
            cbar_ax.set_facecolor('white')
            cbar = fig.colorbar(im_last, cax=cbar_ax)
            cbar_text_color = 'black'
            cbar.set_label(cbar_label, fontsize=12, color=cbar_text_color)
            cbar.ax.tick_params(labelsize=10, colors=cbar_text_color)
            cbar.ax.yaxis.label.set_color(cbar_text_color)
            plt.subplots_adjust(wspace=0, hspace=0)
            out_path = os.path.join(save_dir, f'all_bands_view_{view_i:02d}.png')
            fig.savefig(
                out_path,
                dpi=150,
                facecolor=fig.get_facecolor(),
                bbox_inches='tight',
                pad_inches=0.05,
            )
            plt.close(fig)

    def __plot_sed(self, sed: np.ndarray, waveRange: list, logscale: bool=True,
                   savedFilename: Union[str, None]=None):
        
        minWave = waveRange[0]
        maxWave = waveRange[1]
        
        # to avoid the border effect
        maxWave95 = maxWave * 0.95
        idx = np.where(sed[:, 0] * 10**4 < maxWave95)[0]
        
        plt.figure()
        plt.plot(sed[idx, 0] * 10**4, sed[idx, 1], label='Total')
        if logscale:
            plt.xscale('log')
        # plt.legend(frameon=False)
        plt.xlabel(r'Wavelength $[\AA]$')
        plt.ylabel(r'$F_{\nu}\ [Jy]$')
        plt.xlim(waveRange[0], waveRange[1]) # SED is redshifted.
        plt.title(f'SubhaloID: {self.subhaloID}', fontsize=16)
        plt.tight_layout()
        
        if savedFilename is not None:
            plt.savefig(savedFilename)
            plt.close()
        else:
            plt.show()
        
    def __show_images(self, imgFilenames: list, SEDFilenames: list, 
                      subhaloID: int, survey: str):
        
        for i, (imgname, sedname) in enumerate(zip(imgFilenames, SEDFilenames)):
            img = Image.open(imgname)
            sed = Image.open(sedname)

            w_img, h_img = img.size
            w_sed, h_sed = sed.size

            combined_width = w_img + w_sed
            combined_height = max(h_img, h_sed)

            combined_img = Image.new('RGB', (combined_width, combined_height))

            combined_img.paste(img, (0, 0))
            combined_img.paste(sed, (w_img, 0))

            fig, ax = plt.subplots()
            ax.axis('off')
            plt.imshow(combined_img)
            # plt.suptitle(f'SubhaloID: {subhaloID} ({survey} View {i:02d})', 
            #              fontsize=12, y=0.7)
            plt.tight_layout()
            plt.show()
            
    def __get_environment(self) -> str:
        try:
            from IPython import get_ipython
            ip = get_ipython()
            if ip is None:
                return "CommandLine"
            elif "IPKernelApp" in ip.config:
                return "Jupyter"
            else:
                return "IPythonTerminal"
        except ImportError:
            return "CommandLine"



    def runPostprocess(self, showImages: bool=False):
        
        """
        This method performs the post-processing of the mock observations for the specified subhalo.
        
        Parameters
        ----------
        showImages : bool, optional
            If True, displays the generated images, and can only work in Jupyter environment.
    
        """
        
        if self.config['postProcessing']:
            
            self.logger.info('Run Postprocessing')
        
            surveys = self.config['surveys']
            
            for survey in surveys:
                
                self.logger.info(f'Begin postprocessing for {survey}')
                
                if self.config['outputSEDOnly']:
                    sedFilenames = [self.dataCubeDir + f'/skirt_view_{i:02d}_sed.dat' 
                                    for i in range(self.properties['numViews'])]
                    self.__saveSEDs(sedFilenames, survey)
                    
                else:
                    
                    start_time = time.time()
                    
                    filters = self.properties[f'filters_{survey}']
                    saveBaseDir = f'mock_{survey}/Subhalo_{self.subhaloID}'
                    os.makedirs(saveBaseDir, exist_ok=True)

                    dataCubeFilenames = [os.path.join(self.dataCubeDir, f'skirt_view_{i:02d}_total.fits')
                                        for i in range(self.properties['numViews'])]
                    
                    throughputs = self.__get_throughputs(survey)
                    PSFs = None
                    if self.config[f'includePSF_{survey}']:
                        PSFs = self.__get_PSFs(survey)

                    bkgNoise = None
                    if self.config[f'includeBkg_{survey}']:
                        
                        bkgNoise = {}
                        
                        if self.config[f'noiseType_{survey}'] == 'instrument':
                        
                            skyBkg = self.config[f'skyBkg_{survey}']
                            darkCurrent = self.config[f'darkCurrent_{survey}']
                            readOut = self.config[f'readOut_{survey}']
                            
                            skyBkg = extend(skyBkg, len(filters))
                            darkCurrent = extend(darkCurrent, len(filters))
                            readOut = extend(readOut, len(filters))
                            
                            bkgNoise = {'noiseType': 'instrument',
                                        'skyBkg': skyBkg,
                                        'darkCurrent': darkCurrent,
                                        'readOut': readOut, 
                                        'exposureTime': self.properties[f'exposureTime_{survey}'],
                                        'numExposure': self.properties[f'numExposure_{survey}']}

                            
                        elif self.config[f'noiseType_{survey}'] == 'limitingMagnitude':
                            limitMag = self.config[f'limitMag_{survey}']
                            limitSNR = self.config[f'limitSNR_{survey}']
                            limitAperture = self.config[f'limitAperture_{survey}']
                            zeroPoint = self.config[f'zeroPoint_{survey}']
                            
                            limitMag = extend(limitMag, len(filters))
                            limitSNR = extend(limitSNR, len(filters))
                            limitAperture = extend(limitAperture, len(filters))
                            zeroPoint = extend(zeroPoint, len(filters))
                            
                            bkgNoise = {'noiseType': 'limitingMagnitude',
                                        'limitMag': limitMag,
                                        'limitSNR': limitSNR,
                                        'limitAperture': limitAperture,
                                        'zeroPoint': zeroPoint, 
                                        'exposureTime': self.properties[f'exposureTime_{survey}'],
                                        'numExposure': self.properties[f'numExposure_{survey}']}                   

                    images = []
                    for i in range(self.properties['numViews']):
                        self.logger.info(f'Generating images for view {i}.')
                        images.append(self.__bandpass_images(dataCubeFilenames[i], survey,
                                                        throughputs, PSFs, bkgNoise))

                    self.logger.info('Saving bandpass images.')
                    self.__saveBandpassImages(images, survey)
                    self.__save_all_bands_montage(images, survey, filters)


                    sedFilenames = [self.dataCubeDir + f'/skirt_view_{i:02d}_sed.dat' 
                                    for i in range(self.properties['numViews'])]
                    self.logger.info('Saving SEDs')
                    self.__saveSEDs(sedFilenames, survey)

                    if self.config[f'imgDisplay_{survey}']:
                        enlarge_ratio = self.__survey_enlarge_ratio(survey)

                        if self.config[f'RGBImg_{survey}']:
                            RGBFilters = self.config[f'RGBFilters_{survey}']
                            RGBidx = [filters.index(RGBFilters[i]) for i in range(3)]
                            # res = self.properties[f'resolution_{survey}'][RGBidx[0]]
                            pixelscale = self.properties[f'angleRes_{survey}'][RGBidx[0]].value # pixel scale in arcsec
                            resolution = self.properties[f'resolution_{survey}'][RGBidx[0]].value # resolution in pc

                            for i in range(self.properties['numViews']):
                                
                                RGBImg = convert_to_rgb(images[i], RGBidx)

                                # from PIL import Image
                                # Image.fromarray(RGBImg).save(f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_view_{i:02d}_RGB.png')

                                savedFilename = f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_view_{i:02d}.png'
                                self.__plot_image(RGBImg, pixelscale, resolution, savedFilename=savedFilename,
                                                  enlarge_ratio=enlarge_ratio)
                        
                        else:
                            displayFilter = self.config[f'displayFilter_{survey}']
                            filteridx = filters.index(displayFilter)
                            pixelscale = self.properties[f'angleRes_{survey}'][filteridx].value # pixel scale in arcsec
                            resolution = self.properties[f'resolution_{survey}'][filteridx].value # resolution in pc
                        
                            for i in range(self.properties['numViews']):
                                img = images[i][filteridx]
                                savedFilename = f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_view_{i:02d}.png'
                                self.__plot_image(img, pixelscale, resolution, savedFilename=savedFilename,
                                                  enlarge_ratio=enlarge_ratio)

                    if self.config['displaySED']:
                        
                        for i in range(self.properties['numViews']):
                            
                            sed = np.loadtxt(sedFilenames[i])
                            savedFilename = f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_SED_view_{i:02d}.png'

                            logscale = self.config['displaySEDxlogscale']
                            minWave = self.properties['minWavelength'].to(u.angstrom).value
                            maxWave = self.properties['maxWavelength'].to(u.angstrom).value
                            # redshift = 0 if self.properties['viewRedshift'] == None else self.properties['viewRedshift']
                            # minWave = minWave * (1 + redshift)
                            # maxWave = maxWave * (1 + redshift)
                            
                            self.__plot_sed(sed,  waveRange=[minWave, maxWave], logscale=logscale,
                                        savedFilename=savedFilename)
                    
                    env = self.__get_environment()

                    if showImages:
                    
                        if env == 'CommandLine' or env == 'IPythonTerminal':
                            self.logger.info('Running in command line or IPython terminal, showImages can only be False.')
                            showImages = False
                        elif env == 'Jupyter':
                            showImages = showImages
                    
                    if showImages:
                        
                        imgFilenames = [f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_view_{i:02d}.png' 
                                        for i in range(self.properties['numViews'])]
                        sedFilenames = [f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_SED_view_{i:02d}.png' 
                                        for i in range(self.properties['numViews'])]
                        self.__show_images(imgFilenames, sedFilenames, self.subhaloID, survey)
                    
                    self.logger.info(f'Finish postprocessing for {survey}')
                    end_time = time.time()
                    duration = end_time - start_time
                    self.logger.info(f'Time taken to postprocess for {survey}: {duration:.2f} seconds')
                    