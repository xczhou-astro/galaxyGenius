import numpy as np
from astropy.io import fits
from astropy.convolution import convolve_fft
import astropy.units as u
import astropy.constants as const
from PIL import Image
from scipy.interpolate import interp1d
import os
from skimage.transform import rescale
from scipy.integrate import trapezoid
from .utils import *
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib_scalebar.dimension import _Dimension
from photutils.segmentation import make_2dgaussian_kernel
from typing import Union
from types import NoneType
import json

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
        self.properties = self.__get_properties()
        self.config = self.__read_configs()
        self.dataDir = self.config['dataDir']
        
    def __read_configs(self):
        
        with open(self.dataCubeDir + '/config.json', 'r') as file:
            config = json.load(file)
            
        return config
        
    def __get_properties(self) -> dict:
        path = os.path.join(self.dataCubeDir, 'properties.json')
        with open(path, 'r') as file:
            properties = json.load(file)
        
        return properties
            
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
            wavelength_scale = get_wavelength_scale(throughput_file)
            throughput = np.loadtxt(throughput_file)
            throughput[:, 0] = throughput[:, 0] * wavelength_scale
            throughputs.append(throughput)

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
            fwhm_in_arcsec = extend(fwhm_in_arcsec, numfilters)
            
            pixelScales = self.properties[f'angleRes_{survey}']

            stds = [fwhm / ps for fwhm, ps in zip(fwhm_in_arcsec, pixelScales)]
            
            PSFs = []
            for std in stds:
                size = std * 11
                size = int(size)
                if size % 2 == 0: # make sure the size is odd
                    size += 1
                kernel = make_2dgaussian_kernel(std, size)
                PSFs.append(kernel)
            
        return PSFs

    def __calculate_bandpass(self, img: np.ndarray, tran: np.ndarray, 
                              wave: np.ndarray, factor: float) -> np.ndarray:
        bandpass_img = factor * trapezoid(img * tran.reshape(-1, 1, 1) \
                                          * wave.reshape(-1, 1, 1), wave, axis=0)
        return bandpass_img
    
    def __bandpass_images(self, dataCube: str, survey: str, throughputs: list,
                          PSFs: Union[list, NoneType]=None, bkgNoise: Union[dict, NoneType]=None) -> list:
        
        file = fits.open(dataCube)
        wave_sed = file[1].data['grid_points'] * 10**4
        numExp = self.properties[f'numExposure_{survey}']
        exposureTime = self.properties[f'exposureTime_{survey}']
        areaMirror = np.pi * (self.properties[f'aperture_{survey}'] / 2)**2

        ps_in_sr = self.properties[f'ps_in_sr_{survey}']
        interpolate_ratios = self.properties[f'ratios_{survey}']
        pixelScales = self.properties[f'angleRes_{survey}']

        numfils = len(throughputs)

        image_arrs = []
        trans = []
        waves = []
        factors = []
        conversion_to_Jy = []
        gains = []
        pivots = []
        for i, thr in enumerate(throughputs):
            wave_min = np.min(thr[:, 0])
            wave_max = np.max(thr[:, 0])
            interp = interp1d(thr[:, 0], thr[:, 1])
            idx = np.where((wave_sed > wave_min) & (wave_sed < wave_max))[0]
            wave_in = wave_sed[idx]# * u.angstrom
            trans_in = interp(wave_in)
            image_arr = file[0].data[idx] / wave_in.reshape(-1, 1, 1)**2 # in flambda/u.sr
            converter = (u.MJy/u.sr) * const.c / u.angstrom**2 \
                        * numExp[i] * (exposureTime[i] * u.s) * (areaMirror * u.m**2) / (const.c * const.h) * u.angstrom**2
            converter = converter.to(u.sr ** -1)
            converter = converter.value
            gain = const.h / (areaMirror * u.m**2 * trapezoid(trans_in / wave_in, wave_in))
            gains.append(gain)
            image_arrs.append(image_arr)
            trans.append(trans_in)
            waves.append(wave_in)
            factors.append(converter)
            
            Jy_converter = trapezoid(trans_in * wave_in, wave_in) * u.angstrom**2
            numerator = trapezoid(trans_in, wave_in)
            denominator = trapezoid(trans_in * wave_in ** -2, wave_in)
            pivot = np.sqrt(numerator / denominator)
            pivots.append(pivot)
            Jy_converter = Jy_converter / (const.h * const.c) * areaMirror * u.m**2\
                            * numExp[i] * (exposureTime[i] * u.s)
            Jy_converter = (pivot * u.angstrom)**2 / const.c / Jy_converter
            Jy_converter = Jy_converter.to(u.Jy).value
            conversion_to_Jy.append(Jy_converter)

        bandpass_images = []
        for img, tran, wave, factor in zip(image_arrs, trans, waves, factors):
            bandpass_images.append(self.__calculate_bandpass(img, tran, wave, factor))

        
        # bandpass_images = Parallel(n_jobs=numfils, prefer='threads')(delayed(self.__calculate_bandpass)(img, tran, wave, factor)
        #                                            for img, tran, wave, factor in zip(image_arrs, trans, waves, factors))
        
        resized_imgs = []
        for i, ratio in enumerate(interpolate_ratios):
            resized_imgs.append(rescale(bandpass_images[i], ratio))

        bandpass_images = resized_imgs

        converted_imgs = []
        for i, ps_sr in enumerate(ps_in_sr):
            converted_imgs.append(bandpass_images[i] * ps_sr)
        
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
                
                skyBkg = np.array(bkgNoise['skyBkg'])
                darkCurrent = np.array(bkgNoise['darkCurrent'])
                readOut = np.array(bkgNoise['readOut'])
                exposureTime = np.array(bkgNoise['exposureTime'])
                numExposure = np.array(bkgNoise['numExposure'])

                images_with_bkg_in_electron = []
                for i, img in enumerate(bandpass_images):
                    mean = (skyBkg[i] + darkCurrent[i]) * exposureTime[i] * numExposure[i]
                    img = img + mean
                    img = np.random.poisson(img).astype(np.int32)
                    for _ in range(numExposure[i]):
                        readnoise = np.random.normal(loc=0, scale=readOut[i], size=img.shape)
                        readnoise = np.around(readnoise, 0).astype(np.int32)
                        img = img + readnoise
                    img = img.astype(np.float32) - mean.astype(np.float32)
                    images_with_bkg_in_electron.append(img)
                    
                if self.config['imageUnit'] == 'electron':
                    images_with_bkg = images_with_bkg_in_electron
                elif self.config['imageUnit'] == 'flux':
                    images_with_bkg = []
                    for i, img in enumerate(images_with_bkg_in_electron):
                        img = img * conversion_to_Jy[i]
                        images_with_bkg.append(img)
            
            elif bkgNoise['noiseType'] == 'limitingMagnitude':
                
                limitMag = np.array(bkgNoise['limitMag'])
                limitSNR = np.array(bkgNoise['limitSNR'])
                limitAperture = np.array(bkgNoise['limitAperture'])
                zeroPoint = np.array(bkgNoise['zeroPoint'])
                exposureTime = np.array(bkgNoise['exposureTime'])
                numExposure = np.array(bkgNoise['numExposure'])
                
                pixelScales = np.array(pixelScales)
                
                npix = (np.pi * (limitAperture / 2)**2) / (pixelScales**2)
                
                images_with_bkg_in_Jy = []
                for i, img in enumerate(bandpass_images):
                    img_in_Jy = img * conversion_to_Jy[i]
                    first_term = (gains[i] * img_in_Jy * u.Jy) / (numExposure[i] * exposureTime[i] * u.s)
                    # first_term = (gains[i] * img_in_Jy * u.Jy * pivots[i]**2 * u.angstrom**2) / (numExposure[i] * exposureTime[i] * u.s * const.c)
                    second_term = ((1 / limitSNR[i]) * 10**((zeroPoint[i] - limitMag[i]) / 2.5) * u.Jy)**2 * (1 / npix[i])
                    third_term = (gains[i]) / (numExposure[i] * exposureTime[i] * u.s * npix[i]) * 10**((zeroPoint[i] - limitMag[i]) / 2.5) * u.Jy
                    # third_term = (gains[i] * pivots[i]**2 * u.angstrom**2) / (numExposure[i] * exposureTime[i] * u.s * const.c * npix[i]) \
                    #     * 10**((zeroPoint[i] - limitMag[i]) / 2.5) * u.Jy
                    
                    total = (first_term + second_term - third_term).to(u.Jy**2).value
                    
                    noise = np.sqrt(total) # in Jy
                    noises = np.random.normal(loc=0, scale=noise, size=img.shape)
                    images_with_bkg_in_Jy.append(img + noises)
                    
                if self.config['imageUnit'] == 'electron':
                    images_with_bkg = []
                    for i, img in enumerate(images_with_bkg_in_Jy):
                        img = img / conversion_to_Jy[i]
                        images_with_bkg.append(img)
                        
                elif self.config['imageUnit'] == 'flux':
                    images_with_bkg = images_with_bkg_in_Jy
        
        return images_with_bkg

            
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
            header['MASS'] = (self.properties['stellarMass'], 'Subhalo stellar mass, in log10 scale (Msun)')
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
            header['EXPTIME'] = (self.properties[f'exposureTime_{survey}'][i], 'Exposure time, in s')
            header['EXPNUM'] = (self.properties[f'numExposure_{survey}'][i], 'Number of exposures')
            header['APERTURE'] = (self.properties[f'aperture_{survey}'], 'Aperture size, in meter')
            header['UNIT'] = (imageUnit, imageUnitComment)
            header['REDSHIFT'] = (self.properties['redshift'], 'Redshift')
            header['COSMO'] = (self.properties['cosmology'], 'Cosmology')
            header['BOXSIZE'] = (self.properties['boxLength'], 'Box size, in pc')
            header['FOV'] = (self.properties['fieldOfView'], 'Field of view, in arcsec')
            header['PHYFOV'] = (self.properties['fovSize'], 'Physical field of view, in pc')
            header['LUMIDIS'] = (self.properties['lumiDis'], 'Luminosity distance, in Mpc')
            header['RESOL'] = (self.properties[f'resolution_{survey}'][i], 'Pixel scale, in pc')
            header['PS'] = (self.properties[f'angleRes_{survey}'][i], 'Pixel scale, in arcsec')
            
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
            header['REDSHIFT'] = (self.properties['redshift'], 'Redshift')
            header['lumiDis'] = (self.properties['lumiDis'], 'luminosity distance, in Mpc')
            header['FOV'] = (self.properties['fieldOfView'], 'Field of view, in arcsec')
            header['PHYFOV'] = (self.properties['fovSize'], 'Physical field of view, in pc')
            header['VIEW'] = (i, 'View index')
            
            hdulist.append(tbhdu)
            
        savedSEDName = f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_SEDs.fits'
        hdulist.writeto(savedSEDName, overwrite=True)

    def __plot_image(self, image: np.ndarray, pixelscale: float,
                     savedFilename: Union[str, NoneType]=None):

        numPixels = np.array(image).shape[0]
        
        z = str(np.around(self.properties['redshift'], 2))
        logM = str(np.around(self.properties['stellarMass'], 1))
        
        fig, ax = plt.subplots()
        ax.axis('off')
        im = ax.imshow(image)
        
        scalebarSize = 0.25 * numPixels * pixelscale
        if scalebarSize > 60:
            scalebarUnit = '${^{\prime}}$'
            scalebarSize = np.around(scalebarSize / 60, 2)
            ps = pixelscale / 60
        else:
            scalebarUnit = '${^{\prime\prime}}$'
            scalebarSize = np.around(scalebarSize, 2)
            ps = pixelscale
            
        angle_dim = AngleDimension()
        
        scalebar = ScaleBar(ps, scalebarUnit, dimension=angle_dim, 
                            fixed_value=scalebarSize, fixed_units=scalebarUnit, frameon=False,
                            location='lower right', scale_loc='top',
                            color='white', font_properties={'size': 12})
        ax.add_artist(scalebar)
        ax.text(x=0.05, y=0.15, s=fr'${{\rm log}}M_{{\star}} = {logM}$',
                         fontsize=12, transform=ax.transAxes, color='white')
        ax.text(x=0.05, y=0.1, s=fr'$z$ = {z}', fontsize=12,
                transform=ax.transAxes, color='white')
        ax.text(x=0.05, y=0.05, s=f'ID:{self.subhaloID}', fontsize=12,
                transform=ax.transAxes, color='white')
        # divider = make_axes_locatable(ax)
        # cax = divider.append_axes('right', size='5%', pad=0.05)
        # fig.colorbar(im, cax=cax, label=unit)
        if savedFilename is not None:
            plt.savefig(savedFilename)
            plt.close()
        else:
            plt.show()

    def __plot_sed(self, sed: np.ndarray, waveRange: list, logscale: bool=True,
                   savedFilename: Union[str, NoneType]=None):
        plt.figure()
        plt.plot(sed[:, 0] * 10**4, sed[:, 1], label='Total')
        if logscale:
            plt.xscale('log')
        # plt.legend(frameon=False)
        plt.xlabel(r'Wavelength $[\AA]$')
        plt.ylabel(r'$F_{\nu}\ [Jy]$')
        plt.xlim(waveRange[0], waveRange[1]) # SED is redshifted.
        
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
            plt.suptitle(f'SubhaloID: {subhaloID} ({survey} View {i:02d})', 
                         fontsize=12, y=0.7)
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
            
            print('Run Postprocessing')
        
            surveys = self.config['surveys']
            
            for survey in surveys:
                
                print(f'Begin postprocessing for {survey}')
                
                if self.config['outputSEDOnly']:
                    sedFilenames = [self.dataCubeDir + f'/skirt_view_{i:02d}_sed.dat' 
                                    for i in range(self.properties['numViews'])]
                    self.__saveSEDs(sedFilenames, survey)
                    
                else:
        
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
                        images.append(self.__bandpass_images(dataCubeFilenames[i], survey,
                                                        throughputs, PSFs, bkgNoise))

                    self.__saveBandpassImages(images, survey)

                    sedFilenames = [self.dataCubeDir + f'/skirt_view_{i:02d}_sed.dat' 
                                    for i in range(self.properties['numViews'])]
                    self.__saveSEDs(sedFilenames, survey)

                    if self.config[f'imgDisplay_{survey}']:

                        if self.config[f'RGBImg_{survey}']:
                            RGBFilters = self.config[f'RGBFilters_{survey}']
                            RGBidx = [filters.index(RGBFilters[i]) for i in range(3)]
                            # res = self.properties[f'resolution_{survey}'][RGBidx[0]]
                            pixelscale = self.properties[f'angleRes_{survey}'][RGBidx[0]] # pixel scale in arcsec

                            for i in range(self.properties['numViews']):
                                
                                # print(images[i][RGBidx[0]].shape)
                                # print(images[i][RGBidx[1]].shape)
                                # print(images[i][RGBidx[2]].shape)
                                
                                RGBImg = convert_to_rgb(images[i], RGBidx)
                                savedFilename = f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_view_{i:02d}.png'
                                self.__plot_image(RGBImg, pixelscale, savedFilename=savedFilename)
                        
                        else:
                            displayFilter = self.config[f'displayFilter_{survey}']
                            filteridx = filters.index(displayFilter)
                            # res = self.properties[f'resolution_{survey}'][filteridx]
                            pixelscale = self.properties[f'angleRes_{survey}'][filteridx] # pixel scale in arcsec
                        
                            for i in range(self.properties['numViews']):
                                img = images[i][filteridx]
                                savedFilename = f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_view_{i:02d}.png'
                                self.__plot_image(img, pixelscale, savedFilename=savedFilename)

                    if self.config['displaySED']:
                        
                        for i in range(self.properties['numViews']):
                            
                            sed = np.loadtxt(sedFilenames[i])
                            savedFilename = f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_SED_view_{i:02d}.png'

                            logscale = self.config['displaySEDxlogscale']
                            minWave = self.properties['minWavelength']
                            maxWave = self.properties['maxWavelength']
                            redshift = self.properties['redshift']
                            minWave = minWave * 10**4 # angstrom
                            maxWave = maxWave * 10**4 # angstrom
                            # print(logscale)
                            self.__plot_sed(sed,  waveRange=[minWave, maxWave], logscale=logscale,
                                        savedFilename=savedFilename)
                    
                    env = self.__get_environment()

                    if showImages:
                    
                        if env == 'CommandLine' or env == 'IPythonTerminal':
                            print('Running in command line or IPython terminal, showImages can only be False.')
                            showImages = False
                        elif env == 'Jupyter':
                            showImages = showImages
                    
                    if showImages:
                        
                        imgFilenames = [f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_view_{i:02d}.png' 
                                        for i in range(self.properties['numViews'])]
                        sedFilenames = [f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_SED_view_{i:02d}.png' 
                                        for i in range(self.properties['numViews'])]
                        self.__show_images(imgFilenames, sedFilenames, self.subhaloID, survey)
                    
                    print(f'Finish postprocessing for {survey}')
                    