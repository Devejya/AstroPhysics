import numpy as np
import matplotlib.pyplot as plt
import datetime
from astropy.io import ascii,fits
from scipy.stats import linregress
from astropy import constants


class A2:
    def __init__(self, spectrum_file: str, bands: np.array, vega_file: str, 
                temp_range: np.array, d_lambda: np.array):
        self.spectrum_file = spectrum_file
        self.bands = bands
        self.vega_file = vega_file
        self.temp_range = temp_range
        self.d_lambda = d_lambda
        self.ll_range = np.arange(1000,20000,d_lambda)
        self.h = constants.h.value
        self.c = constants.c.value
        self.k = constants.k_B.value
        self.spectrum_data = ascii.read(spectrum_file)
        self.wavelengths = self.spectrum_data['Wavelength']
        self.flux = self.spectrum_data['Flux']
        self.vega_data = ascii.read(vega_file)
        self.vega_wavelength = self.vega_data['Wavelength']
        self.vega_flux = self.vega_data['Flux']


    def find_flux(self, wavelengths: np.array, flux: np.array, bands: np.array)-> dict:
        """
        returns the arrray of photon flux given a list wavelengths, fluxes and bands
        input args:
            wavelengths: array of wavelengths (Angstrom)
            flux: array of fluxes (W/m^2/A)
            bands: lower, upper wavelength limits for N filters
        """
        output = np.array([])
        n_filter = np.shape(bands)[0]
        for i in np.arange(0, n_filter):
            wavelength_range = np.where((wavelengths > bands[i][0]) & (wavelengths < bands[i][1]))
            ll = wavelengths[wavelength_range]
            dlambda = ll[1::]-ll[0:-1]
            dlambda = np.append(dlambda,dlambda[-1])
            photon_flux = 1./(h*c) * np.sum(flux[wavelength_range] * flux[wavelength_range] 
                        * dlambda)*1.e-10
            output = np.append(output,photon_flux)
        return dict(['g','r','i'],output)

    
    def find_magnitudes(self):
        "return the g,r,i Vega magnitudes"
        self.photon_flux = self.find_flux(self.wavelengths, self.flux, self.bands)
        self.vega_flux = self.find_flux(self.vega_wavelength, self.vega_flux, self.bands)
        filters = ['g', 'r', 'i']
        magnitudes = []
        for i in range(len(photon_flux)):
            magnitudes.append(-2.5*np.log10(self.photon_flux[i], self.vega_flux[i]))
            
        return dict(zip(filters, magnitudes))

    
    def blackbody_flux_mag(self) -> tuple:
        "return the blackbody flux and magnitudes"
        ll_range_m = self.ll_range * 1.e-10
        bb_gmag = np.arrau([])
        bb_rmag = np.array([])
        bb_imag = np.array([])

        for T in self.temp_range():
            B_flux = 1.e-10*2*h*c*c/(ll_range_m**5)/(np.exp((h*c)/(ll_range_m*k*T))-1)
            bb_fluxs = self.find_flux(self.ll_range, B_flux, self.bands)
            bb_gmag = np.append(bb_gmag, -2.5*np.log10(bb_fluxs[0]/self.vega_flux[0]))
            bb_rmag = np.append(bb_rmag, -2.5*np.log10(bb_fluxs[1]/self.vega_flux[1]))
            bb_imag = np.append(bb_imag, -2.5*np.log10(bb_fluxs[2]/self.vega_flux[2]))

        self.bb_flux = bb_fluxs
        self.bb_mags = dict(zip(['g', 'r', 'i'], [bb_gmag, bb_rmag, bb_imag]))
        return self.bb_flux, self.bb_mags

    
    def plot_data(self):
        ''
        fig, ax = plt.subplots()
        plt.plot(self.temp_range, self.bb_mags['g']-self.bb_mags['r'])
        plt.ylabel('(g-r)')
        plt.xlabel('T (K)')
        plt.show()

        fig,ax=plt.subplots()
        plt.plot(self.temp_range, self.bb_mags['r']-self.bb_mags['i'])
        plt.ylabel('(r-i)')
        plt.xlabel('T (K)')
        plt.show()

if __name__ == "__main__":

    Ass2 = A2("phy270_ass2_spectrum.txt", np.array([[4000,5400],[5500,6850],[6700,8250]]), 
                "vega.txt", np.arange(5000,20000,500), 5)

    photon_flux = Ass2.photon_flux
    print("Photon Flux: ", photon_flux)

    vega_flux = Ass2.vega_flux
    print("Vega Flux: ", vega_flux)

    vega_mags = Ass2.find_magnitudes()
    print("Vega magnitudes: ", vega_mags)

    bb_flux_mags = Ass2.blackbody_flux_mag()
    print("BB Flux: ", bb_flux_mags[0], "\n", "BB Magnitudes: ", bb_flux_mags[1])

    Ass2.plot_data()










    

def q1():
	h=constants.h.value
	k=constants.k_B.value
	c=constants.c.value
# 2D array that has (lambda_min,lambda_max) for each filter under consideration
	band_definitions=np.array([[4000,5400],[5500,6850],[6700,8250]])
	data=ascii.read('phy270_ass2_spectrum.txt')
	lam=data['Wavelength'] 
	flux=data['Flux']
	fluxes=calc_photon_flux(lam,flux,band_definitions)
	print ('Photon flux: g=',fluxes[0],'r=',fluxes[1],'i=',fluxes[2])

	vega=ascii.read('vega.txt')
	Vega_lambda=vega['Wavelength']
	Vega_flambda=vega['Flux'] # convert to W/m^2/s
	Vega_fluxes=calc_photon_flux(Vega_lambda,Vega_flambda,band_definitions)
	print (Vega_fluxes)
	gmag=-2.5*np.log10(fluxes[0]/Vega_fluxes[0])
	rmag=-2.5*np.log10(fluxes[1]/Vega_fluxes[1])
	imag=-2.5*np.log10(fluxes[2]/Vega_fluxes[2])
	print (gmag,rmag,imag)
	Tarr=np.arange(5000,20000,500)
	dl=5
	ll=np.arange(1000,20000,dl)
	llm=ll*1.e-10
	BBgmag=np.array([])
	BBrmag=np.array([])
	BBimag=np.array([])
	for T in Tarr:
		B=1.e-10*2*h*c*c/(llm**5)/(np.exp((h*c)/(llm*k*T))-1) # in W/m^2/A
		BB_fluxes=calc_photon_flux(ll,B,band_definitions)
		BBgmag=np.append(BBgmag,-2.5*np.log10(BB_fluxes[0]/Vega_fluxes[0]))
		BBrmag=np.append(BBrmag,-2.5*np.log10(BB_fluxes[1]/Vega_fluxes[1]))
		BBimag=np.append(BBimag,-2.5*np.log10(BB_fluxes[2]/Vega_fluxes[2]))
	fig,ax=plt.subplots()
	print (np.size(Tarr),Tarr)
	print (np.size(BBgmag-BBrmag),BBgmag-BBrmag)
	plt.plot(Tarr,BBgmag-BBrmag)
	plt.ylabel('(g-r)')
	plt.xlabel('T (K)')
	plt.show()
	fig,ax=plt.subplots()
	plt.plot(Tarr,BBrmag-BBimag)
	plt.ylabel('(r-i)')
	plt.xlabel('T (K)')
	plt.show()

def photon_flux(wavelengths: list, fluxs: np.array, bands: np.array):
    '''
    returns the arrray of photon flux given a list wavelengths, fluxes and bands
    input args:
        wavelengths: array of wavelengths (Angstrom)
        fluxs: array of fluxes (W/m^2/A)
        bands: lower, upper wavelength limits for N filters
    '''
    h = constants.h.value
    c=constants.c.value
    output=np.array([])
# Determine number of fluxes to compute based on shape of band_definitions.	
    nfilter=np.shape(band_definitions)[0]
# Calculate wavelength interval.  
    for i in np.arange(0,nfilter):
        lamrange=np.where((l>band_definitions[i][0]) & (l<band_definitions[i][1]))
        ll=l[lamrange]
        dlambda=ll[1::]-ll[0:-1]
        dlambda=np.append(dlambda,dlambda[-1])
        # Calcualte photon_flux in photons/m2/s
        photon_flux=1./(h*c)*np.sum(l[lamrange]*f[lamrange]*dlambda)*1.e-10
        output=np.append(output,photon_flux)
    return output

if __name__ == "__main__":

    #define constants
