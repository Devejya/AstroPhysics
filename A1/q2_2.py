import datetime
import numpy as np
from astropy.io import ascii
from scipy.stats import linregress


def proper_motion(filename: str) -> float:
    data = ascii.read(filename, names=['Days', 'RA', 'Dec'])

    days_sample = np.arange(4)*365.25

    sampleRA = np.interp(days_sample,data['Days'],data['RA'])
    sampleDec = np.interp(days_sample,data['Days'],data['Dec'])

    ra = linregress(days_sample,sampleRA)
    dec = linregress(days_sample,sampleDec)

    slope_ra = ra[0]
    slope_dec = dec[0]

    proper_motion_value = np.sqrt(slope_ra**2+slope_dec**2)*365.25
    
    return proper_motion_value

def parallax(filename: str) -> float:

    data = ascii.read(filename, names=['Days', 'RA', 'Dec'])

    days_sample = np.arange(4)*365.25

    sampleRA = np.interp(days_sample,data['Days'],data['RA'])
    sampleDec = np.interp(days_sample,data['Days'],data['Dec'])

    ra = linregress(days_sample,sampleRA)
    dec = linregress(days_sample,sampleDec)

    # Fit a line to RA as a function of Dec for the annually-separated points.  Then calculaate 
    # the maximum perpendicular distance from this line to get the parallax	
    ra_v_dec =linregress(ra,dec)
    ra_v_dec_slope = ra_v_dec[0]
    ra_v_dec_intercept = ra_v_dec[1]

    x_min = (ra_v_dec_slope*(data['Dec']-ra_v_dec_intercept)+data['RA'])/(1+ra_v_dec_slope**2)
    y_min = x_min*ra_v_dec_slope + ra_v_dec_intercept

    return np.max(np.sqrt((x_min-data['RA'])*(x_min-data['RA']) + (y_min-data['Dec'])*(y_min-data['Dec'])))


if __name__ == "__main__":
    filename 
