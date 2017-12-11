import os
import glob
import logging
import operator as op

from mpi4py import MPI

import numpy as np

import plio.utils
from plio.utils import log
from plio.utils.utils import check_file_exists
from plio.io import io_gdal, io_hdf, io_json
from plio.date import astrodate, julian2ls, julian2season

import pvl

logger = logging.getLogger(__name__)


def enum(*sequential, **named):
    """Handy way to fake an enumerated type in Python
    http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)


def getnearest(iterable, value):
    """
    Given an iterable, get the index nearest to the input value

    Parameters
    ----------
    iterable : iterable
               An iterable to search

    value : int, float
            The value to search for

    Returns
    -------
        : int
          The index into the list
    """
    return min(enumerate(iterable), key=lambda i: abs(i[1] - value))


def checkbandnumbers(bands, checkbands):
    """
    Given a list of input bands, check that the passed
    tuple contains those bands.

    In case of THEMIS, we check for band 9 as band 9 is the temperature
    band required to derive thermal temperature.  We also check for band 10
    which is required for TES atmosphere calculations.

    Parameters
    -----------
    bands       tuple of bands in the input image
    checkbands  list of bands to check against

    Returns
    --------
    Boolean     True if the bands are present, else False
    """
    for c in checkbands:
        if c not in bands:
            return False
    return True


def checkdeplaid(incidence):
    """
    Given an incidence angle, select the appropriate deplaid method.

    Parameters
    -----------
    incidence       float incidence angle extracted from the campt results.

    """
    if incidence >= 95 and incidence <= 180:
        return 'night'
    elif incidence >=90 and incidence < 95:
        logger.error("Incidence angle is {}.  This is a twilight image, using night time deplaid.".format(incidence))
        return 'night'
    elif incidence >= 85 and incidence < 90:
        logger.error("Incidence angle is {}.  This is a twilight image, using daytime deplaid".format(incidence))
        return 'day'
    elif incidence >= 0 and incidence < 85:
        logger.error("Incidence angle is {}.  This is a daytime image, you may not want to use this product.".format(incidence))
        return 'day'
    else:
        logger.error("Incidence does not fall between 0 and 180.")
        return False


def check_change_by(iterable, by=1, piecewise=False):
    """
    Check that a given iterable increases by one with each index

    Parameters
    ----------
    iterable : iterable
               Any Python iterable object

    by : int
         The value by which the iterables should be increasing

    piecewise : boolean
                If false, return a boolean for the entire iterable,
                else return a list with elementwise monotinicy checks

    Returns
    -------
    monotonic : bool/list
                A boolean list of all True if monotonic, or including
                an inflection point
    """
    increasing = [True]  + [(by == i[1] - i[0]) for i in zip(iterable,iterable[1:])]
    if piecewise:
        return increasing
    else:
        return all(increasing)


def checkmonotonic(iterable, op=op.gt, piecewise=False):
    """
    Check if a given iterable is monotonically increasing.

    Parameters
    ----------
    iterable : iterable
                Any Python iterable object

    op : object
         An operator.operator object, e.g. op.gt (>) or op.geq(>=)

    piecewise : boolean
                If false, return a boolean for the entire iterable,
                else return a list with elementwise monotinicy checks

    Returns
    -------
    monotonic : bool/list
                A boolean list of all True if monotonic, or including
                an inflection point
    """
    monotonic =  [True] + [x < y for x, y in zip(iterable, iterable[1:])]
    if piecewise == True:
        return monotonic
    else:
        return all(monotonic)

def convert_mean_pressure(elevation):
    """
    Convert from raw elevation, in km, to pressure in Pascals using
    Hugh Kieffer's algorithm.

    689.7 is the constant pressure at sea level

    Parameters
    -----------
    elevation : float or ndarray
                elevation in km

    Returns
    --------
      : float
        Pressure in Pascals
    """
    return 689.7 * np.exp(-elevation / 10.8)


def find_in_dict(obj, key):
    """
    Recursively find an entry in a dictionary

    Parameters
    ----------
    obj : dict
          The dictionary to search
    key : str
          The key to find in the dictionary

    Returns
    -------
    item : obj
           The value from the dictionary
    """
    if key in obj:
        return obj[key]
    for k, v in obj.items():
        if isinstance(v,dict):
            item = find_in_dict(v, key)
            if item is not None:
                return item


# note that this decorator ignores **kwargs
def memoize(obj):
    cache = obj.cache = {}
    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        if args not in cache:
            cache[args] = obj(*args, **kwargs)
        return cache[args]
    return memoizer


def extract_ancillary_data(job, temperature,parameters, workingpath, shape, reference_dataset):
    """
    For all ancillary data sets, extract the requested spatial extent

    Parameters
    ----------
    job : dict
        Job specification dictionary

    temperature : object
                Plio GeoDataset object

    parameters : dict
                of extent and time parameters
    """
    ancillarydata = job['ancillarydata']
    for k, v in ancillarydata.items():
        ancillarydata[k] = v

    # Flag to punt on the image if some ancillary data load fails.
    ancillary_failure = False

    #Iterate through the ancillary data.  Clip and resample to the input image
    for k, v in ancillarydata.items():
        if isinstance(v, int) or isinstance(v, float):
            #The user wants a constant value to be used
            arr = np.empty(shape, dtype=np.float32)
            logger.debug('{} set to a constant value, {}'.format(k, v))
            del arr
        else:
            basename = os.path.basename(v)
            root, extension = os.path.splitext(basename)
            tif = os.path.join(workingpath, root + '.tif')

            if v == 'montone':  # Custom dust opacity
                startls = parameters['startlsubs'][0]
                startmartianyear = int(parameters['startmartianyear'][0])
                files = glob.glob('/scratch/jlaura/KRC/basemaps/tau_geotiff/MY{}*'.format(startmartianyear))

                if len(files) == 0:
                    logger.error('Requested an image with Mars Year {}.  No Montabone data available for that Mars Year'.format(startmartianyear))
                    ancillary_failure = True
                    continue
                ls = {}
                for t, f in enumerate(files):
                    base, ext = os.path.splitext(f)
                    ls[float(base.split('_')[2])] = t
                keys = []
                for key in ls.keys():
                    try:
                        keys.append(float(key))
                    except: pass
                key = min(keys, key=lambda x: abs(x-startls))
                v = files[ls[key]]

            if v == 'tes':
                startls = startlsubs[0]
                files = glob.glob('/scratch/jlaura/KRC/basemaps/tes_opac/*')
                ls = {}
                for t, f in enumerate(files):
                    base, ext = os.path.splitext(f)
                    ls[float(base.split('_')[-1][2:])] = t
                keys = []
                for key in ls.keys():
                    try:
                        keys.append(float(key))
                    except: pass

                key = min(keys, key=lambda x: abs(x-startls))
                v = files[ls[key]]

            #Clip and resample the image to the correct resolution
            v = io_gdal.GeoDataset(v)
            io_gdal.match_rasters(reference_dataset, v, tif)

            #Read the resampled tif and extract the array
            ancillarydata[k] = io_gdal.GeoDataset(tif)
            logger.debug('Dataset {} extract.'.format(v))

    if ancillary_failure:
        print('FAILED TO EXTRACT ANCILLARY DATA')
        MPI.COMM_WORLD.Abort(1)
        sys.exit()

    return ancillarydata


def extract_temperature(isiscube, reference_dataset=None):
    """
    Extract the temperature data from the processed ISIS cube.

    Parameters
    ----------
    isiscube : str
               PATH to an ISIS cube to extract
    """
    temperature = io_gdal.GeoDataset(isiscube)
    processing_resolution = temperature.pixel_width
    tempshape = list(temperature.raster_size)[::-1]
    logger.info('Themis temperature data has {} lines and {} samples'.format(tempshape[0], tempshape[1]))
    srs = temperature.spatial_reference.ExportToWkt()
    logger.info('The input temperature image projection is: {}'.format(srs))
    return temperature

def extract_metadata(isiscube, parameters):
    """
    Given a Davinci processed, level 2 cube, extract the necessary
    metadata from the header to support clipping supplemental data sets.
    """
    header = pvl.load(isiscube)

    ulx = maxlat = find_in_dict(header, 'MaximumLatitude')
    uly = find_in_dict(header, 'MinimumLongitude')
    lrx = minlat = find_in_dict(header, 'MinimumLatitude')
    lry = find_in_dict(header, 'MaximumLongitude')
    logger.debug('Input TI image LAT range is {} to {}'.format(minlat, maxlat))

    starttime = find_in_dict(header, 'StartTime')
    stoptime = find_in_dict(header, 'StopTime')

    #Convert UTC to Julian
    starttime_julian = astrodate.utc2jd(starttime)
    stoptime_julian = astrodate.utc2jd(stoptime)
    logger.debug('Input TI image time range is {} to {} (Julian)'.format(starttime_julian, stoptime_julian))

    #LsubS
    startlsubs, startmartianyear = julian2ls.julian2ls(starttime_julian)
    stoplsubs, stopmartianyear = julian2ls.julian2ls(stoptime_julian)
    season, startseason, stopseason = julian2season.j2season(startlsubs)

    logger.debug('Input TI image time range is {} / {} to {} / {} (LsubS)'.format(startlsubs[0],startmartianyear[0],
                                            stoplsubs[0], stopmartianyear[0]))
    season = season[0]
    logger.debug('Season: {}, Start Season: {}, Stop Season {}'.format(season, startseason, stopseason))

    parameters['startseason'] = startseason
    parameters['stopseason'] = stopseason
    parameters['season'] = season
    parameters['startlsubs'] = startlsubs
    parameters['starttime'] = starttime
    parameters['stoptime'] = stoptime
    parameters['startlatitude'] = lrx
    parameters['stoplatitude'] = ulx
    parameters['startmartianyear'] = startmartianyear
    parameters['stopmartianyear'] = stopmartianyear

    return parameters


def cost_func(variables, obs3, obs9, rock3, rock9):
    """ Basic cost function used in differential evolution technique to find the
        correct balance between fine component and rock.
    Parameters
    ----------
    variables: list of tuples
               A list of values that indicate bounds for alpha, fine3, and fine9

    obs3:      float
               The floating point value representing the observed value of a
                 pixel from the third band of an input image.

    obs9:      float
               The floating point value representing the observed value of a
                 pixel from the ninth band of an input image.

    rock3:     float
               The expected value of rock in band 3.

    rock9:     float
               The expected value of rock in band 9.

    Returns
    -------
    z:         float
               The cost of the current configuration as described by the
                 mathematical cost function.
    """

    alpha = variables[0]
    fine3 = variables[1]
    fine9 = variables[2]

    if fine9 < fine3:
        return 1e10

    delta_obs = abs(obs3-obs9)
    delta_rock = abs(rock3 - rock9)
    delta_fine = abs(fine3 - fine9)

    z = abs((delta_rock * alpha + delta_fine * (1-alpha)) - delta_obs)

    return z

def generate_rad_image(temp_image, band_num):
    """
    Given a temperature image, generate a new image that is filled with the radiance value for
    rocks from a given band number. The new radiance image is the same size as the
    temperature image.

    Parameters
    ----------
    temp_image : ndarray
               Array representation of the image

    band_num : int
               Band number that needs to be extracted

    Returns
    ----------
    rad_image : ndarray
               Array representation of the new rad image

    """
    band_values = {
    "1": 0.000173866,
    "2": 0.000173866,
    "3": 0.000266925,
    "4": 0.000310375,
    "5": 0.000352625,
    "6": 0.000382479,
    "7": 0.000396245,
    "8": 0.000398744,
    "9": 0.000393225,
    "10": 0.000348427}

    if int is type(band_num):
        band_num = str(band_num)

    radiance_value = band_values[band_num]
    rad_image = np.full(temp_image.shape, radiance_value)
    return rad_image

def extract_band(job, image, band_num):
    """
    Extract the temperature data from the processed ISIS cube.

    Parameters
    ----------
    job : dict
               Job specification dictionary

    image : str
               PATH to an ISIS cube to extract bands from

    band_num : int
               Band number that needs to be extracted

    Returns
    ----------
    : ndarray
               Array representation of the extracted band
    """
    header = pvl.load(job['images'])
    bands = find_in_dict(header, 'BAND_BIN_BAND_NUMBER')

    for i, band in enumerate(bands):
        if band_num == band:
            geo_image = io_gdal.GeoDataset(image)
            return geo_image.read_array(band = i + 1)

def extract_latlon_transform(isiscube, job):
    """
    Given an ISIS cube, extract the upper left (x, y) coords and the height and width

    Parameters
    ----------
    isiscube : str
               PATH to an ISIS cube to use for latlon to pixel translation

    job : dict
               Job specification dictionary

    Returns
    ----------
    xoff : int
               x coordinate of the upper left pixel

    yoff : int
               y coordinate of the upper left pixel

    width : int
               Width of the new area

    height : int
               Height of the new area
    """
    isiscube_geodata = io_gdal.GeoDataset(isiscube)
    lry, uly = job["lat_extent"]
    ulx, lrx = job["lon_extent"]

    ul_coords = isiscube_geodata.latlon_to_pixel(uly, ulx)
    lr_coords = isiscube_geodata.latlon_to_pixel(lry, lrx)

    xoff = ul_coords[0]
    yoff = ul_coords[1]

    width = abs(lr_coords[0] - xoff)
    height = abs(lr_coords[1] - yoff)
    return xoff, yoff, width, height
