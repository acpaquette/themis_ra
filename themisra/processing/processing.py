import os
import logging
import operator as op

from mpi4py import MPI

import numpy as np

import plio.utils
from plio.utils import log
from plio.utils.utils import check_file_exists

import pvl

import themisra.utils.utils as util
from themisra.wrappers import pipelinewrapper, isiswrapper

def process_header(job):
    """
    Given the input image and job instructions, check that the necessary
    header information is present to process.

    Parameters
    ----------
    job : dict
          Containing the PATH to an images

    Returns
    -------
    job : dict
          With updated, image header specific values

    """

    header = pvl.load(job['images'])
    bands = util.find_in_dict(header, 'BAND_BIN_BAND_NUMBER')
    #bands = header['BAND_BIN_BAND_NUMBER']

    #Extract the instrument name
    if not 'name' in job.keys() or job['name'] == None:
        instrument = util.find_in_dict(header, 'INSTRUMENT_NAME')
        job['name'] = util.instrumentmap[instrument]

    #Check that the required bands are present
    if not util.checkbandnumbers(bands, job['bands']):
        logger.error("Image {} contains bands {}.  Band(s) {} must be present.\n".format(i, bands, job['bands']))
        return

    if 'kerneluri' in job['projection'].keys():
        kernel = job['projection']['kerneluri']
    else:
        kernel = None

    return job

def preprocessimage(job, workingpath, parameters):
    """
    Preprocess a THEMIS EDR for Davinci using ISIS and write files into the
    workingpath.

    Parameters
    ----------
    job : dict
          A dictionary of job containing an image list and
          processing parameters

    workingpath : str
                  The working directory for intermediate files

    Returns
    -------

    outcube: str
               PATH to the preprocessed images.
    """
    # Check that the image exists
    image = job['images']
    if not check_file_exists(image):
        MPI.COMM_WORLD.Abort(1)

    #TODO Check with Trent - is this dumb to offer reprojection?
    #TODO Add reprojection attempt to get ancillary data into the right form.

    logger = logging.getLogger(__name__)
    logger.info('Reading image {}'.format(image))

    # Process the image header
    job = process_header(job)
    if job is None:
        MPI.COMM_WORLD.Abort(1)

    basepath, fname = os.path.split(image)
    fname, _ = os.path.splitext(fname)

    #Convert to ISIS
    outcube = os.path.join(workingpath, '{}.cub'.format(fname))
    kernel = job.get('kernel', None)
    isiswrapper.preprocess_for_davinci(image, outcube, kernel)

    return outcube
