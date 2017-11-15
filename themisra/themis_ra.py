import glob
import logging
import os
#import re
#import subprocess
import sys
import time

#import numpy as np
from mpi4py import MPI

from plio.io import io_gdal, io_hdf, io_json
from plio.date import astrodate, julian2ls, julian2season
import plio.utils
from plio.utils import log
from plio.utils.utils import check_file_exists, find_in_dict

import themisra.utils.utils as util

from themisra.wrappers import pipelinewrapper, isiswrapper
#import gdal

#import pysis
#import pvl

#from krc.wrappers import pipelinewrapper, isiswrapper
#from krc.utils import utils
#from krc.interpolation import interpolator as interp
#from krc import config

#Constants
instrumentmap = {'THERMAL EMISSION IMAGING SYSTEM':'THEMIS'}  #Mapping of instrument names as stored in the header to short names
processingpipelines = {'themis_davinci':pipelinewrapper.themis_davinci}


#Get MPI to abort cleanly in the case of an error
#sys_excepthook = sys.excepthook
def mpi_excepthook(v, t, tb):
    sys_excepthook(v, t, tb)
    MPI.COMM_WORLD.Abort(1)
sys.excepthook = mpi_excepthook

    #Setup logging
#log.setup_logging(level=config.LOG_LEVEL)
logger = logging.getLogger(__name__)

def main():
    comm = MPI.COMM_WORLD
    rank = comm.rank
    if rank == 0:
        t_start = time.time()
        #Parse the job input
        if len(sys.argv) < 2:
            logger.error("Please supply an input configuration file.")
            sys.exit()
        logger.info("Processing using {} cores".format(comm.size))
        job = io_json.read_json(sys.argv[1])

        #Create a temporary working directory
        workingpath = plio.utils.utils.create_dir(basedir=job['workingdir'])
        # Storage for image / observation parameters
        parameters = {}

        # ISIS preprocessing
        util.preprocessimage(job, workingpath, job['images'])


if __name__ == '__main__':
    main()

