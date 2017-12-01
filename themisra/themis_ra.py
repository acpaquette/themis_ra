import glob
import logging
import os
import subprocess
import sys
import time
import pvl

from mpi4py import MPI

from plio.io import io_gdal, io_hdf, io_json
from plio.date import astrodate, julian2ls, julian2season
import plio.utils
from plio.utils import log
from plio.utils.utils import check_file_exists, find_in_dict
import themisra.utils.utils as util

import themisra.processing.processing as processing

from themisra.wrappers import pipelinewrapper, isiswrapper

#Constants
instrumentmap = {'THERMAL EMISSION IMAGING SYSTEM':'THEMIS'}  #Mapping of instrument names as stored in the header to short names


#Get MPI to abort cleanly in the case of an error
sys_excepthook = sys.excepthook
def mpi_excepthook(v, t, tb):
    sys_excepthook(v, t, tb)
    MPI.COMM_WORLD.Abort(1)
sys.excepthook = mpi_excepthook

    #Setup logging
#log.setup_logging(level=config.LOG_LEVEL)
logger = logging.getLogger(__name__)

def cost(temp):
    print()

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
        working_path = plio.utils.utils.create_dir(basedir=job['workingdir'])

        # ISIS preprocessing
        processing.preprocess_image(job, working_path)

        # DaVinci processing
        isistemp, isisrad = processing.process_image(job,working_path)
        processing.map_ancillary(isistemp, job)

        band_three = util.extract_band(job, isistemp, 3)
        band_nine = util.extract_band(job, isistemp, 9)

        rock_three = util.generate_rad_image(band_three, 3)
        rock_nine = util.generate_rad_image(band_nine, 9)

if __name__ == '__main__':
    main()
