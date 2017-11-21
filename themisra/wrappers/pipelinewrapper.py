import subprocess
import os
import logging

from pysis.exceptions import ProcessError

logger = logging.getLogger('ThemisTi')


def themis_davinci(inpath, outpath, deplaid=False, uddw=False, tesatm=False, rtilt=False, force=False):
    """
    Run the images through a processing pipeline to convert to temperature.
    Parameters
    ----------
    inpath : str
             PATH to the image to be processed
    outpath : str
              PATH to write the resulting Davinci processed cube
    deplaid : bool
              Determines whether or not to apply deplaid
    uddw : bool
           Determines whether or not to apply uddw
    tesatm : bool
             Determines whether or not to apply tesatm
    rtilt : bool
            Determines whether or not to apply rtilt
    force : bool
            Determines whether or not to apply force
    Returns
    -------
    outpath : str
              PATH to the Davinci processed cube
    """
    basepath = os.path.dirname(__file__)
    davincipath = '../davinci/'
    fpath = os.path.join(basepath, davincipath)
    fpath = os.path.join(fpath, 'ti_pipeline.dv')
    basepath, fname = os.path.split(inpath)
    outname, ext = os.path.splitext(fname)
    outpath = os.path.join(outpath, outname) + '_dvprocessed'
    cmd = r'{} {} {} {} {} {} {} {}'.format(fpath, int(uddw), int(tesatm),
                                            int(deplaid), int(rtilt), int(force),
                                            inpath, outpath)
    cmd = cmd.split()
    try:
        response = subprocess.check_output(cmd, shell=False)
        #logger.debug(response)
    except ProcessError as e:
        logger.error(e.stdout)
        logger.error(e.stderr)
    return outpath
