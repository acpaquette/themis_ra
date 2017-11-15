import os

from plio.utils.utils import find_in_dict
import pvl
import pysis
from pysis import isis


def preprocess_for_davinci(image, outcube, kernel):
    """
    Run spiceinit and campt to determine incidence and local time.  These
    parameters are required by the davinci pipeline in order to
    determine whether or not deplaid should be run.
    """
    isis.thm2isis(from_=image, to=outcube)
    if kernel is not None:
        print(kernel)
        isis.spiceinit(from_=outcube, ck=kernel)
    else:
        isis.spiceinit(from_=outcube)

def postprocess_for_davinci(incube, kernel=None):
    workingpath, fname = os.path.split(incube)
    fname = os.path.splitext(fname)[0]
    #Processing the temperature to a level2 image
    if kernel:
        isis.spiceinit(from_=incube, ck=kernel)
    else:
        isis.spiceinit(from_=incube)

    isiscube = os.path.join(workingpath, '{}_proj.cub'.format(fname))
    isis.cam2map(from_=incube, to=isiscube,
                 map='$base/templates/maps/simplecylindrical.map')

    return isiscube

def campt_header(outcube):
    """
    Compute the incidence angle at the center of the image and the local
    solar time.  These are required by the Davinci processing pipeline to
    determine what processing to perform.
    """

    workingpath, fname = os.path.split(outcube)
    fname = os.path.splitext(fname)[0]

    header = pvl.load(outcube)
    samples = find_in_dict(header, 'Samples')
    lines = find_in_dict(header, 'Lines')

    coordinatelist = os.path.join(workingpath, 'coordinatelist.lis')
    with open(coordinatelist, 'w') as f:
        f.write('{},{}\n'.format(samples/2, lines/2))
        f.write('1,1\n') #UpperLeft
        f.write('{},{}\n'.format(samples-1, lines-1)) #LowerRight
    campt = pvl.loads(isis.campt(from_=outcube, to=os.path.join(workingpath, fname + '_campt.pvl'),
                      usecoordlist='yes', coordlist=coordinatelist,
                      coordtype='image'))
    for j, g in enumerate(campt.items()):
        if j == 0:
            #Incidence at the center of the image
            try:
                incidence = g[1]['Incidence'].value
            except:
                incidence = g[1]['Incidence']
        elif j == 1:
            #Upper Left Corner Pixel
            stoplocaltime = g[1]['LocalSolarTime'].value
        elif j == 2:
            #Lower Right Corner Pixel
            startlocaltime = g[1]['LocalSolarTime'].value
    return incidence, stoplocaltime, startlocaltime
