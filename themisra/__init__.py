__name__ = 'themis_ra'
__version__ = '0.1.1'
__minimum_isis_version__ = (3,4)

# Hard dependency on ISIS3
import pysis
pysis.check_isis_version(*__minimum_isis_version__)
