I'm placing data here which can be used by collaborators to investigate the phase-lag in the periodic variables.

Data is stored in `hdf5` format (.h5), readily accessible via tools such as astropy.table.

Column descriptions:

SOURCEID: unique identifier
MEANMJDOBS: the Modified Julian Date (MJD) of the observation
JAPERMAG3: the 3 arcsecond aperture magnitude of the source at J band (and so on for H, K).
JPPERRBITS: Quality/error flags for photometry (see http://wsa.roe.ac.uk/ppErrBits.html)
j_corr: Residual J magnitude after a 0, 2, or 4th degree polynomial was subtracted during the periodicity search.
