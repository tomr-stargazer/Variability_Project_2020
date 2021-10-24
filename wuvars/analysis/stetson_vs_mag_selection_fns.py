#######################
# Selection functions #
#######################

# These are here 

def onc_selection_fn(ds):
    """
    Selects a vertical strip to exclude from the center of the ONC field, i.e.:

    xoox
    xoox
    xoox
    xoox

    x: used
    o: not used

    """    
    
    ra = ds['mean']['RA']
    dec = ds['mean']['DEC']

    ra_range = max(ra) - min(ra)
    selection_region = (ra > (min(ra) + ra_range/4)) & (ra < (max(ra) - ra_range/4))
    
    return selection_region


def ngc_selection_fn(ds):
    """
    Selects the central 25% box to exclude from the NGC field, i.e.:

    xxxx
    xoox
    xoox
    xxxx

    x: used
    o: not used

    """    
    
    ra = ds['mean']['RA']
    dec = ds['mean']['DEC']

    ra_range = max(ra) - min(ra)
    dec_range = max(dec) - min(dec)
    selection_region = (
        (ra > (min(ra) + ra_range/4)) & (ra < (max(ra) - ra_range/4)) &
        (dec > (min(dec) + dec_range/4)) & (dec < (max(dec) - dec_range/4))
                       )
    return selection_region


def ic_selection_fn(ds):
    """
    Selects the southwest 50% box to exclude from the IC field, i.e.:

    xxxx
    xooo
    xooo
    xooo

    x: used
    o: not used

    """        
    ra = ds['mean']['RA']
    dec = ds['mean']['DEC']

    ra_range = max(ra) - min(ra)
    dec_range = max(dec) - min(dec)
    selection_region = (
        (ra > min(ra)) & (ra < (max(ra) - ra_range/4)) &
        (dec > min(dec)) & (dec < (max(dec) - dec_range/4))
        )
                       
    return selection_region
