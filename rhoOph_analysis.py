""" This module performs a likelihood analysis of diffuse sources in the 
    Rho Ophiuchus region. It is sensitive to inputs (directories, etc.) so
    cannot easily be ported elsewhere, without needing to specify the input
    files. 
    To call on the command line:
    > ipython rhoOph_analysis key
    where key = "Tau353", "NICER", "NICEST", "DobDSS", or "Dob2MASS"
"""
def MakeLike(predir,direc,fit='NewMinuit'):
    """ Make observation and likelihood objects. 
        Use: (obs,like) = MakeLike(predir,direc,optimizer)
             optimizer = "DRMNFB", "Minuit", or "NewMinuit"
             predir = working directory of Rho Oph analysis
             direc = directory for diffuse source results/sourcemaps
        Returns: obs = observations object. Used for later fitting.
                 like = likelihood object. Used for later fitting.
    """
    if ( direc == 'fermi' ):
        mdl = predir + '/MODEL_3FGL.xml'
    else:
        mdl = predir + direc + '/mdl_%s.xml'%direc.lower()

    bexpmap = predir + 'bexpmap.fits'
    if ( direc == 'fermi' ):
        srcmap = 'srcmap_fermi.fits'
    else:
        srcmap = predir + '/srcmap_allDiff_oldIC.fits'

    # observation inputs, con't
    expcube = predir + 'ltCube.fits'
    irfs = 'P8R2_ULTRACLEANVETO_V6'
    fit_good = ['newminuit','drmnfb','minuit']

    obs = BinnedObs(binnedExpMap=bexpmap, srcMaps=srcmap, expCube=expcube, 
                    irfs=irfs)
    if ( fit.lower() not in fit_good ):
        raise KeyError("Improper optimizer. Use one of: %s"%fit_good)
    else:
        like = BinnedAnalysis(obs,srcModel=mdl, optimizer=fit)

    # For analyses NOT using Fermi diffuse model: change isotropic emission.
    if ( direc != 'fermi' ):
        eg_v05 = like.model['eg_v05'].funcs['Spectrum']
        try:
            eg_v05.params['Normalization'].setBounds([0.2,1.3])
        except:
            eg_v05.params['Normalization'].setBound([0.2,1.3])
            print ("Incorrect command -- ...setBound(), not ...setBounds()")
        else:
            pass

    return (obs,like)

def LikeFit(like,fit='NewMinuit'):
    """ Perform likelihood fit with chosen optimizer.
        Use: (log like, like object) = LikeFit(optimizer)
             optimizer = "DRMNFB", "Minuit", or "NewMinuit"
        Returns: log likelihood value
                 likelihood object, used for further analysis
    """
    # Fit model to gamma-ray data with covariance, if applicable
    if ( fit.lower() == 'newminuit' ):
        likeobj = pyLike.NewMinuit(like.logLike)
        lglike = like.fit(covar=True,optObject=likeobj)
    elif ( fit.lower() == 'minuit' ):
        likeobj = pyLike.Minuit(like.logLike)
        lglike = like.fit(covar=True,optObject=likeobj)
    elif ( fit.lower() == 'drmnfb' ):
        lglike = like.fit()
    else:
        raise KeyError("Improper optimizer. Use 'NewMinuit', 'Minuit', or"
                       + "'DRMNFB'")
    return (lglike,like)

def IterateFit(like,PS,diff_src):
    """ Iterative fit: PS fit with DRMNFB, then diffuse fit with DRMNFB, then
        diffuse with NewMinuit for more accurate fit + covariance matrix.
        Do this twice over: PS + diff, then another PS + diffuse.
        Use: IterateFit(likelihood object, PS names, diffuse source names)
             PS names, diffuse source names = list
        Returns: likelihood object, log(likelihood) of final fit
    """
    # Free point sources and fit.
    like = ParamChange(like,PS,change='thaw')
    print "Fitting %s point source normalizations\n"%len(PS)
    lglike1_ps,like = LikeFit(like,fit='DRMNFB')
    
    # Freeze point sources, free diffuse sources, and fit
    like = ParamChange(like,diff_src,change='thaw')
    print "Fitting diffuse sources: %s"%diff_src
    lglike1_diff,like = LikeFit(like,fit='DRMNFB')

    # Get more accurate fit on diffuse sources.
    print "NewMinuit fit on diffuse sources"
    lglike_diff_newmin, like = LikeFit(like,fit="NewMinuit")
    like = ParamChange(like,PS,change='thaw')
    print "Fitting point sources again"
    lglike_ps2, like = LikeFit(like,fit='DRMNFB')

    # Final fit: diffuse sources only.
    like = ParamChange(like,diff_src,change='thaw')
    lglike_diff_newmin, like = LikeFit(like,fit='NewMinuit')
    return (like,lglike_diff_newmin)

def ParamChange(like,sr,change='freeze'):
    """ Free up point source normalizations in order to fit them. 
        Use: like = PSchange(like,sr,change='freeze')
        > change = 'freeze' OR 'thaw'
        > sr = PS, or diff_src
               This specifies which sources we want to change,
               point source or diffuse sources.
    """
    if (len(sr) == 0):
        raise ValueError("No source specified")

    # Default freeze everything.
    for i in xrange( len(like.model.params) ):
        like.freeze(i)
    
    # Free up specified sources ("sr")
    if (change.lower() == 'thaw'):
        for src in sr:
            pars = like.model[src].funcs['Spectrum'].paramNames
            par_ind = min(map(like.par_index,[src]*len(pars),pars))
            like.thaw(par_ind)
    return like


def main(A):
    import time
    start = time.clock()
    ## Globals ##
    predir = ('/home/abrahams/HICO_survey/SourceSearch/RhoOph/l353b17/' 
              + 'ro_pass8_files/CleanData/')
    diff_src = ['HI','bubble','eg_v05','galprop']
    # Nearby or 'strong' or highly curved point sources
    PS = ['3FGL J1621.1-2331', '3FGL J1628.0-3203', '3FGL J1614.5-2231', 
          '3FGL J1617.3-2519', '3FGL J1714.6-3327', '3FGL J1616.8-2300',
          '3FGL J1626.0-2951', '3FGL J1625.7-2527', '3FGL J1553.3-2421',
          '3FGL J1645.7-2149', '3FGL J1625.6-2058']

    ## Make likelihood objects: define directories, create obs, like objects
    if ( A.lower() == 'tau353' ):
        direc = 'Tau353'
        diff_name = 'Tau353'
    elif ( A.lower() == 'nicer' ):
        direc = 'NICER'
        diff_name = 'NICER'
    elif ( A.lower() == 'nicest' ):
        direc = 'NICEST'
        diff_name = 'NICEST'
    elif ( A.lower() == 'dobdss' ):
        direc = 'DobDSS'
        diff_name = 'Dobashi DSS'
    elif ( A.lower() == 'dob2mass' ):
        direc = 'Dob2MASS'
        diff_name = 'Dobashi 2MASS'
    elif ( A.lower() == 'fermi' ):
        direc = 'fermi'
        diff_name = 'fermi'
    else:
        raise KeyError("Incorrect analysis source. Call:\n"
           + "ipython rhoOph_analysis key\n"
           + "where key = nicest, dobdss, dob2mass, tau353, or nicer")

    # Fermi diffuse model is does NOT have own directory.
    diff_src.append(diff_name)
    if ( direc == 'fermi' ):
        f = open('LogLike.dat','a')
        g = open('README','a')
    else:
        f = open(direc+'/LogLike.dat','a')
        g = open(direc+'/README','a')

    # Round 1: fit ALL diffuse sources and write results.
    g.write("Starting oldIC analysis of %s\n"%diff_name)
    obs,like = MakeLike(predir,direc,fit='DRMNFB')
    like,lglike1 = IterateFit(like,PS,diff_src)
    f.write('Log(like) for ALL = %s\n'%round(lglike1,4))
    if ( direc == 'fermi' ):
        like.writeXml('mdl_all.xml')
    else:
        like.writeXml(direc+'/mdl_all.xml')

    # Round 2: remove Fermi Bubbles (from like AND diff_src) and fit again 
    like.deleteSource('bubble')
    g.write("Starting oldIC analysis of %s ... removing bubbles\n"%diff_name)
    diff_src.pop( diff_src.index('bubble') )
    (like,lglike2) = IterateFit(like,PS,diff_src)
    f.write('Log(like) for NO bubble = %s\n'%round(lglike2,4))
    if ( direc == 'fermi' ):
        like.writeXml('mdl_noBubble.xml')
    else:
        like.writeXml(direc+'/mdl_noBubble.xml')

    # Round 3: remove HI, too.
    like.deleteSource('HI')
    g.write("Starting oldIC analysis of %s .... removing HI\n"%diff_name)
    diff_src.pop( diff_src.index('HI') )
    (like,lglike3) = IterateFit(like,PS,diff_src)
    f.write('Log(like) for NO bubble OR HI = %s\n'%round(lglike3,4))
    if ( direc == 'fermi' ):
        like.writeXml('mdl_noBubble_noHI.xml')
    else:
        like.writeXml(direc+'/mdl_noBubble_noHI.xml')

    g.write("Finished oldIC analysis of %s\n"%diff_name)
    g.write("Timing: %s"%(time.clock()-start))

import pyLikelihood
from BinnedAnalysis import *
import sys

#
# Fit the following dust tracer.
#
gas = ['DobDSS','Dob2MASS']
#gas = ['Tau353', 'NICER', 'NICEST', 'DobDSS', 'Dob2MASS']
for template in gas:
    print "Starting %s"%template
    main(template)
    print "Done %s"%template

print "Done"
