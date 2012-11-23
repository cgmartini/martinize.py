##############################
## 2 # COMMAND LINE PARSING ##  -> @CMD <-
##############################
import sys,logging
import DOC

# Helper function to parse atom strings given on the command line:
#   resid
#   resname/resid
#   chain/resname/resid
#   resname/resid/atom
#   chain/resname/resid/atom
#   chain//resid
#   chain/resname/atom
def str2atom(a):
    a = a.split("/")   
    if len(a) == 1: # Only a residue number:
        return (None,None,int(a[0]),None)
    if len(a) == 2: # Residue name and number (CYS/123):
        return (None,a[0],int(a[1]),None)
    if len(a) == 3:
        if a[2].isdigit(): # Chain, residue name, residue number
            return (None,a[1],int(a[2]),a[0])
        else: # Residue name, residue number, atom name
            return (a[2],a[0],int(a[1]),None)
    return (a[3],a[1],int(a[2]),a[0])

def option_parser(args,options,lists,version=0):

    # Check whether there is a request for help
    if '-h' in args or '--help' in args:
        DOC.help()

    # Convert the option list to a dictionary, discarding all comments
    options = dict([i for i in options if not type(i) == str])

    # This information we would like to print to some files, so let's put it in our information class
    options['Version']             = version
    options['Arguments']           = args[:]
   
    while args:
        ar = args.pop(0)
        options[ar].setvalue([args.pop(0) for i in range(options[ar].num)])

    ## LOGGING ##
    # Set the log level and communicate which options are set and what is happening
    # If 'Verbose' is set, change the logger level
    logLevel = options["-v"] and logging.DEBUG or logging.INFO
    logging.basicConfig(format='%(levelname)-7s    %(message)s',level=logLevel)

    logging.info('MARTINIZE, script version %s'%version)
    logging.info('If you use this script please cite:')
    logging.info('de Jong et al., J. Chem. Theory Comput., 2013, DOI:10.1021/ct300646g')
    
    # Process the raw options from the command line
    # Boolean options are set to more intuitive variables
    options['Collagen']            = options['-collagen']
    options['ChargesAtBreaks']     = options['-cb']
    options['NeutralTermini']      = options['-nt']
    options['ExtendedDihedrals']   = options['-ed']
    options['RetainHETATM']        = False # options['-hetatm']
    options['SeparateTop']         = options['-sep']
    options['MixedChains']         = False # options['-mixed']
    options['ElasticNetwork']      = options['-elastic']
    
    # The make the program flexible, the forcefield parameters are defined
    # for multiple forcefield. Check if a existing one is defined:
    ###_tmp  = __import__(options['-ff'].value.lower())
    ###options['ForceField']  = getattr(_tmp,options['-ff'].value.lower())()
    try:
        try:
            # Try to load the forcefield class from a different file
            _tmp  = __import__(options['-ff'].value.lower())
            options['ForceField']  = getattr(_tmp,options['-ff'].value.lower())()
        except:
            # Try to load the forcefield class from the current file
            options['ForceField']  = globals()[options['-ff'].value.lower()]()
    except:
        logging.error("Forcefield '%s' can not be found."%(options['-ff']))
        sys.exit()

 
    # Parsing of some other options into variables
    options['ElasticMaximumForce'] = options['-ef'].value 
    options['ElasticMinimumForce'] = options['-em'].value
    options['ElasticLowerBound']   = options['-el'].value
    options['ElasticUpperBound']   = options['-eu'].value
    options['ElasticDecayFactor']  = options['-ea'].value
    options['ElasticDecayPower']   = options['-ep'].value
    options['ElasticBeads']        = options['-eb'].value.split(',')
    options['PosResForce']         = options['-pf'].value

    options['PosRes']              = [i.lower() for i in options['-p'].value.split(",")]
    if "none"     in options['PosRes']: options['PosRes'] = []
    if "backbone" in options['PosRes']: options['PosRes'].append("BB")
    
    
    if options['ForceField'].ElasticNetwork:
        #
        # Some forcefields, like elnedyn, always use an elatic network. This is set in the 
        # forcefield file, with the parameter ElasticNetwork.
        #
        options['ElasticNetwork']  = True
#        # Unless explicitly told not to, using Elnedyn will 
#        # merge all chains into a single moleculetype to 
#        # allow a global Elnedyn network.
#        if "no" in lists['merges']:
#            lists['merges'].remove("no")
#        else:
#            lists['merges'] = ["all"]
    
    
    # Merges, links and cystines
    options['mergeList'] = "all" in lists['merges'] and ["all"] or [i.split(",") for i in lists['merges']]

    # Process links
    linkList   = []
    linkListCG = []
    for i in lists['links']:
        ln     = i.split(",")
        a, b   = str2atom(ln[0]), str2atom(ln[1])
        if len(ln) > 3: # Bond with given length and force constant
            bl, fc = (ln[2] and float(ln[2]) or None, float(ln[3]))
        elif len(a) == 3: # Constraint at given distance
            bl, fc = float(a[2]), None
        else: # Constraint at distance in structure
            bl, fc = None, None
        # Store the link, but do not list the atom name in the
        # atomistic link list. Otherwise it will not get noticed 
        # as a valid link when checking for merging chains
        linkList.append(((None,a[1],a[2],a[3]),(None,b[1],b[2],b[3])))
        linkListCG.append((a,b,bl,fc))
    
    
    # Cystines -- This should be done for all special bonds listed in the _special_ dictionary
    CystineCheckBonds = False   # By default, do not detect cystine bridges
    CystineMaxDist2   = (10*0.22)**2 # Maximum distance (A) for detection of SS bonds
    for i in lists['cystines']:
        if i.lower() == "auto":
            CystineCheckBonds = True
        elif i.replace(".","").isdigit():
            CystineCheckBonds = True
            CystineMaxDist2   = (10*float(i))**2
        else:
            # This item should be a pair of cysteines
            cysA, cysB = [str2atom(j) for j in i.split(",")]
            # Internally we handle the residue number shifted by ord(' ')<<20. We have to add this to the
            # cys-residue numbers given here as well.
            constant = 32<<20
            linkList.append((("SG","CYS",cysA[2]+constant,cysA[3]),("SG","CYS",cysB[2]+constant,cysB[3])))
            linkListCG.append((("SC1","CYS",cysA[2]+constant,cysA[3]),("SC1","CYS",cysB[2]+constant,cysB[3]),-1,-1))

    # Now we have done everything to it, we can add Link/cystine related stuff to options
    # 'multi' is not stored anywhere else, so that we also add
    options['linkList']          = linkList
    options['linkListCG']        = linkListCG
    options['CystineCheckBonds'] = CystineCheckBonds
    options['CystineMaxDist2']   = CystineMaxDist2
    options['multi']             = lists['multi']

    logging.info("Chain termini will%s be charged"%(options['NeutralTermini'] and " not" or ""))
    
    logging.info("Residues at chain brakes will%s be charged"%((not options['ChargesAtBreaks']) and " not" or ""))

    if options.has_key('ForceField'):
        logging.info("The %s forcefield will be used."%(options['ForceField'].name))
    else:
        logging.error("Forcefield '%s' has not been implemented."%(options['-ff']))
        sys.exit()
    
    if options['ExtendedDihedrals']:  
        logging.info('Dihedrals will be used for extended regions. (Elastic bonds may be more stable)')
    else:                  
        logging.info('Local elastic bonds will be used for extended regions.')
    
    
    if options['PosRes']:
        logging.info("Position restraints will be generated.")
        logging.warning("Position restraints are only enabled if -DPOSRES is set in the MDP file")
    
    
    if options['MixedChains']:
        logging.warning("So far no parameters for mixed chains are available. This might crash the program!")
    
    
    if options['RetainHETATM']:
        logging.warning("I don't know how to handle HETATMs. This will probably crash the program.")

    return options 
