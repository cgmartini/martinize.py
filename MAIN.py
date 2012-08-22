#############
## 8 # MAIN #  -> @MAIN <-
#############
import sys,logging,random,math,os,re
import IO,TOP,DOC,ELN,FUNC

def main(options):
    # Check whether to read from a gro/pdb file or from stdin
    # We use an iterator to wrap around the stream to allow
    # inferring the file type, without consuming lines already
    inStream = IO.streamTag(options["-f"] and options["-f"].value or sys.stdin)
    
    
    # The streamTag iterator first yields the file type, which 
    # is used to specify the function for reading frames
    fileType = inStream.next()
    if fileType == "GRO":
        frameIterator = IO.groFrameIterator
    else:
        frameIterator = IO.pdbFrameIterator
    
    
    ## ITERATE OVER FRAMES IN STRUCTURE FILE ##
    
    # Now iterate over the frames in the stream
    # This should become a StructureFile class with a nice .next method
    model     = 1
    cgOutPDB  = None
    ssTotal   = []
    cysteines = []
    for title,atoms,box in frameIterator(inStream):
    
        if fileType == "PDB":
            # The PDB file can have chains, in which case we list and process them specifically
            # TER statements are interpreted as chain separators
            # A chain may have breaks in which case the breaking residues are flagged
            chains = [ IO.Chain(options,[i for i in IO.residues(chain)]) for chain in IO.pdbChains(atoms) ]        
        else:
            # The GRO file does not define chains. Here breaks in the backbone are
            # interpreted as chain separators. 
            residuelist = [residue for residue in IO.residues(atoms)]
            # The breaks are indices to residues
            broken = breaks(residuelist)
            # Reorder, such that each chain is specified with (i,j,k)
            # where i and j are the start and end of the chain, and 
            # k is a chain identifier
            chains = zip([0]+broken,broken+[len(residuelist)],range(len(broken)+1))
            chains = [ IO.Chain(residuelist[i:j],name=chr(65+k)) for i,j,k in chains ]
    
        for chain in chains:
            chain.multiscale = "all" in options['multi'] or chain.id in options['multi']
    
        # Check the chain identifiers
        if model == 1 and len(chains) != len(set([i.id for i in chains])):
            # Ending down here means that non-consecutive blocks of atoms in the 
            # PDB file have the same chain ID. The warning pertains to PDB files only, 
            # since chains from GRO files get a unique chain identifier assigned.
            logging.warning("Several chains have identical chain identifiers in the PDB file.")
    
        # Check if chains are of mixed type. If so, split them.
        # Note that in some cases HETATM residues are part of a 
        # chain. This will get problematic. But we cannot cover
        # all, probably.
        if not options['MixedChains']:
            demixedChains = []
            for chain in chains:
                demixedChains.extend(chain.split())
            chains = demixedChains
    
        n = 1
        logging.info("Found %d chains:"%len(chains))
        for chain in chains:
            logging.info("  %2d:   %s (%s), %d atoms in %d residues."%(n,chain.id,chain.type(),chain.natoms,len(chain)))
            n += 1
    
        # Check all chains
        keep = []
        for chain in chains:
            if chain.type() == "Water":
                logging.info("Removing %d water molecules (chain %s)."%(len(chain),chain.id))
            elif chain.type() in ("Protein","Nucleic"):
                keep.append(chain)
            elif options['RetainHETATM']:
                keep.append(chain)
            else:
                logging.info("Removing HETATM chain %s consisting of %d residues."%(chain.id,len(chain)))
        chains = keep
    
        
        # Check which chains need merging
        if model == 1:
            order, merge = IO.check_merge(chains, options['mergeList'], options['linkList'], options['CystineCheckBonds'] and options['CystineMaxDist2'])
    
    
        # Get the total length of the sequence
        seqlength = sum([len(chain) for chain in chains])
        logging.info('Total size of the system: %s residues.'%seqlength)
    
    
        ## SECONDARY STRUCTURE
        ss = '' 
        if options['Collagen']:
            for chain in chains:
                chain.set_ss("F")
                ss += chain.ss
        elif options["-ss"]:
            # If the string given for the sequence consists strictly of upper case letters
            # and does not appear to be a file, assume it is the secondary structure
            # Is that safe or silly?
            ss = options["-ss"].value.replace('~','L').replace(' ','L')
            if ss.isalnum() and ss.isupper() and not os.path.exists(options["-ss"].value):
                ss = options["-ss"].value
                logging.info('Secondary structure read from command-line:\n'+ss)
            else:
                # There ought to be a file with the name specified
                ssfile = [ i.strip() for i in open(options["-ss"].value) ]
        
                # Try to read the file as a Gromacs Secondary Structure Dump
                # Those have an integer as first line
                if ssfile[0].isdigit():
                    logging.info('Will read secondary structure from file (assuming Gromacs ssdump).')
                    ss = "".join([ i for i in ssfile[1:] ])
                else:
                    # Get the secondary structure type from DSSP output
                    logging.info('Will read secondary structure from file (assuming DSSP output).')
                    pss = re.compile(r"^([ 0-9]{4}[0-9]){2}")
                    ss  = "".join([i[16] for i in open(options["-ss"].value) if re.match(pss,i)])        
            
            # Now set the secondary structure for each of the chains
            sstmp = ss
            for chain in chains:
                ln = min(len(sstmp),len(chain)) 
                chain.set_ss(sstmp[:ln])
                sstmp = ss[:ln]                         
        else:
            if options["-dssp"]:
                method, executable = "dssp", options["-dssp"].value
            #elif options["-pymol"]:
            #    method, executable = "pymol", options["-pymol"].value
            else:
                logging.warning("No secondary structure or determination method speficied. Protein chains will be set to 'COIL'.")
                method, executable = None, None
        
            for chain in chains:
                ss += chain.dss(method, executable)
        
            # Used to be: if method in ("dssp","pymol"): but pymol is not supported
            if method in ["dssp"]:
                logging.debug('%s determined secondary structure:\n'%method.upper()+ss)
        
        # Collect the secondary structure classifications for different frames
        ssTotal.append(ss)    
    
        # Write the coarse grained structure if requested
        if options["-x"].value:
            logging.info("Writing coarse grained structure.")
            if cgOutPDB == None:
                cgOutPDB = open(options["-x"].value,"w")
            cgOutPDB.write("MODEL %8d\n"%model)
            cgOutPDB.write(title)
            cgOutPDB.write(IO.pdbBoxString(box))
            atid = 1
            for i in order:
                ci = chains[i]
                if ci.multiscale:
                    for r in ci.residues:
                        for name,resn,resi,chain,x,y,z in r:
                            insc  = resi>>20
                            resi -= insc<<20
                            cgOutPDB.write(IO.pdbAtomLine%(atid,name,resn[:3],chain,resi,chr(insc),x,y,z,1,0))
                            atid += 1
                coarseGrained = ci.cg()
                if coarseGrained:
                    for name,resn,resi,chain,x,y,z,ssid in coarseGrained:
                        insc  = resi>>20
                        resi -= insc<<20
                        if ci.multiscale:
                            name = "v"+name
                        cgOutPDB.write(IO.pdbAtomLine%(atid,name,resn[:3],chain,resi,chr(insc),x,y,z,1,ssid))
                        atid += 1 
                    cgOutPDB.write("TER\n")          
                else:
                    logging.warning("No mapping for coarse graining chain %s (%s); chain is skipped."%(ci.id,ci.type()))
            cgOutPDB.write("ENDMDL\n")
    
        # Gather cysteine sulphur coordinates
        cyslist = [cys["SG"] for chain in chains for cys in chain["CYS"]]
        cysteines.append([cys for cys in cyslist if cys])
    
        model += 1
    
    
    # Write the index file if requested
    if options["-n"].value:
        logging.info("Writing index file.")
        NAA,NVZ,NCG = [],[],[]
        atid = 1
        for i in order:
            ci = chains[i]
            coarseGrained = ci.cg()
            if ci.multiscale:
                NAA.extend([" %5d"%(a+atid) for a in range(ci.natoms)]) 
                atid += ci.natoms
            if coarseGrained:
                if ci.multiscale:
                    NVZ.extend([" %5d"%(a+atid) for a in range(len(coarseGrained))])
                else:
                    NCG.extend([" %5d"%(a+atid) for a in range(len(coarseGrained))])
                atid += len(coarseGrained)               
        outNDX   = open(options["-n"].value,"w")
        outNDX.write("\n[ AA ]\n"+"\n".join([" ".join(NAA[i:i+15]) for i in range(0,len(NAA),15)]))
        outNDX.write("\n[ VZ ]\n"+"\n".join([" ".join(NVZ[i:i+15]) for i in range(0,len(NVZ),15)]))
        outNDX.write("\n[ CG ]\n"+"\n".join([" ".join(NCG[i:i+15]) for i in range(0,len(NCG),15)]))
        outNDX.close()
    
    # Evertything below here we only need, if we need to write a Topology
    if options['-o']:

        # Collect the secondary structure stuff and decide what to do with it
        # First rearrange by the residue
        ssTotal = zip(*ssTotal)
        ssAver  = []
        for i in ssTotal:
            si = list(set(i))
            if len(si) == 1:
                # Only one type -- consensus
                ssAver.append(si[0])
            else:
                # Transitions between secondary structure types
                i = list(i)
                si = [(1.0*i.count(j)/len(i),j) for j in si]
                si.sort()
                if si[-1][0] > options["-ssc"].value:
                    ssAver.append(si[-1][1])
                else:
                    ssAver.append(" ")
        
        ssAver = "".join(ssAver)
        logging.info('(Average) Secondary structure has been determined (see head of .itp-file).')
        
        
        # Divide the secondary structure according to the division in chains
        # This will set the secondary structure types to be used for the 
        # topology.
        for chain in chains:
            chain.set_ss(ssAver[:len(chain)])
            ssAver = ssAver[len(chain):]
        
        
        # Now the chains are complete, each consisting of a residuelist, 
        # and a secondary structure designation if the chain is of type 'Protein'.
        # There may be mixed chains, there may be HETATM things. 
        # Water has been discarded. Maybe this has to be changed at some point.
        # The order in the coarse grained files matches the order in the set of chains.
        #
        # If there are no merges to be done, i.e. no global Elnedyn network, no 
        # disulphide bridges, no links, no distance restraints and no explicit merges,
        # then we can write out the topology, which will match the coarse grained file.
        #
        # If there are merges to be done, the order of things may be changed, in which
        # case the coarse grained structure will not match with the topology...
        
        
        ## CYSTINE BRIDGES ##
        
        # Extract the cysteine coordinates (for all frames) and the cysteine identifiers
        if options['CystineCheckBonds']:
            logging.info("Checking for cystine bridges, based on sulphur (SG) atoms lying closer than %.4f nm"%math.sqrt(options['CystineMaxDist2']/100))
        
            cyscoord  = zip(*[[j[4:7] for j in i] for i in cysteines])
            cysteines = [i[:4] for i in cysteines[0]]
        
            bl, kb    = options['ForceField'].special[(("SC1","CYS"),("SC1","CYS"))]
        
            # Check the distances and add the cysteines to the link list if the 
            # SG atoms have a distance smaller than the cutoff.
            rlc = range(len(cysteines))
            for i in rlc[:-1]:
                for j in rlc[i+1:]:
                    # Checking the minimum distance over all frames
                    # But we could also take the maximum, or the mean
                    d2 = min([FUNC.distance2(a,b) for a,b in zip(cyscoord[i],cyscoord[j])])
                    if d2 <= options['CystineMaxDist2']:
                        a, b = cysteines[i], cysteines[j]
                        options['linkListCG'].append((("SC1","CYS",a[2],a[3]),("SC1","CYS",b[2]-(32<<20),b[3]),bl,kb))
                        a,b = (a[0],a[1],a[2]-(32<<20),a[3]),(b[0],b[1],b[2]-(32<<20),b[3])
                        logging.info("Detected SS bridge between %s and %s (%f nm)"%(a,b,math.sqrt(d2)/10))
        
        
        ## REAL ITP STUFF ##
        
        # Check whether we have identical chains, in which case we 
        # only write the ITP for one...
        # This means making a distinction between chains and 
        # moleculetypes.
        
        molecules = [tuple([chains[i] for i in j]) for j in merge]
        
        # At this point we should have a list or dictionary of chains
        # Each chain should be given a unique name, based on the value
        # of options["-o"] combined with the chain identifier and possibly
        # a number if there are chains with identical identifiers.
        # For each chain we then write an ITP file using the name for 
        # moleculetype and name + ".itp" for the topology include file.
        # In addition we write a master topology file, using the value of
        # options["-o"], with an added extension ".top" if not given.
        
        # *NOTE*: This should probably be gathered in a 'Universe' class
        itp = 0
        moleculeTypes = {}
        for mi in range(len(molecules)):
            mol = molecules[mi]
            # Check if the moleculetype is already listed
            # If not, generate the topology from the chain definition
            if not mol in moleculeTypes or options['SeparateTop']:
                # Name of the moleculetype
                # NOTE: The naming should be changed; now it becomes Protein_X+Protein_Y+...
                name = "+".join([chain.getname(options['-name'].value) for chain in mol])
                moleculeTypes[mol] = name
    
                # Write the molecule type topology
                top = TOP.Topology(mol[0],options=options,name=name)
                for m in mol[1:]:
                    top += TOP.Topology(m,options=options)
    
                # Have to add the connections, like the connecting network
                # Gather coordinates
                mcg, coords = zip(*[(j[:4],j[4:7]) for m in mol for j in m.cg()])
                mcg         = list(mcg)
        
                # Run through the link list and add connections
                for atomA,atomB,bondlength,forceconst in options['linkListCG']:
                    if bondlength == -1 and forceconst == -1:
                        bondlength, forceconst = options['ForceField'].special[(atomA[:2],atomB[:2])]
                    # Check whether this link applies to this group
                    atomA = atomA in mcg and mcg.index(atomA)+1
                    atomB = atomB in mcg and mcg.index(atomB)+1
                    if atomA and atomB:
                        cat = (mcg[atomA][1] == "CYS" and mcg[atomB][1] == "CYS") and "Cystine" or "Link"
                        top.bonds.append(TOP.Bond((atomA,atomB),options=options,type=1,parameters=(bondlength,forceconst),category=cat))
        
                # Elastic Network
                # The elastic network is added after the topology is constructed, since that
                # is where the correct atom list with numbering and the full set of 
                # coordinates for the merged chains are available. 
                if options['ElasticNetwork']:
                    rubberType = options['ForceField'].EBondType
                    rubberList = ELN.rubberBands(
                        [(i[0],j) for i,j in zip(top.atoms,coords) if i[4] in options['ElasticBeads']],
                        options['ElasticLowerBound'],options['ElasticUpperBound'],
                        options['ElasticDecayFactor'],options['ElasticDecayPower'],
                        options['ElasticMaximumForce'],options['ElasticMinimumForce'])
                    top.bonds.extend([TOP.Bond(i,options=options,type=rubberType,category="Rubber band") for i in rubberList])
    
                # Write out the MoleculeType topology
                destination = options["-o"] and open(moleculeTypes[mol]+".itp",'w') or sys.stdout
                destination.write(str(top))        
        
                itp += 1
        
            # Check whether other chains are equal to this one 
            # Skip this step if we are to write all chains to separate moleculetypes
            if not options['SeparateTop']:
                for j in range(mi+1,len(molecules)):
                    if not molecules[j] in moleculeTypes and mol == molecules[j]:
                        # Molecule j is equal to a molecule mi
                        # Set the name of the moleculetype to the one of that molecule
                        moleculeTypes[molecules[j]] = moleculeTypes[mol]
        
        logging.info('Written %d ITP file%s'%(itp,itp>1 and "s" or ""))
                
        
        
        #----+--------------------------------------------
        ## B # MORE WORK -- WRITING THE MASTER TOPOLOGY ##
        #----+--------------------------------------------
        
        
        # Processing stuff
        
        # Output stream
        top  = options["-o"] and open(options['-o'].value,'w') or sys.stdout
        
        # ITP file listing
        itps = '\n'.join(['#include "%s.itp"'%molecule for molecule in set(moleculeTypes.values())])
        
        # Molecule listing
        logging.info("Output contains %d molecules:"%len(molecules))
        n = 1
        for molecule in molecules:
            stuff = (n, moleculeTypes[molecule], len(molecule)>1 and "s" or " ", " ".join([i.id for i in molecule]))
            logging.info("  %2d->  %s (chain%s %s)"%stuff)
            n += 1
        molecules   = '\n'.join(['%s \t 1'%moleculeTypes[molecule] for molecule in molecules])
        
        # Set a define if we are to use rubber bands
        useRubber   = options['ElasticNetwork'] and "#define RUBBER_BANDS" or ""
        
        # Do not set a define for position restrains here, as people are more used to do it in mdp file?
        
        top.write(
'''#include "martini.itp"
    
%s
  
%s
    
[ system ]
; name
Martini system from %s
    
[ molecules ]
; name        number
%s
        ''' % (useRubber, itps, options["-f"] and options["-f"].value or "stdin", molecules))
    
        logging.info('Written topology files')
    
    # Maybe there are forcefield specific log messages?
    options['ForceField'].messages()

    # The following lines are always printed (if no errors occur).
    print "\n\tThere you are. One MARTINI. Shaken, not stirred.\n"
    
    Q = DOC.martiniq.pop(random.randint(0,len(DOC.martiniq)-1))
    print "\n", Q[1], "\n%80s"%("--"+Q[0]), "\n"
    
    
