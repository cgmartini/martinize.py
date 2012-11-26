##################
## 7 # TOPOLOGY ##  -> @TOP <-
##################
import IO,FUNC,MAP
import logging,math

# This is a generic class for Topology Bonded Type definitions
class Bonded:
    # The init method is generic to the bonded types,
    # but may call the set method if atoms are given
    # as (ID, ResidueName, SecondaryStructure) tuples
    # The set method is specific to the different types.
    def __init__(self,other=None,options=None,**kwargs):
        self.atoms = []
        self.type = -1
        self.parameters = []
        self.comments = []
        self.category = None 

        if options and type(options) == dict:
            self.options = options
        if other:
            # If other is given, then copy the attributes
            # if it is of the same class or set the 
            # attributes according to the key names if
            # it is a dictionary
            if other.__class__ == self.__class__:
                for attr in dir(other):
                    if not attr[0] == "_":
                        setattr(self,attr,getattr(other,attr))
            elif type(other) == dict:
                for attr in other.keys():
                    setattr(self,attr,other[attr])
            elif type(other) in (list,tuple):
                self.atoms = other

        # For every item in the kwargs keys, set the attribute
        # with the same name. This can be used to specify the 
        # attributes directly or to override attributes 
        # copied from the 'other' argument.
        for key in kwargs:
            setattr(self,key,kwargs[key])
        
        # If atoms are given as tuples of
        # (ID, ResidueName[, SecondaryStructure])
        # then determine the corresponding parameters 
        # from the lists above
        if self.atoms and type(self.atoms[0]) == tuple:
            self.set(self.atoms,**kwargs)          

    def __nonzero__(self):
        return bool(self.atoms) 

    def __str__(self):
        if not self.atoms or not self.parameters:
            return ""
        s = ["%5d" % i for i in self.atoms]
        # For exclusions, no type is defined, which equals -1
        if self.type != -1: s.append(" %5d " % self.type)
        # Print integers and floats in proper format and neglect None terms
        s.extend([FUNC.formatString(i) for i in self.parameters if i != None])
        if self.comments:
            s.append(';')
            if type(self.comments) == str:
                s.append(self.comments)
            else:
                s.extend([str(i) for i in self.comments])
        return " ".join(s)

    def __iadd__(self,num):
        self.atoms = [i+int(num) for i in self.atoms]
        return self

    def __add__(self,num):
        out  = self.__class__(self)
        out += num
        return out

    def __eq__(self,other):
        if type(other) in (list,tuple):
            return self.atoms == other
        else:
            return self.atoms == other.atoms and self.type == other.type and self.parameters == other.parameters

    # This function needs to be overridden for descendents
    def set(self,atoms,**kwargs):
        pass


# The set method of this class will look up parameters for backbone beads
# Side chain bonds ought to be set directly, using the constructor
# providing atom numbers, bond type, and parameters
# Constraints are bonds with kb = None, which can be extracted 
# using the category
class Bond(Bonded):
    def set(self,atoms,**kwargs):
        ids,r,ss,ca     = zip(*atoms)     
        self.atoms      = ids
        self.type       = 1
        self.positionCa = ca
        self.comments   = "%s(%s)-%s(%s)" % (r[0],ss[0],r[1],ss[1])
        # The category can be used to keep bonds sorted
        self.category   = kwargs.get("category")

        self.parameters = self.options['ForceField'].bbGetBond(r,ca,ss)
        # Backbone bonds also can be constraints. We could change the type further on, but this is more general.
        # Even better would be to add a new type: BB-Constraint
        if self.parameters[1] == None:
            self.category = 'Constraint'

    # Overriding __str__ method to suppress printing of bonds with Fc of 0
    def __str__(self):
        if len(self.parameters) > 1 and self.parameters[1] == 0:
            return ""
        return Bonded.__str__(self)


# Similar to the preceding class
class Angle(Bonded):
    def set(self,atoms,**kwargs):
        ids,r,ss,ca     = zip(*atoms)
        self.atoms      = ids
        self.type       = 2
        self.positionCa = ca
        self.comments   = "%s(%s)-%s(%s)-%s(%s)" % (r[0],ss[0],r[1],ss[1],r[2],ss[2])
        self.category   = kwargs.get("category")
        
        self.parameters = self.options['ForceField'].bbGetAngle(r,ca,ss)

# Similar to the preceding class
class Vsite(Bonded):
    def set(self,atoms,**kwargs):
        ids,r,ss,ca     = zip(*atoms)
        self.atoms      = ids
        self.type       = 1
        self.positionCa = ca
        self.comments   = "%s"% (r[0])
        self.category   = kwargs.get("category")
        self.parameters = kwargs.get("parameters") 

# Similar to the preceding class
class Exclusion(Bonded):
    def set(self,atoms,**kwargs):
        ids,r,ss,ca     = zip(*atoms)
        self.atoms      = ids
        self.positionCa = ca
        self.comments   = "%s"% (r[0])
        self.category   = kwargs.get("category")
        self.parameters = kwargs.get("parameters") 

# Similar to the preceding class
class Dihedral(Bonded):
    def set(self,atoms,**kwargs):
        ids,r,ss,ca     = zip(*atoms)
        self.atoms      = ids
        self.type       = 1
        self.positionCa = ca
        self.comments   = "%s(%s)-%s(%s)-%s(%s)-%s(%s)" % (r[0],ss[0],r[1],ss[1],r[2],ss[2],r[3],ss[3])
        self.category   = kwargs.get("category")

        if ''.join(i for i in ss) == 'FFFF':
            # Collagen
            self.parameters = self.options['ForceField'].bbDihedDictD['F']
        elif ''.join(i for i in ss) == 'EEEE' and self.options['ExtendedDihedrals']:
            # Use dihedrals
            self.parameters = self.options['ForceField'].bbDihedDictD['E']
        elif set(ss).issubset("H123"):
            # Helix
            self.parameters = self.options['ForceField'].bbDihedDictD['H']
        else:
            self.parameters = None


# This list allows to retrieve Bonded class items based on the category
# If standard, dictionary type indexing is used, only exact matches are
# returned. Alternatively, partial matching can be achieved by setting
# a second 'True' argument. 
class CategorizedList(list):
    def __getitem__(self,tag): 
        if type(tag) == int:
            # Call the parent class __getitem__
            return list.__getitem__(self,tag)

        if type(tag) == str:
            return [i for i in self if i.category == tag]

        if tag[1]:
            return [i for i in self if tag[0] in i.category]
        else:
            return [i for i in self if i.category == tag[0]]


class Topology:
    def __init__(self,other=None,options=None,name=""):
        self.name        = ''
        self.nrexcl      = 1
        self.atoms       = CategorizedList()
        self.vsites      = CategorizedList() 
        self.exclusions  = CategorizedList() 
        self.bonds       = CategorizedList()
        self.angles      = CategorizedList()
        self.dihedrals   = CategorizedList()
        self.impropers   = CategorizedList()
        self.constraints = CategorizedList()
        self.posres      = CategorizedList()
        self.sequence    = []
        self.secstruc    = ""
        # Okay, this is sort of funny; we will add a 
        #   #define mapping virtual_sitesn
        # to the topology file, followed by a header
        #   [ mapping ]
        self.mapping     = []
        # For multiscaling we have to keep track of the number of 
        # real atoms that correspond to the beads in the topology
        self.natoms      = 0        
        self.multiscale  = options['multi']

        if options:
            self.options = options
        else:
            self.options = {}

        if not other:
            # Returning an empty instance
            return
        elif isinstance(other,Topology):
            for attrib in ["atoms","vsites","bonds","angles","dihedrals","impropers","constraints","posres"]:
                setattr(self,attrib,getattr(other,attrib,[]))
        elif isinstance(other,IO.Chain):
            if other.type() == "Protein":
                self.fromAminoAcidSequence(other)
            elif other.type == "Nucleic":
                # Currently there are no Martini Nucleic Acids
                self.fromNucleicAcidSequence(other)
            elif other.type == "Mixed":
                logging.warning('Mixed Amino Acid /Nucleic Acid chains are not yet implemented')
                # How can you have a mixed chain?
                # Well, you could get a covalently bound lipid or piece of DNA to a protein :S
                # But how to deal with that?
                # Probably one should separate the chains into blocks of specified type,
                # determine the locations of links, then construct the topologies for the 
                # blocks and combine them according to the links.
                pass
            else:
                # This chain should not be polymeric, but a collection of molecules
                # For each unique residue type fetch the proper moleculetype 
                self.fromMoleculeList(other)
        if name:
            self.name = name

    def __iadd__(self,other):
        if not isinstance(other,Topology):
            other = Topology(other)
        shift     = len(self.atoms)
        last      = self.atoms[-1]
        atoms     = zip(*other.atoms)
        atoms[0]  = [i+shift for i in atoms[0]]   # Update atom numbers
        atoms[2]  = [i+last[2] for i in atoms[2]] # Update residue numbers
        atoms[5]  = [i+last[5] for i in atoms[5]] # Update charge group numbers
        self.atoms.extend(zip(*atoms))
        for attrib in ["bonds","vsites","angles","dihedrals","impropers","constraints"]:
            getattr(self,attrib).extend([source+shift for source in getattr(other,attrib)])
        return self

    def __add__(self,other):
        out = Topology(self)
        if not isinstance(other,Topology):
            other = Topology(other)
        out += other
        return out

    def __str__(self):
        if self.multiscale:
             out  = [ '; MARTINI (%s) Multiscale virtual sites topology section for "%s"' %(self.options['ForceField'].name,self.name) ]
        else:
             string  = '; MARTINI (%s) Coarse Grained topology file for "%s"' %(self.options['ForceField'].name, self.name)
             string += '\n; Created by martinize.py version %s \n; Using the following options:  ' %(self.options['Version'])
             string += ' '.join(self.options['Arguments'])
             out  = [ string ]
        if self.sequence:
            out += [
                '; Sequence:',
                '; ' + ''.join([ MAP.AA321.get(AA) for AA in self.sequence ]),
                '; Secondary Structure:',
                '; ' + self.secstruc,
                ]
        
        # Do not print a molecule name when multiscaling
        # In that case, the topology created here needs to be appended
        # at the end of an atomistic moleculetype
        if not self.multiscale:
            out += [ '\n[ moleculetype ]',
                     '; Name         Exclusions',  
                     '%-15s %3d' % (self.name,self.nrexcl)]

        out.append('\n[ atoms ]')

        # For virtual sites and dummy beads we have to be able to specify the mass.
        # Thus we need two different format strings:
        fs8 = '%5d %5s %5d %5s %5s %5d %7.4f ; %s'  
        fs9 = '%5d %5s %5d %5s %5s %5d %7.4f %7.4f ; %s'  
        out.extend([len(i)==9 and fs9%i or fs8%i for i in self.atoms])

        # Print out the vsites only if they excist. Right now it can only be type 1 virual sites.
        vsites = [str(i) for i in self.vsites]
        if vsites:
            out.append('\n[ virtual_sites2 ]')
            out.extend(vsites)

        # Print out the exclusions only if they excist.
        exclusions = [str(i) for i in self.exclusions]
        if exclusions:
            out.append('\n[ exclusions ]')
            out.extend(exclusions)

        if self.multiscale:
            out += ['\n;\n; Coarse grained to atomistic mapping\n;',
                    '#define mapping virtual_sitesn',
                    '[ mapping ]']
            for i,j in self.mapping:
                out.append( ("%5d     2 "%i)+" ".join(["%5d"%k for k in j]) )
            
            logging.info('Created virtual sites section for multiscaled topology')
            return "\n".join(out)

        # Bonds in order: backbone, backbone-sidechain, sidechain, short elastic, long elastic        
        out.append("\n[ bonds ]")       
        # Backbone-backbone
        bonds = [str(i) for i in self.bonds["BB"]]
        if bonds:
            out.append("; Backbone bonds")
            out.extend(bonds)
        # Rubber Bands
        bonds = [str(i) for i in self.bonds["Rubber",True]]
        if bonds:
            # Add a CPP style directive to allow control over the elastic network
            out.append("#ifdef RUBBER_BANDS")
            out.append("#ifndef RUBBER_FC\n#define RUBBER_FC %f\n#endif"%self.options['ElasticMaximumForce'])
            out.extend(bonds)
            out.append("#endif")
        # Backbone-Sidechain/Sidechain-Sidechain
        bonds = [str(i) for i in self.bonds["SC"]]
        if bonds:
            out.append("; Sidechain bonds")
            out.extend(bonds)
        # Short elastic/Long elastic
        bonds = [str(i) for i in self.bonds["Elastic short"]]
        if bonds:
            out.append("; Short elastic bonds for extended regions")
            out.extend(bonds)
        bonds = [str(i) for i in self.bonds["Elastic long"]]
        if bonds:
            out.append("; Long elastic bonds for extended regions")
            out.extend(bonds)
        # Cystine bridges
        bonds = [str(i) for i in self.bonds["Cystine"]]
        if bonds:
            out.append("; Cystine bridges")
            out.extend(bonds)
        # Other links
        bonds = [str(i) for i in self.bonds["Link"]]
        if bonds:
            out.append("; Links/Cystine bridges")
            out.extend(bonds)

        # Constraints
        out.append("\n[ constraints ]")
        out.extend([str(i) for i in self.bonds["Constraint"]])

        # Angles
        out.append("\n[ angles ]")
        out.append("; Backbone angles")
        out.extend([str(i) for i in self.angles["BBB"]])
        out.append("; Backbone-sidechain angles")
        out.extend([str(i) for i in self.angles["BBS"]])
        out.append("; Sidechain angles")
        out.extend([str(i) for i in self.angles["SC"]])

        # Dihedrals
        out.append("\n[ dihedrals ]")
        out.append("; Backbone dihedrals")
        out.extend([str(i) for i in self.dihedrals["BBBB"] if i.parameters])
        out.append("; Sidechain improper dihedrals")
        out.extend([str(i) for i in self.dihedrals["SC"] if i.parameters])

        # Postition Restraints
        if self.posres:
            out.append("\n#ifdef POSRES")
            out.append("#ifndef POSRES_FC\n#define POSRES_FC %.2f\n#endif"%self.options['PosResForce'])
            out.append(" [ position_restraints ]")
            out.extend(['  %5d    1    POSRES_FC    POSRES_FC    POSRES_FC'%i for i in self.posres])
            out.append("#endif")

        logging.info('Created coarsegrained topology')
        return "\n".join(out)

  
    # The sequence function can be used to generate the topology for 
    # a sequence :) either given as sequence or as chain
    def fromAminoAcidSequence(self,sequence,secstruc=None,links=None,breaks=None,
                              mapping=None,rubber=False,multi=False):
        # Shift for the atom numbers of the atomistic part in a chain 
        # that is being multiscaled
        shift = 0
        # First check if we get a sequence or a Chain instance
        if isinstance(sequence, IO.Chain):
            chain         = sequence
            links         = chain.links
            breaks        = chain.breaks
            # If the mapping is not specified, the actual mapping is taken,
            # used to construct the coarse grained system from the atomistic one.
            # The function argument "mapping" could be used to use a default 
            # mapping scheme in stead, like the mapping for the GROMOS96 force field.
            mapping = mapping           or chain.mapping
            multi   = self.options['multi']  or chain.multiscale
            self.secstruc = chain.sstypes or len(chain)*"C"
            self.sequence = chain.sequence
            # If anything hints towards multiscaling, do multiscaling
            self.multiscale = self.multiscale or chain.multiscale or multi
            if self.multiscale:
                shift        = self.natoms
                self.natoms += len(chain.atoms())
        elif not secstruc:
            # If no secondary structure is provided, set all to coil
            chain         = None
            self.secstruc = len(self.sequence)*"C"
        else:
            # If a secondary structure is provided, use that. chain is none.
            chain         = None
            self.secstruc = secstruc

        logging.debug(self.secstruc)
        logging.debug(self.sequence)

        # Fetch the sidechains
        # Pad with empty lists for atoms, bonds, angles 
        # and dihedrals, and take the first four lists out
        # This will avoid errors for residues for which 
        # these are not defined.

        sc = [(self.options['ForceField'].sidechains[res]+5*[[]])[:5] for res in self.sequence]

        # ID of the first atom/residue
        # The atom number and residue number follow from the last 
        # atom c.q. residue id in the list processed in the topology
        # thus far. In the case of multiscaling, the real atoms need 
        # also be accounted for.
        startAtom = self.natoms + 1 
        startResi = self.atoms and self.atoms[-1][2]+1 or 1

        # Backbone bead atom IDs
        bbid = [startAtom]
        for i in zip(*sc)[0]:
            bbid.append(bbid[-1]+len(i)+1)

        # Calpha positions, to get Elnedyn BBB-angles and BB-bond lengths
        # positionCa = [residue[1][4:] for residue in chain.residues]
        # The old method (line above) assumed no hydrogens: Ca would always be
        # the second atom of the residue. Now we look at the name.
        positionCa = []
        for residue in chain.residues:
            for atom in residue:
                if atom[0] == "CA":
                    positionCa.append(atom[4:])

        # Residue numbers for this moleculetype topology
        resid = range(startResi,startResi+len(self.sequence))     
        
        # This contains the information for deriving backbone bead types,
        # bb bond types, bbb/bbs angle types, and bbbb dihedral types and
        # Elnedyn BB-bondlength BBB-angles
        seqss = zip(bbid,self.sequence,self.secstruc,positionCa)

        # Fetch the proper backbone beads          
        bb = [self.options['ForceField'].bbGetBead(res,typ) for num,res,typ,Ca in seqss]

        # If termini need to be charged, change the bead types
        if not self.options['NeutralTermini']:
            bb[0]  ="Qd"
            bb[-1] = "Qa"

        # If breaks need to be charged, change the bead types 
        if self.options['ChargesAtBreaks']:
            for i in breaks:
                bb[i]   = "Qd"
                bb[i-1] = "Qa"

        # For backbone parameters, iterate over fragments, inferred from breaks
        for i,j in zip([0]+breaks,breaks+[-1]):
            # Extract the fragment
            frg = j==-1 and seqss[i:] or seqss[i:j]

            # Iterate over backbone bonds
            self.bonds.extend([Bond(pair,category="BB",options=self.options,) for pair in zip(frg,frg[1:])])

            # Iterate over backbone angles
            # Don't skip the first and last residue in the fragment
            self.angles.extend([Angle(triple,options=self.options,category="BBB") for triple in zip(frg,frg[1:],frg[2:])])

            # Get backbone quadruples
            quadruples = zip(frg,frg[1:],frg[2:],frg[3:])

            # No i-1,i,i+1,i+2 interactions defined for Elnedyn
            if self.options['ForceField'].UseBBBBDihedrals:
                # Process dihedrals
                for q in quadruples:
                    id,rn,ss,ca = zip(*q)
                    # Maybe do local elastic networks
                    if ss == ("E","E","E","E") and not self.options['ExtendedDihedrals']:
                        # This one may already be listed as the 2-4 bond of a previous one
                        if not (id[0],id[2]) in self.bonds:
                            self.bonds.append(Bond(options=self.options,atoms=(id[0],id[2]),parameters=self.options['ForceField'].ebonds['short'],type=1,
                                                   comments="%s(%s)-%s(%s) 1-3"%(rn[0],id[0],rn[2],id[2]),
                                                   category="Elastic short"))
                        self.bonds.append(Bond(options=self.options,atoms=(id[1],id[3]),parameters=self.options['ForceField'].ebonds['short'],type=1,
                                               comments="%s(%s)-%s(%s) 2-4"%(rn[1],id[1],rn[3],id[3]),
                                               category="Elastic short"))
                        self.bonds.append(Bond(options=self.options,atoms=(id[0],id[3]),parameters=self.options['ForceField'].ebonds['long'],type=1,
                                               comments="%s(%s)-%s(%s) 1-4"%(rn[0],id[0],rn[3],id[3]),
                                               category="Elastic long"))
                    else:
                        # Since dihedrals can return None, we first collect them separately and then
                        # add the non-None ones to the list
                        dihed = Dihedral(q,options=self.options,category="BBBB")
                        if dihed:
                            self.dihedrals.append(dihed)

            # Elnedyn does not use backbone-backbone-sidechain-angles
            if self.options['ForceField'].UseBBSAngles:
                # Backbone-Backbone-Sidechain angles
                # If the first residue has a sidechain, we take SBB, otherwise we skip it
                # For other sidechains, we 'just' take BBS
                if len(frg) > 1 and frg[1][0]-frg[0][0] > 1:
                    self.angles.append(Angle(options=self.options,atoms=(frg[0][0]+1,frg[0][0],frg[1][0]),parameters=self.options['ForceField'].bbsangle,type=2,
                                            comments="%s(%s)-%s(%s) SBB"%(frg[0][1],frg[0][2],frg[1][1],frg[1][2]),
                                            category="BBS"))
    
                # Start from first residue: connects sidechain of second residue
                for (ai,ni,si,ci),(aj,nj,sj,cj),s in zip(frg[0:],frg[1:],sc[1:]):
                    if s[0]:
                        self.angles.append(Angle(options=self.options,atoms=(ai,aj,aj+1),parameters=self.options['ForceField'].bbsangle,type=2,
                                                comments="%s(%s)-%s(%s) SBB"%(ni,si,nj,sj),
                                                category="BBS"))
           
        # Now do the atom list, and take the sidechains along
        #
        # AtomID AtomType ResidueID ResidueName AtomName ChargeGroup Charge ; Comments
        # 
        atid  = startAtom
        for resi,resname,bbb,sidechn,ss in zip(resid,self.sequence,bb,sc,self.secstruc):
            scatoms, bon_par, ang_par, dih_par, vsite_par = sidechn

            # Side chain bonded terms
            # Collect bond, angle and dihedral connectivity
            bon_con,ang_con,dih_con,vsite_con = (self.options['ForceField'].connectivity[resname]+4*[[]])[:4]

            # Side Chain Bonds/Constraints
            for atids,par in zip(bon_con,bon_par):
                if par[1] == None:
                    self.bonds.append(Bond(options=self.options,atoms=atids,parameters=[par[0]],type=1,
                                           comments=resname,category="Constraint"))
                else:
                    self.bonds.append(Bond(options=self.options,atoms=atids,parameters=par,type=1,
                                           comments=resname,category="SC"))
                # Shift the atom numbers
                self.bonds[-1] += atid

            # Side Chain Angles
            for atids,par in zip(ang_con,ang_par):
                self.angles.append(Angle(options=self.options,atoms=atids,parameters=par,type=2,
                                         comments=resname,category="SC"))
                # Shift the atom numbers
                self.angles[-1] += atid

            # Side Chain Dihedrals
            for atids,par in zip(dih_con,dih_par):
                self.dihedrals.append(Dihedral(options=self.options,atoms=atids,parameters=par,type=2,
                                               comments=resname,category="SC"))
                # Shift the atom numbers
                self.dihedrals[-1] += atid

            # Side Chain V-Sites
            for atids,par in zip(vsite_con,vsite_par):
                self.vsites.append(Vsite(options=self.options,atoms=atids,parameters=par,type=1,
                                               comments=resname,category="SC"))
                # Shift the atom numbers
                self.vsites[-1] += atid
            
            # Side Chain exclusions
            # The new polarizable forcefield give problems with the charges in the sidechain, if the backbone is also charged.
            # To avoid that, we add explicit exclusions
            if bbb in self.options['ForceField'].charges.keys() and resname in self.options['ForceField'].mass_charge.keys():
                for i in [i for i, d in enumerate(scatoms) if d=='D']:
                    self.exclusions.append(Exclusion(options=self.options,atoms=(atid,i+atid+1),comments='%s(%s)'%(resname,resi),parameters=(None,)))

            # All residue atoms
            counter = 0  # Counts over beads
            for atype,aname in zip([bbb]+list(scatoms),MAP.CoarseGrained.residue_bead_names):
                if self.multiscale:
                    atype,aname = "v"+atype,"v"+aname
                # If mass or charge diverse, we adopt it here. 
                # We don't want to do this for BB beads because of charged termini.
                if resname in self.options['ForceField'].mass_charge.keys() and counter != 0:
                    M,Q = self.options['ForceField'].mass_charge[resname]
                    aname = Q[counter-1]>0 and 'SCP' or Q[counter-1]<0 and 'SCN' or aname
                    self.atoms.append((atid,atype,resi,resname,aname,atid,Q[counter-1],M[counter-1],ss))
                else:
                    self.atoms.append((atid,atype,resi,resname,aname,atid,self.options['ForceField'].charges.get(atype,0),ss))
                # Doing this here save going over all the atoms onesmore.
                # Generate position restraints for all atoms or Backbone beads only.
                if 'all' in self.options['PosRes']:
                    self.posres.append((atid)) 
                elif aname in self.options['PosRes']:
                    self.posres.append((atid))
                if mapping:
                    self.mapping.append((atid,[i+shift for i in mapping[counter]]))
                atid    += 1
                counter += 1

        # The rubber bands are best applied outside of the chain class, as that gives
        # more control when chains need to be merged. The possibility to do it on the 
        # chain level is retained to allow building a complete chain topology in 
        # a straightforward manner after importing this script as module.
        if rubber and chain:
            rubberList = rubberBands(
                [(i[0],j[4:7]) for i,j in zip(self.atoms,chain.cg()) if i[4] in ElasticBeads],
                ElasticLowerBound,ElasticUpperBound,
                ElasticDecayFactor,ElasticDecayPower,
                ElasticMaximumForce,ElasticMinimumForce)
            self.bonds.extend([Bond(i,options=self.options,type=6,category="Rubber band") for i in rubberList])
        
        # Note the equivalent of atomistic atoms that have been processed 
        if chain and self.multiscale:
            self.natoms += len(chain.atoms())

    def fromNucleicAcidSequence(self,other):
        logging.warning('Nucleic Acid parameters are not available in MARTINI. Maybe *you* should create them?')
        pass

    def fromMoleculeList(self,other):
        pass

