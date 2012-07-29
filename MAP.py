##########################
## 4 # FG -> CG MAPPING ##  -> @MAP <-
##########################
import FUNC

# Amino acid codes:                                                                                  
AA3     = FUNC.spl("TRP TYR PHE HIS ARG LYS CYS ASP GLU ILE LEU MET ASN PRO HYP GLN SER THR VAL ALA GLY") #@#
AA1     = FUNC.spl("  W   Y   F   H   R   K   C   D   E   I   L   M   N   P   O   Q   S   T   V   A   G") #@#


# Dictionaries for conversion from one letter code to three letter code v.v.                         
AA123, AA321 = FUNC.hash(AA1,AA3),FUNC.hash(AA3,AA1)                                                           


# Residue classes:
protein = AA3
water   = FUNC.spl("HOH SOL TIP")
lipids  = FUNC.spl("DPP DHP DLP DMP DSP POP DOP DAP DUP DPP DHP DLP DMP DSP PPC DSM DSD DSS")
nucleic = FUNC.spl("DAD DCY DGU DTH ADE CYT GUA THY URA DA DC DG DT")


residueTypes = dict(
    [(i,"Protein") for i in protein ]+
    [(i,"Water")   for i in water   ]+
    [(i,"Lipid")   for i in lipids  ]+
    [(i,"Nucleic") for i in nucleic ]
    )

class CoarseGrained:
    # Class for mapping an atomistic residue list to a coarsegrained one
    # Should get an __init__ function taking a residuelist, atomlist, Pymol selection or ChemPy model
    # The result should be stored in a list-type attribute
    # The class should have pdbstr and grostr methods

    # Standard mapping groups
    # Protein backbone
    bb        = "N CA C O H H1 H2 H3 O1 O2"                                                                    #@#  
    # Lipid tails
    palmitoyl1    = FUNC.nsplit("C1B C1C C1D C1E","C1F C1G C1H C1I","C1J C1K C1L C1M","C1N C1O C1P")                #@#
    palmitoyl2    = FUNC.nsplit("C2B C2C C2D C2E","C2F C2G C2H C2I","C2J C2K C2L C2M","C2N C2O C2P")                #@#
    oleyl1        = FUNC.nsplit("C1B C1C C1D C1E","C1F C1G C1H","C1I C1J","C1K C1L C1M C1N","C1O C1P C1Q C1R")      #@#
    oleyl2        = FUNC.nsplit("C2B C2C C2D C2E","C2F C2G C2H","C2I C2J","C2K C2L C2M C2N","C2O C2P C2Q C2R")      #@#
    #lauroyl1      = []
    #stearoyl1     = []
    #arachidonoyl1 = []
    #linoleyl1     = []
    #hexanoyl1     = []
    # Lipid head groups      
    #phoshpatidylcholine      = 
    phosphatydilethanolamine = FUNC.nsplit("N H1 H2 H3 CA","CB P OA OB OC OD","CC CD OG C2A OH","CE OE C1A OF")     #@#
    phosphatidylglycerol     = FUNC.nsplit("H1 O1 CA H2 O2 CB","CC P OA OB OC OD","CD CE OG C2A OH","CF OE C1A OF") #@#
    #phosphatidylserine       =

    dna_bb = "P OP1 OP2 O5' O3'","C5' O4' C4'","C3' O3' C2' C1'"

    # This is the mapping dictionary
    # For each residue it returns a list, each element of which
    # lists the atom names to be mapped to the corresponding bead.
    # The order should be the standard order of the coarse grained
    # beads for the residue. Only atom names matching with those 
    # present in the list of atoms for the residue will be used
    # to determine the bead position. This adds flexibility to the
    # approach, as a single definition can be used for different 
    # states of a residue (e.g., GLU/GLUH).
    # For convenience, the list can be specified as a set of strings,
    # converted into a list of lists by 'FUNC.nsplit' defined above.
    mapping = {
        "ALA":  FUNC.nsplit(bb + " CB"),
        "CYS":  FUNC.nsplit(bb,"CB SG"),
        "ASP":  FUNC.nsplit(bb,"CB CG OD1 OD2"),
        "GLU":  FUNC.nsplit(bb,"CB CG CD OE1 OE2"),
        "PHE":  FUNC.nsplit(bb,"CB CG CD1 HD1","CD2 HD2 CE2 HE2","CE1 HE1 CZ HZ"),
        "GLY":  FUNC.nsplit(bb),
        "HIS":  FUNC.nsplit(bb,"CB CG","CD2 HD2 NE2 HE2","ND1 HD1 CE1 HE1"),
        "ILE":  FUNC.nsplit(bb,"CB CG1 CG2 CD CD1"),
        "LYS":  FUNC.nsplit(bb,"CB CG CD","CE NZ HZ1 HZ2 HZ3"),
        "LEU":  FUNC.nsplit(bb,"CB CG CD1 CD2"),
        "MET":  FUNC.nsplit(bb,"CB CG SD CE"),
        "ASN":  FUNC.nsplit(bb,"CB CG ND1 ND2 OD1 OD2 HD11 HD12 HD21 HD22"),
        "PRO":  FUNC.nsplit(bb,"CB CG CD"),
        "HYP":  FUNC.nsplit(bb,"CB CG CD OD"),
        "GLN":  FUNC.nsplit(bb,"CB CG CD OE1 OE2 NE1 NE2 HE11 HE12 HE21 HE22"),
        "ARG":  FUNC.nsplit(bb,"CB CG CD","NE HE CZ NH1 NH2 HH11 HH12 HH21 HH22"),    
        "SER":  FUNC.nsplit(bb,"CB OG HG"),
        "THR":  FUNC.nsplit(bb,"CB OG1 HG1 CG2"),
        "VAL":  FUNC.nsplit(bb,"CB CG1 CG2"),
        "TRP":  FUNC.nsplit(bb,"CB CG CD2","CD1 HD1 NE1 HE1 CE2","CE3 HE3 CZ3 HZ3","CZ2 HZ2 CH2 HH2"),
        "TYR":  FUNC.nsplit(bb,"CB CG CD1 HD1","CD2 HD2 CE2 HE2","CE1 HE1 CZ OH HH"),
        "POPE": phosphatydilethanolamine + palmitoyl1 + oleyl2,
        "DOPE": phosphatydilethanolamine + oleyl1     + oleyl2,
        "DPPE": phosphatydilethanolamine + palmitoyl1 + palmitoyl2,
        "POPG": phosphatidylglycerol     + palmitoyl1 + oleyl2,
        "DOPG": phosphatidylglycerol     + oleyl1     + oleyl2,
        "DPPG": phosphatidylglycerol     + palmitoyl1 + palmitoyl2,
        "DA": FUNC.nsplit("P OP1 OP2 O5' O3'","C5' O4' C4'","C3' O3' C2' C1'","N9 C4","C8 N7 C5","C6 N6 N1","C2 N3"),
        "DG": FUNC.nsplit("P OP1 OP2 O5' O3'","C5' O4' C4'","C3' O3' C2' C1'","N9 C4","C8 N7 C5","C6 O6 N1","C2 N3 N2"),
        "DC": FUNC.nsplit("P OP1 OP2 O5' O3'","C5' O4' C4'","C3' O3' C2' C1'","N1 C6","C2 O2 N3","C4 O4 C5"),
        "DT": FUNC.nsplit("P OP1 OP2 O5' O3'","C5' O4' C4'","C3' O3' C2' C1'","N1 C6","C2 O2 N3","C4 O4 C5 C7"),
        }

    # Generic names for side chain beads
    residue_bead_names = FUNC.spl("BB SC1 SC2 SC3 SC4")

    # This dictionary contains the bead names for all residues,
    # following the order in 'mapping'
    names  = {
        "POPE": "NH3 PO4 GL1 GL2 C1A C2A C3A C4A C1B C2B D3B C4B C5B".split(),
        "POPG": "GLC PO4 GL1 GL2 C1A C2A C3A C4A C1B C2B D3B C4B C5B".split() 
        }
    # Add default bead names for all amino acids
    names.update([(i,("BB","SC1","SC2","SC3","SC4")) for i in AA3])

    # Add the default bead names for all DNA nucleic acids
    names.update([(i,("BB1","BB2","BB3","SC1","SC2","SC3","SC4")) for i in nucleic])

    # This dictionary allows determining four letter residue names
    # for ones specified with three letters, e.g., resulting from
    # truncation to adhere to the PDB format.
    # Each entry returns a prototypical test, given as a string,
    # and the residue name to be applied if eval(test) is True.
    # This is particularly handy to determine lipid types.
    # The test assumes there is a local or global array 'atoms'
    # containing the atom names of the residue in correct order.
    restest = {
        "POP": [('atoms[0] == "CA"', "POPG"),
                ('atoms[0] == "N"',  "POPE")]
        }

    # Crude mass for weighted average. No consideration of united atoms.
    # This will probably give only minor deviations, while also giving less headache
    mass = {'H': 1,'C': 12,'N': 14,'O': 16,'S': 32,'P': 31,'M': 0}

# Determine average position for a set of weights and coordinates
# This is a rather specific function that requires a list of items
# [(m,(x,y,z),id),..] and returns the weighted average of the 
# coordinates and the list of ids mapped to this bead
def aver(b):
    mwx,ids = zip(*[((m*x,m*y,m*z),i) for m,(x,y,z),i in b])              # Weighted coordinates     
    tm  = sum(zip(*b)[0])                                                 # Sum of weights           
    return [sum(i)/tm for i in zip(*mwx)],ids                             # Centre of mass           

# Return the CG beads for an atomistic residue, using the mapping specified above
# The residue 'r' is simply a list of atoms, and each atom is a list:
# [ name, resname, resid, chain, x, y, z ]
def map(r,ca2bb = False):
    p = CoarseGrained.mapping[r[0][1]]                                             # Mapping for this residue 
    if ca2bb: p[0] = ["CA"]                                                        # Elnedyn maps BB to CA, ca2bb is False or True
    # Get the name, mass and coordinates for all atoms in the residue
    a = [(i[0],CoarseGrained.mass.get(i[0][0],0),i[4:]) for i in r]                    
    # Store weight, coordinate and index for atoms that match a bead
    q = [[(m,coord,a.index((atom,m,coord))) for atom,m,coord in a if atom in i] for i in p]

    # Bead positions      
    return zip(*[aver(i) for i in q])

