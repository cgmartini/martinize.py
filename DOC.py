###################################
## 1 # OPTIONS AND DOCUMENTATION ##  -> @DOC <-
###################################

import martinize 
    
# This is a simple and versatily option class that allows easy
# definition and parsing of options. 
class Option:
    def __init__(self,func=str,num=1,default=None,description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description
    def __nonzero__(self): 
        if self.func == bool:
            return self.value != False
        return bool(self.value)
    def __str__(self):
        return self.value and str(self.value) or ""
    def setvalue(self,v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [ self.func(i) for i in v ]
    

# Lists for gathering arguments to options that can be specified 
# multiple times on the command line.
lists = {
    'cystines': [],
    'merges'  : [],
    'links'   : [],
    'multi'   : [],
    }

# List of 
options = [
#   NOTE: Options marked with (+) can be given multiple times on the command line
#   option              type number default description
    """
Primary input/output
--------------------
The input file (-f) should be a coordinate file in PDB or GROMOS
format. The format is inferred from the structure of the file, also
allowing reading MEAD/GRASP (.pqr) files. The input can also be
provided through stdin, allowing piping of structures, e.g. using grep
or sed to make a selection of atoms. This can be useful for removing
HETATM records: 

grep -v '^HETATM' input.pdb | martinize.py

The input structure can have multiple frames/models. If an output
structure file (-x) is given, each frame will be coarsegrained,
resulting in a multimodel output structure. Having multiple frames may
also affect the topology. If secondary structure is determined
internally, the structure will be averaged over the frames. Likewise,
interatomic distances, as used for backbone bond lengths in Elnedyn
and in elastic networks, are also averaged over the frames available.

If an output file (-o) is indicated for the topology, that file will
be used for the master topology, using #include statements to link the
moleculetype definitions, which are written to separate files. If no
output filename is given, the topology and the moleculetype
definitions are written to stdout.
""",
    ("-f",        Option(str,             1,     None, "Input file (PDB|GRO)")),
    ("-o",        Option(str,             1,     None, "Output topology (TOP)")),
    ("-x",        Option(str,             1,     None, "Output coarse grained structure (PDB)")),
    ("-n",        Option(str,             1,     None, "Output index file")),
    ("-v",        Option(bool,            0,    False, "Verbose. Be load and noisy.")), 
    ("-h",        Option(bool,            0,    False, "Display this help.")),
    """
Secondary structure
-------------------
The secondary structure plays a central role in the assignment of atom
types and bonded interactions in MARTINI. Martinize allows
specification of the secondary structure as a string (-ss), or as a
file containing a specification in GROMACS' ssdump format
(-ss). Alternatively, DSSP can be used for an on-the-fly assignment of
the secondary structure. For this, the option -dssp has to be used
giving the location of the executable as the argument. 
The option -collagen will set the whole structure to collagen. If this
is not what you want (eg only part of the structure is collagen, you
can give a secondary structure file/string (-ss) and specifiy collagen
as "F". Parameters for collagen are taken from: Gautieri et al., 
J. Chem. Theory Comput, 2010, 6, 1210-1218. 
With multimodel input files, the secondary structure as determined with
DSSP or PyMOL will be averaged over the frames. In this case, a cutoff
can be specified (-ssc) indicating the fraction of frames to match a
certain secondary structure type for designation.
""",
    ("-ss",       Option(str,             1,     None, "Secondary structure (File or string)")),
    ("-ssc",      Option(float,           1,      0.5, "Cutoff fraction for ss in case of ambiguity (default: 0.5).")),
    ("-dssp",     Option(str,             1,     None, "DSSP executable for determining structure")),
#    ("-pymol",    Option(str,             1,     None, "PyMOL executable for determining structure")),
    ("-collagen", Option(bool,            0,    False, "Use collagen parameters")),
    """
Topology
--------

Several options are available to tune the resulting topology. By
default, termini are charged, and chain breaks are kept neutral. This
behaviour can be changed using -nt and -cb, respectively.

Disulphide bridges can be specified using -cys. This option can be
given multiple times on the command line. The argument is a pair of
cysteine residues, using the format
chain/resn/resi,chain/resn/resi. For disulphide bridges, the residue
name is not required, and the chain identifier is optional. If no
chain identifier is given, all matching residue pairs will be checked,
and pairs within the cutoff distance (0.22 nm) will be linked. It is
also possible to let martinize detect cysteine pairs based on this
cut-off distance, by giving the keyword 'auto' as argument to -cys.
Alternatively, a different cut-off distance can be specified, which
will also trigger a search of pairs satisfying the distance
criterion.

In addition to cystine bridges, links between other atoms can be
specified using -link. This requires specification of the atoms, using
the format
chain/resi/resn/atom,chain/resi/resn/atom,bondlength,forceconstant.
If only two atoms are given, a constraint will be added with length
equal to the (average) distance in the coordinate file. If a bond
length is added, but no force constant, then the bondlength will be
used to set a constraint.

Linking atoms requires that the atoms are part of the same
moleculetype. Therefore any link between chains will cause the chains
to be merged. Merges can also be specified explicitly, using the
option -merge with a comma-separated list of chain identifiers to be
joined into one moleculetype. The option -merge can be used several
times. Note that specifying a chain in several merge groups will cause
all chains involved to be merged into a single moleculetype.

The moleculetype definitions are written to topology include (.itp)
files, using a name consisting of the molecule class (e.g. Protein)
and the chain identifier. With -name a name can be specified instead.
By default, martinize only writes a moleculetype for each unique
molecule, inferred from the sequence and the secondary structure
definition. It is possible to force writing a moleculetype definition
for every single molecule, using -sep.

The option -p can be used to write position restraints, using the 
force constant specified with -pf, which is set to 1000 kJ/mol 
by default.

For stability, elastic bonds are used to retain the structure of 
extended strands. The option -ed causes dihedrals to be used 
instead.

Different forcefields can be specified with -ff. All the parameters and
options belonging to that forcefield  will be set (eg. bonded interactions,
BB-bead positions, Elastic Network, etc.). By default martini 2.1 is
used.
""",
    ("-nt",       Option(bool,                     0,      False, "Set neutral termini (charged is default)")), 
    ("-cb",       Option(bool,                     0,      False, "Set charges at chain breaks (neutral is default)")), 
    ("-cys",      Option(lists['cystines'].append, 1,       None, "Disulphide bond (+)")),
    ("-link",     Option(lists['links'].append,    1,       None, "Link (+)")),
    ("-merge",    Option(lists['merges'].append,   1,       None, "Merge chains: e.g. -merge A,B,C (+)")),
#    ("-mixed",    Option(bool,                     0,      False, "Allow chains of mixed type (default: False)")),
    ("-name",     Option(str,                      1,       None, "Moleculetype name")),
    ("-p",        Option(str,                      1,     'None', "Output position restraints (None/All/Backbone) (default: None)")),
    ("-pf",       Option(float,                    1,       1000, "Position restraints force constant (default: 1000 kJ/mol/nm^2)")),
    ("-ed",       Option(bool,                     0,      False, "Use dihedrals for extended regions rather than elastic bonds)")),
    ("-sep",      Option(bool,                     0,      False, "Write separate topologies for identical chains.")),
    ("-ff",       Option(str,                      1,'martini21', "Which forcefield to use: "+' ,'.join(n for n in martinize.forcefields))),
    """
Elastic network
---------------

Martinize can write an elastic network for atom pairs within a cutoff
distance. The force constant (-ef) and the upper distance bound (-eu) 
can be speficied. If a force field with an intrinsic Elastic
network is specified (eg. Elnedyn) with -ff, -elastic in implied and
the default values for the force constant and upper cutoff are used.
However, these can be overwritten.
""",
# Fij = Fc exp( -a (rij - lo)**p )
    ("-elastic",  Option(bool,            0,    False, "Write elastic bonds")),
#    ("-elnedyn",  Option(bool,            0,    False, "Use Elnedyn mapping and parameters with elastic network")),
    ("-ef",       Option(float,           1,      500, "Elastic bond force constant Fc")),
    ("-el",       Option(float,           1,        0, "Elastic bond lower cutoff: F = Fc if rij < lo")),
    ("-eu",       Option(float,           1,     0.90, "Elastic bond upper cutoff: F = 0  if rij > up")),
    ("-ea",       Option(float,           1,        0, "Elastic bond decay factor a")),
    ("-ep",       Option(float,           1,        1, "Elastic bond decay power p")),
    ("-em",       Option(float,           1,        0, "Remove elastic bonds with force constant lower than this")),
    ("-eb",       Option(str,             1,     'BB', "Comma separated list of bead names for elastic bonds")),
#    ("-cgo",      Option(str,             1,     None, "PyMOL CGO file for elastic network visualization")), 
#    ("-hetatm",   Option(bool,            0,    False, "Include HETATM records from PDB file (Use with care!)")),
    """
Multiscaling
------------

Martinize can process a structure to yield a multiscale system,
consisting of a coordinate file with atomistic parts and
corresponding, overlaid coarsegrained parts. For chains that are
multiscaled, rather than writing a full moleculetype definition, 
additional [atoms] and [virtual_sitesn] sections are written, to 
be appended to the atomistic moleculetype definitions. 
The option -multi can be specified multiple times, and takes a chain
identifier as argument. Alternatively, the keyword 'all' can be given
as argument, causing all chains to be multiscaled.
""",
    ("-multi",    Option(lists['multi'].append,    1,     None, "Chain to be set up for multiscaling (+)")),
"========================================================================\n"
    ]

## Martini Quotes
martiniq = [
    ("Robert Benchley",
     "Why don't you get out of that wet coat and into a dry martini?"),
    ("James Thurber",
     "One martini is all right, two is two many, three is not enough"),
    ("Philip Larkin",
     "The chromatic scale is what you use to give the effect of drinking a quinine martini and having an enema simultaneously."),
    ("William Emerson, Jr.",
     "And when that first martini hits the liver like a silver bullet, there is a sigh of contentment that can be heard in Dubuque."),
    ("Alec Waugh",
     "I am prepared to believe that a dry martini slightly impairs the palate, but think what it does for the soul."),
    ("Gerald R. Ford",
     "The three-martini lunch is the epitome of American efficiency. Where else can you get an earful, a bellyful and a snootful at the same time?"),
    ("P. G. Wodehouse",
     "He was white and shaken, like a dry martini."),
    ]

desc = ""
    
def help():
    import sys
    from martinize import forcefields
    for item in options:
        if type(item) == str:
            print item
    for item in options:
        if type(item) != str:
            print "%10s  %s"%(item[0],item[1].description)
    print
    sys.exit()
    
