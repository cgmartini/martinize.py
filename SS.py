#############################
## 5 # SECONDARY STRUCTURE ##  -> @SS <-
#############################
import logging, os, sys
import subprocess as subp
import FUNC, IO

#----+--------------------------------------+
## A | SECONDARY STRUCTURE TYPE DEFINITIONS |
#----+--------------------------------------+

# This table lists all coarse grained secondary structure types
# The following are matched lists. Make sure they stay matched.
# The lists do not need to be of the same length. The longer list
# will be truncated when combined with a shorter list, e.g. with
# dihedral definitions, which are not present for coil and termini
#
ss_names = {
 "F": "Collagenous Fiber",                                                                  #@#
 "E": "Extended structure (beta sheet)",                                                    #@#
 "H": "Helix structure",                                                                    #@#
 "1": "Helix start (H-bond donor)",                                                         #@#
 "2": "Helix end (H-bond acceptor)",                                                        #@#
 "3": "Ambivalent helix type (short helices)",                                              #@#
 "T": "Turn",                                                                               #@#
 "S": "Bend",                                                                               #@#
 "C": "Coil",                                                                               #@#
}

bbss = list(ss_names.keys())
bbss = FUNC.spl("  F     E     H     1     2     3     T     S     C")  # SS one letter


# The following dictionary contains secondary structure types as assigned by
# different programs. The corresponding Martini secondary structure types are
# listed in cgss
#
# NOTE:
#  Each list of letters in the dictionary ss should exactly match the list
#  in cgss.
#
ssdefs = {
    "dssp":  list(".HGIBETSC~"),             # DSSP one letter secondary structure code     #@#
    "pymol": list(".H...S...L"),             # Pymol one letter secondary structure code    #@#
    "gmx":   list(".H...ETS.C"),             # Gromacs secondary structure dump code        #@#
    "self":  list("FHHHEETSCC")              # Internal CG secondary structure codes        #@#
}
cgss     =   list("FHHHEETSCC")              # Corresponding CG secondary structure types   #@#


#----+-------------------------------------------+
## B | SECONDARY STRUCTURE PATTERN SUBSTITUTIONS |
#----+-------------------------------------------+


# For all structure types specific dihedrals may be used if four or
# more consecutive residues are assigned that type.

# Helix start and end regions are special and require assignment of
# specific types. The following pattern substitutions are applied
# (in the given order). A dot matches any other type.
# Patterns can be added to the dictionaries. This only makes sense
# if for each key in patterns there is a matching key in pattypes.
patterns = {
    "H": FUNC.pat(".H. .HH. .HHH. .HHHH. .HHHHH. .HHHHHH. .HHHHHHH. .HHHH HHHH.")                #@#
}
pattypes = {
    "H": FUNC.pat(".3. .33. .333. .3333. .13332. .113322. .1113222. .1111 2222.")                #@#
}


#----+----------+
## C | INTERNAL |
#----+----------+


# Pymol Colors
#          F   E   H   1   2   3   T   S   C
ssnum  = (13,  4,  2,  2,  2,  2,  6, 22,  0)                                             #@#

# Dictionary returning a number for a given type of secondary structure
# This can be used for setting the b-factor field for coloring
ss2num = FUNC.hash(bbss, ssnum)


# List of programs for which secondary structure definitions can be processed
programs = list(ssdefs.keys())


# Dictionaries mapping ss types to the CG ss types
ssd = dict([(i, FUNC.hash(ssdefs[i], cgss)) for i in programs])


# From the secondary structure dictionaries we create translation tables
# with which all secondary structure types can be processed. Anything
# not listed above will be mapped to C (coil).
# Note, a translation table is a list of 256 characters to map standard
# ascii characters to.
def tt(program):
    return "".join([ssd[program].get(chr(i), "C") for i in range(256)])


# The translation table depends on the program used to obtain the
# secondary structure definitions
sstt = dict([(i, tt(i)) for i in programs])


# The following translation tables are used to identify stretches of
# a certain type of secondary structure. These translation tables have
# every character, except for the indicated secondary structure, set to
# \x00. This allows summing the sequences after processing to obtain
# a single sequence summarizing all the features.
null = "\x00"
sstd = dict([(i, ord(i)*null+i+(255-ord(i))*null) for i in cgss])


# Pattern substitutions
def typesub(seq, patterns, types):
    seq = null+seq+null
    for i, j in zip(patterns, types):
        seq = seq.replace(i, j)
    return seq[1:-1]


# The following function translates a string encoding the secondary structure
# to a string of corresponding Martini types, taking the origin of the
# secondary structure into account, and replacing termini if requested.
def ssClassification(ss, program="dssp"):
    # Translate dssp/pymol/gmx ss to Martini ss
    ss  = ss.translate(sstt[program])
    # Separate the different secondary structure types
    sep = dict([(i, ss.translate(sstd[i])) for i in list(sstd.keys())])
    # Do type substitutions based on patterns
    # If the ss type is not in the patterns lists, do not substitute
    # (use empty lists for substitutions)
    typ = [typesub(sep[i], patterns.get(i, []), pattypes.get(i, [])) for i in list(sstd.keys())]
    # Translate all types to numerical values
    typ = [[ord(j) for j in list(i)] for i in typ]
    # Sum characters back to get a full typed sequence
    typ = "".join([chr(sum(i)) for i in zip(*typ)])
    # Return both the actual as well as the fully typed sequence
    return ss, typ


# The following functions are for determination of secondary structure,
# given a list of atoms. The atom format is generic and can be written out
# as PDB or GRO. The coordinates are in Angstrom.
# NOTE: There is the *OLD* DSSP and the *NEW* DSSP, which require
# different calls. The old version uses '--' to indicate reading from stdin
# whereas the new version uses '-i /dev/stdin'
def call_dssp(chain, atomlist, executable='dsspcmbi'):
    '''Get the secondary structure, by calling to dssp'''
    ssdfile = 'chain_%s.ssd' % chain.id

    try:
        if os.system(executable+" -V 2>/dev/null"):
            logging.debug("New version of DSSP; Executing '%s -i /dev/stdin -o %s'" % (executable, ssdfile))
            p = subp.Popen([executable, "-i", "/dev/stdin", "-o", ssdfile], stderr=subp.PIPE, stdout=subp.PIPE, stdin=subp.PIPE)
        else:
            logging.debug("Old version of DSSP; Executing '%s -- %s'" % (executable, ssdfile))
            p = subp.Popen([executable, "--", ssdfile], stderr=subp.PIPE, stdout=subp.PIPE, stdin=subp.PIPE)
    except OSError:
        logging.error("A problem occured calling %s." % executable)
        sys.exit(1)

    for atom in atomlist:
        if atom[0][:2] == 'O1': atom = ('O',)+atom[1:]
        if atom[0][0] != 'H' and atom[0][:2] != 'O2': p.stdin.write(IO.pdbOut(atom))
    p.stdin.write('TER\n')
    data = p.communicate()
    p.wait()
    main, ss = 0, ''
    for line in open(ssdfile).readlines():
        if main and not line[13] == "!": ss += line[16]
        if line[:15] == '  #  RESIDUE AA': main = 1
    return ss

ssDetermination = {
    "dssp": call_dssp
    }
