#######################
## 8 # STRUCTURE I/O ##  -> @IO <-
#######################
import logging, math, random, sys
import MAP, SS, FUNC

#----+---------+
## A | PDB I/O |
#----+---------+

d2r = 3.14159265358979323846264338327950288/180

# Reformatting of lines in structure file
pdbBoxLine  = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n"


def pdbBoxString(box):
    # Box vectors
    u, v, w  = box[0:3], box[3:6], box[6:9]

    # Box vector lengths
    nu, nv, nw = [math.sqrt(FUNC.norm2(i)) for i in (u, v, w)]

    # Box vector angles
    alpha = nv*nw == 0 and 90 or math.acos(FUNC.cos_angle(v, w))/d2r
    beta  = nu*nw == 0 and 90 or math.acos(FUNC.cos_angle(u, w))/d2r
    gamma = nu*nv == 0 and 90 or math.acos(FUNC.cos_angle(u, v))/d2r

    return pdbBoxLine % (10*FUNC.norm(u), 10*FUNC.norm(v), 10*FUNC.norm(w), alpha, beta, gamma)


def pdbAtom(a):
    ##01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ##ATOM   2155 HH11 ARG C 203     116.140  48.800   6.280  1.00  0.00
    if a.startswith("TER"):
        return 0
    # NOTE: The 27th field of an ATOM line in the PDB definition can contain an
    #       insertion code. We shift that 20 bits and add it to the residue number
    #       to ensure that residue numbers will be unique.
    ## ===> atom name,       res name,        res id,                        chain,
    atom = [a[12:16].strip(), a[17:20].strip(), int(a[22:26])+(ord(a[26])<<20), a[21],
    #             x,              y,              z
            float(a[30:38]), float(a[38:46]), float(a[46:54])]
    # If the chain identifier is empty, the chain is to None
    if atom[3].strip() == '':
        atom[3] = None
    return tuple(atom) 


def pdbOut(atom, i=1, **kwargs):
    # insc contains the insertion code, shifted by 20-bitwise.
    # This means there are multiple residues with the same "resi",
    # which we circumvent by subtracting "insc" from "resi".
    # At other places this subtraction as to be inverted.
    insc = atom[2] >> 20
    resi = atom[2]-(insc << 20)
    if atom[3] == None:
        chain = ' '
    else:
        chain = atom[3]
    pdbline = "ATOM  %5i  %-3s %3s%2s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f           %1s \n"
    if "ssid" in kwargs and type(kwargs["ssid"]) == type(int()):
        occupancy = kwargs["ssid"]
    else:
        occupancy = 40
    return pdbline % ((i, atom[0][:3], atom[1], chain, resi, chr(insc)) + atom[4:] + (1, occupancy, atom[0][0]))


def isPdbAtom(a):
    return a.startswith("ATOM") or (options["-hetatm"] and a.startswith("HETATM")) or a.startswith("TER")


def pdbBoxRead(a):
    fa, fb, fc, aa, ab, ac = [float(i) for i in a.split()[1:7]]
    ca, cb, cg, sg         = math.cos(d2r*aa), math.cos(d2r*ab), math.cos(d2r*ac), math.sin(d2r*ac)
    wx, wy                 = 0.1*fc*cb, 0.1*fc*(ca-cb*cg)/sg
    wz                     = math.sqrt(0.01*fc*fc - wx*wx - wy*wy)
    return [0.1*fa, 0, 0, 0.1*fb*cg, 0.1*fb*sg, 0, wx, wy, wz]


# Function for splitting a PDB file in chains, based
# on chain identifiers and TER statements
def pdbChains(pdbAtomList):
    chain = []
    for atom in pdbAtomList:
        if not atom:  # Was a "TER" statement
            if chain:
                yield chain
            else:
                logging.info("Skipping empty chain definition")
            chain = []
            continue
        if not chain or chain[-1][3] == atom[3]:
            chain.append(atom)
        else:
            yield chain
            chain = [atom]
    if chain:
        yield chain


# Simple PDB iterator
def pdbFrameIterator(streamIterator):
    title, atoms, box = [], [], []
    for i in streamIterator:
        if i.startswith("ENDMDL"):
            yield "".join(title), atoms, box
            title, atoms, box = [], [], []
        elif i.startswith("TITLE"):
            title.append(i)
        elif i.startswith("CRYST1"):
            box = pdbBoxRead(i)
        elif i.startswith("ATOM") or i.startswith("HETATM"):
            atoms.append(pdbAtom(i))
    if atoms:
        yield "".join(title), atoms, box


#----+---------+
## B | GRO I/O |
#----+---------+

groline = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"


def groBoxRead(a):
    b = [float(i) for i in a.split()] + 6*[0]                    # Padding for rectangular boxes
    return b[0], b[3], b[4], b[5], b[1], b[6], b[7], b[8], b[2]  # Return full definition xx,xy,xz,yx,yy,yz,zx,zy,zz


def groAtom(a):
    # In PDB files, there might by an insertion code. To handle this, we internally add
    # constant to all resids. To be consistent, we have to do the same for gro files.
    # 32 equal ord(' '), eg an empty insertion code
    constant = 32 << 20
    #012345678901234567890123456789012345678901234567890
    #    1PRN      N    1   4.168  11.132   5.291
    #  ===> atom name,        res name,          res id,             chain,
    return (a[10:15].strip(), a[5:10].strip(),   int(a[:5])+constant, None,
    #                x,                 y,                 z
            10*float(a[20:28]), 10*float(a[28:36]), 10*float(a[36:44]))


# Simple GRO iterator
def groFrameIterator(streamIterator):
    while True:
        try:
            title = streamIterator.next()
        except StopIteration:
            break
        natoms = streamIterator.next().strip()
        if not natoms:
            break
        natoms = int(natoms)
        atoms  = [groAtom(streamIterator.next()) for i in range(natoms)]
        box    = groBoxRead(streamIterator.next())
        yield title, atoms, box


#----+-------------+
## C | GENERAL I/O |
#----+-------------+

# It is not entirely clear where this fits in best.
# Called from main.
def getChargeType(resname, resid, choices):
    '''Get user input for the charge of residues, based on list with choises.'''
    print 'Which %s type do you want for residue %s:' % (resname, resid+1)
    for i, choice in choices.iteritems():
        print '%s. %s' % (i, choice)
    choice = None
    while choice not in choices.keys():
        choice = input('Type a number:')
    return choices[choice]


# *NOTE*: This should probably be a CheckableStream class that
# reads in lines until either of a set of specified conditions
# is met, then setting the type and from thereon functioning as
# a normal stream.
def streamTag(stream):
    # Tag the stream with the type of structure file
    # If necessary, open the stream, taking care of
    # opening using gzip for gzipped files

    # First check whether we have have an open stream or a file
    # If it's a file, check whether it's zipped and open it
    if type(stream) == str:
        if stream.endswith("gz"):
            logging.info('Read input structure from zipped file.')
            s = gzip.open(stream)
        else:
            logging.info('Read input structure from file.')
            s = open(stream)
    else:
        logging.info('Read input structure from command-line')
        s = stream

    # Read a few lines, but save them
    x = [s.readline(), s.readline()]
    if x[-1].strip().isdigit():
        # Must be a GRO file
        logging.info("Input structure is a GRO file. Chains will be labeled consecutively.")
        yield "GRO"
    else:
        # Must be a PDB file then
        # Could wind further to see if we encounter an "ATOM" record
        logging.info("Input structure is a PDB file.")
        yield "PDB"

    # Hand over the lines that were stored
    for i in x:
        yield i

    # Now give the rest of the lines from the stream
    for i in s:
        yield i


#----+-----------------+
## D | STRUCTURE STUFF |
#----+-----------------+

# This list allows to retrieve atoms based on the name or the index
# If standard, dictionary type indexing is used, only exact matches are
# returned. Alternatively, partial matching can be achieved by setting
# a second 'True' argument.
class Residue(list):
    def __getitem__(self, tag):
        if type(tag) == int:
            # Call the parent class __getitem__
            return list.__getitem__(self, tag)
        if type(tag) == str:
            for i in self:
                if i[0] == tag:
                    return i
            else:
                return
        if tag[1]:
            return [i for i in self if tag[0] in i[0]]  # Return partial matches
        else:
            return [i for i in self if i[0] == tag[0]]  # Return exact matches only


def residues(atomList):
    residue = [atomList[0]]
    for atom in atomList[1:]:
        if (atom[1] == residue[-1][1] and  # Residue name check
            atom[2] == residue[-1][2] and  # Residue id check
            atom[3] == residue[-1][3]):    # Chain id check
            residue.append(atom)
        else:
            yield Residue(residue)
            residue = [atom]
    yield Residue(residue)


def residueDistance2(r1, r2):
    return min([FUNC.distance2(i, j) for i in r1 for j in r2])


def breaks(residuelist, selection=("N", "CA", "C"), cutoff=2.5):
    # Extract backbone atoms coordinates
    bb = [[atom[4:] for atom in residue if atom[0] in selection] for residue in residuelist]
    # Needed to remove waters residues from mixed residues.
    bb = [res for res in bb if res != []]

    # We cannot rely on some standard order for the backbone atoms.
    # Therefore breaks are inferred from the minimal distance between
    # backbone atoms from adjacent residues.
    return [i+1 for i in range(len(bb)-1) if residueDistance2(bb[i], bb[i+1]) > cutoff]


def contacts(atoms, cutoff=5):
    rla = range(len(atoms))
    crd = [atom[4:] for atom in atoms]
    return [(i, j) for i in rla[:-1] for j in rla[i+1:]
            if FUNC.distance2(crd[i], crd[j]) < cutoff]


def add_dummy(beads, dist=0.11, n=2):
    # Generate a random vector in a sphere of -1 to +1, to add to the bead position
    v    = [random.random()*2.-1, random.random()*2.-1, random.random()*2.-1]
    # Calculated the length of the vector and divide by the final distance of the dummy bead
    norm_v = FUNC.norm(v)/dist
    # Resize the vector
    vn   = [i/norm_v for i in v]
    # m sets the direction of the added vector, currently only works when adding one or two beads.
    m = 1
    for j in range(n):
        newName = 'SCD'
        newBead = (newName, tuple([i+(m*j) for i, j in zip(beads[-1][1], vn)]), beads[-1][2])
        beads.append(newBead)
        m *= -2
    return beads


def check_merge(chains, m_list=[], l_list=[], ss_cutoff=0):
    chainIndex = range(len(chains))

    if 'all' in m_list:
        logging.info("All chains will be merged in a single moleculetype.")
        return chainIndex, [chainIndex]

    chainID = [chain.id for chain in chains]

    # Mark the combinations of chains that need to be merged
    merges = []
    if m_list:
        # Build a dictionary of chain IDs versus index
        # To give higher priority to top chains the lists are reversed
        # before building the dictionary
        chainIndex.reverse()
        chainID.reverse()
        dct = dict(zip(chainID, chainIndex))
        chainIndex.reverse()
        # Convert chains in the merge_list to numeric, if necessary
        # NOTE The internal numbering is zero-based, while the
        # command line chain indexing is one-based. We have to add
        # one to the number in the dictionary to bring it on par with
        # the numbering from the command line, but then from the
        # result we need to subtract one again to make indexing
        # zero-based
        merges = [[(i.isdigit() and int(i) or dct[i]+1)-1 for i in j] for j in m_list]
        for i in merges:
            i.sort()

    # Rearrange merge list to a list of pairs
    pairs = [(i[j], i[k]) for i in merges for j in range(len(i)-1) for k in range(j+1, len(i))]

    # Check each combination of chains for connections based on
    # ss-bridges, links and distance restraints
    for i in chainIndex[:-1]:
        for j in chainIndex[i+1:]:
            if (i, j) in pairs:
                continue
            # Check whether any link links these two groups
            for a, b in l_list:
                if ((a in chains[i] and b in chains[j]) or (a in chains[j] and b in chains[i])):
                    logging.info("Merging chains %d and %d to allow link %s" % (i+1, j+1, str((a, b))))
                    pairs.append(i < j and (i, j) or (j, i))
                    break
            if (i, j) in pairs:
                continue
            # Check whether any cystine bond given links these two groups
            #for a,b in s_list:
            #    if ((a in chains[i] and b in chains[j]) or
            #        (a in chains[j] and b in chains[i])):
            #        logging.info("Merging chains %d and %d to allow cystine bridge"%(i+1,j+1))
            #        pairs.append( i<j and (i,j) or (j,i) )
            #        break
            #if (i,j) in pairs:
            #    continue
            # Check for cystine bridges based on distance
            if not ss_cutoff:
                continue
            # Get SG atoms from cysteines from either chain
            # Check this pair of chains
            for cysA in chains[i]["CYS"]:
                for cysB in chains[j]["CYS"]:
                    d2 = FUNC.distance2(cysA["SG"][4:7], cysB["SG"][4:7])
                    if d2 <= ss_cutoff:
                        logging.info("Found SS contact linking chains %d and %d (%f nm)" % (i+1, j+1, math.sqrt(d2)/10))
                        pairs.append((i, j))
                    break
                if (i, j) in pairs:
                    break

    # Sort the combinations
    pairs.sort(reverse=True)

    merges = []
    while pairs:
        merges.append(set([pairs[-1][0]]))
        for i in range(len(pairs)-1, -1, -1):
            if pairs[i][0] in merges[-1]:
                merges[-1].add(pairs.pop(i)[1])
            elif pairs[i][1] in merges[-1]:
                merges[-1].add(pairs.pop(i)[0])
    merges = [list(i) for i in merges]
    for i in merges:
        i.sort()

    order = [j for i in merges for j in i]

    if merges:
        logging.warning("Merging chains.")
        logging.warning("This may change the order of atoms and will change the number of topology files.")
        logging.info("Merges: " + ", ".join([str([j+1 for j in i]) for i in merges]))

    if len(merges) == 1 and len(merges[0]) > 1 and set(merges[0]) == set(chainIndex):
        logging.info("All chains will be merged in a single moleculetype")

    # Determine the order for writing; merged chains go first
    merges.extend([[j] for j in chainIndex if j not in order])
    order.extend([j for j in chainIndex if j not in order])

    return order, merges


## !! NOTE !! ##
## XXX The chain class needs to be simplified by extracting things to separate functions/classes
class Chain:
    # Attributes defining a chain
    # When copying a chain, or slicing, the attributes in this list have to
    # be handled accordingly.
    _attributes = ("residues", "sequence", "seq", "ss", "ssclass", "sstypes")

    def __init__(self, options, residuelist=[], name=None, multiscale=False):
        self.residues   = residuelist
        self._atoms     = [atom[:3] for residue in residuelist for atom in residue]
        self.sequence   = [residue[0][1] for residue in residuelist]
        # *NOTE*: Check for unknown residues and remove them if requested
        #         before proceeding.
        self.seq        = "".join([MAP.AA321.get(i, "X") for i in self.sequence])
        self.ss         = ""
        self.ssclass    = ""
        self.sstypes    = ""
        self.mapping    = []
        self.multiscale = multiscale
        self.options    = options

        # Unknown residues
        self.unknowns   = "X" in self.seq

        # Determine the type of chain
        self._type      = ""
        self.type()

        # Determine number of atoms
        self.natoms     = len(self._atoms)

        # BREAKS: List of indices of residues where a new fragment starts
        # Only when polymeric (protein, DNA, RNA, ...)
        # For now, let's remove it for the Nucleic acids...
        self.breaks     = self.type() in ("Protein", "Mixed") and breaks(self.residues) or []

        # LINKS:  List of pairs of pairs of indices of linked residues/atoms
        # This list is used for cysteine bridges and peptide bonds involving side chains
        # The list has items like ((#resi1, #atid1), (#resi2, #atid2))
        # When merging chains, the residue number needs ot be update, but the atom id
        # remains unchanged.
        # For the coarse grained system, it needs to be checked which beads the respective
        # atoms fall in, and bonded terms need to be added there.
        self.links      = []

        # Chain identifier; try to read from residue definition if no name is given
        self.id         = name or residuelist and residuelist[0][0][3] or ""

        # Container for coarse grained beads
        self._cg        = None

    def __len__(self):
        # Return the number of residues
        # DNA/RNA contain non-CAP d/r to indicate type. We remove those first.
        return len(''.join(i for i in self.seq if i.isupper()))

    def __add__(self, other):
        newchain = Chain(name=self.id+"+"+other.id)
        # Combine the chain items that can be simply added
        for attr in self._attributes:
            setattr(newchain, attr, getattr(self, attr) + getattr(other, attr))
        # Set chain items, shifting the residue numbers
        shift  = len(self)
        newchain.breaks     = self.breaks + [shift] + [i+shift for i in other.breaks]
        newchain.links      = self.links + [((i[0]+shift, i[1]), (j[0]+shift, j[1])) for i, j in other.links]
        newchain.natoms     = len(newchain.atoms())
        newchain.multiscale = self.multiscale or other.multiscale
        # Return the merged chain
        return newchain

    def __eq__(self, other):
        return (self.seq        == other.seq    and
                self.ss         == other.ss     and
                self.breaks     == other.breaks and
                self.links      == other.links  and
                self.multiscale == other.multiscale)

    # Extract a residue by number or the list of residues of a given type
    # This facilitates selecting residues for links, like chain["CYS"]
    def __getitem__(self, other):
        if type(other) == str:
            if other not in self.sequence:
                return []
            return [i for i in self.residues if i[0][1] == other]
        elif type(other) == tuple:
            # This functionality is set up for links
            # between coarse grained beads. So these are
            # checked first,
            for i in self.cg():
                if other == i[:4]:
                    return i
            else:
                for i in self.atoms():
                    if other[:3] == i[:3]:
                        return i
                else:
                    return []
        return self.sequence[other]

    # Extract a piece of a chain as a new chain
    def __getslice__(self, i, j):
        newchain = Chain(self.options, name=self.id)
        # Extract the slices from all lists
        for attr in self._attributes:
            setattr(newchain, attr, getattr(self, attr)[i:j])
        # Breaks that fall within the start and end of this chain need to be passed on.
        # Residue numbering is increased by 20 bits!!
        # XXX I don't know if this works.
        ch_sta, ch_end      = newchain.residues[0][0][2], newchain.residues[-1][0][2]
        newchain.breaks     = [crack for crack in self.breaks if ch_sta < (crack << 20) < ch_end]
        newchain.links      = [link for link in self.links if ch_sta < (link << 20) < ch_end]
        newchain.multiscale = self.multiscale
        newchain.natoms     = len(newchain.atoms())
        newchain.type()
        # Return the chain slice
        return newchain

    def _contains(self, atomlist, atom):
        atnm, resn, resi, chn = atom

        # If the chain does not match, bail out
        if chn != self.id:
            return False

        # Check if the whole tuple is in
        if atnm and resn and resi:
            return (atnm, resn, resi) in self.atoms()

        # Fetch atoms with matching residue id
        match = (not resi) and atomlist or [j for j in atomlist if j[2] == resi]
        if not match:
            return False

        # Select atoms with matching residue name
        match = (not resn) and match or [j for j in match if j[1] == resn]
        if not match:
            return False

        # Check whether the atom is given and listed
        if not atnm or [j for j in match if j[0] == atnm]:
            return True

        # It just is not in the list!
        return False

    def __contains__(self, other):
        return self._contains(self.atoms(), other) or self._contains(self.cg(), other)

    def __hash__(self):
        return id(self)

    def atoms(self):
        if not self._atoms:
            self._atoms = [atom[:3] for residue in self.residues for atom in residue]
        return self._atoms

    # Split a chain based on residue types; each subchain can have only one type
    def split(self):
        chains = []
        chainStart = 0
        for i in range(len(self.sequence)-1):
            if MAP.residueTypes.get(self.sequence[i], "Unknown") != MAP.residueTypes.get(self.sequence[i+1], "Unknown"):
                # Use the __getslice__ method to take a part of the chain.
                chains.append(self[chainStart:i+1])
                chainStart = i+1
        if chains:
            logging.debug('Splitting chain %s in %s chains' % (self.id, len(chains)+1))
        return chains + [self[chainStart:]]

    def getname(self, basename=None):
        name = []
        if basename:                      name.append(basename)
        if self.type() and not basename:  name.append(self.type())
        if type(self.id) == int:
            name.append(chr(64+self.id))
        elif self.id.strip():
            name.append(str(self.id))
        return "_".join(name)

    def set_ss(self, ss, source="self"):
        if len(ss) == 1:
            self.ss = len(self)*ss
        else:
            self.ss = ss
        # Infer the Martini backbone secondary structure types
        self.ssclass, self.sstypes = SS.ssClassification(self.ss, source)

    def dss(self, method=None, executable=None):
        # The method should take a list of atoms and return a
        # string of secondary structure classifications
        if self.type() == "Protein":
            if method:
                atomlist = [atom for residue in self.residues for atom in residue]
                self.set_ss(SS.ssDetermination[method](self, atomlist, executable), source=method)
            else:
                self.set_ss(len(self)*"C")
        else:
            self.set_ss(len(self.sequence)*"-")
        return self.ss

    def type(self, other=None):
        if other:
            self._type = other
        elif not self._type and len(self):
            # Determine the type of chain
            self._type     = set([MAP.residueTypes.get(i, "Unknown") for i in set(self.sequence)])
            self._type     = len(self._type) > 1 and "Mixed" or list(self._type)[0]
        return self._type

    # XXX The following (at least the greater part of it) should be made a separate function, put under "MAPPING"
    def cg(self, force=False, com=False):
        # Generate the coarse grained structure
        # Set the b-factor field to something that reflects the secondary structure

        # If the coarse grained structure is set already, just return,
        # unless regeneration is forced.
        if self._cg and not force:
            return self._cg
        self._cg = []
        atid     = 1
        bb       = [1]
        fail     = False
        previous = ''
        for residue, rss, resname in zip(self.residues, self.sstypes, self.sequence):
            # For DNA we need to get the O3' to the following residue when calculating COM
            # The force and com options ensure that this part does not affect itp generation or anything else
            if com:
                # Just an initialization, this should complain if it isn't updated in the loop
                store = 0
                for ind, i in enumerate(residue):
                    if i[0] == "O3'":
                        if previous != '':
                            residue[ind] = previous
                            previous = i
                        else:
                            store = ind
                            previous = i
                # We couldn't remove the O3' from the 5' end residue during the loop so we do it now
                if store > 0:
                    del residue[store]

            # Check if residues names has changed, for example because user has set residues interactively.
            residue = [(atom[0], resname)+atom[2:] for atom in residue]
            if residue[0][1] in ("SOL", "HOH", "TIP"):
                continue
            if not residue[0][1] in MAP.CoarseGrained.mapping.keys():
                logging.warning("Skipped unknown residue %s\n" % residue[0][1])
                continue
            # Get the mapping for this residue
            # CG.map returns bead coordinates and mapped atoms
            # This will fail if there are (too many) atoms missing, which is
            # only problematic if a mapped structure is written; the topology
            # is inferred from the sequence. So this is the best place to raise
            # an error
            try:
                beads, ids = MAP.map(residue, ca2bb=self.options['ForceField'].ca2bb)
                beads      = zip(MAP.CoarseGrained.names[residue[0][1]], beads, ids)
                if residue[0][1] in self.options['ForceField'].polar:
                    beads = add_dummy(beads, dist=0.14, n=2)
                elif residue[0][1] in self.options['ForceField'].charged:
                    beads = add_dummy(beads, dist=0.11, n=1)
            except ValueError:
                logging.error("Too many atoms missing from residue %s %d(ch:%s):",
                              residue[0][1], residue[0][2]-(32 << 20), residue[0][3])
                logging.error(repr([i[0] for i in residue]))
                fail = True

            for name, (x, y, z), ids in beads:
                # Add the bead with coordinates and secondary structure id to the list
                self._cg.append((name, residue[0][1][:3], residue[0][2], residue[0][3], x, y, z, SS.ss2num[rss]))
                # Add the ids to the list, after converting them to indices to the list of atoms
                self.mapping.append([atid+i for i in ids])

            # Increment the atom id; This pertains to the atoms that are included in the output.
            atid += len(residue)

            # Keep track of the numbers for CONECTing
            bb.append(bb[-1]+len(beads))

        if fail:
            logging.error("Unable to generate coarse grained structure due to missing atoms.")
            sys.exit(1)

        return self._cg

    def conect(self):
        # Return pairs of numbers that should be CONECTed
        # First extract the backbone IDs
        cg = self.cg()
        bb = [i+1 for i, j in zip(range(len(cg)), cg) if j[0] == "BB"]
        bb = zip(bb, bb[1:]+[len(bb)])
        # Set the backbone CONECTs (check whether the distance is consistent with binding)
        conect = [(i, j) for i, j in bb[:-1] if FUNC.distance2(cg[i-1][4:7], cg[j-1][4:7]) < 14]
        # Now add CONECTs for sidechains
        for i, j in bb:
            nsc = j-i-1
