################################
## 6 # FORCE FIELD PARAMETERS ##  -> @FF <-
################################

class polbb:
    '''The forcefield has been implemented with some changes compared to the published parameters:
    - Backbone-Backbone bonds are constraints in stead of strong bonds.
    - Trp has an extra constrain added to the sidechain
    - The Backbone-Sidechain bonds with high force constants are replaced by constraints except for Trp and His.
    '''
    def __init__(self):
        import SS,FUNC,IO 

        # parameters are defined here for the following (protein) forcefields:
        self.name = 'polbb'
        
        # Charged types:
        self.charges      = {"Qd":1, "Qa":-1, "SQd":1, "SQa":-1, "RQd":1, "AQa":-1}                              #@#
        self.dummycharges = {"DIM":-0.40, "DIP":0.40,}                                                           #@#
        
        #----+---------------------+
        ## A | BACKBONE PARAMETERS |
        #----+---------------------+
        #
        # -=NOTE=- 
        #  if the secondary structure types differ between bonded atoms
        #  the bond is assigned the lowest corresponding force constant 
        #
        # -=NOTE=-
        # if proline is anywhere in the helix, the BBB angle changes for 
        # all residues
        #
        
        ###############################################################################################
        #                               F     E     H     1     2     3     T     S     C    # SS one letter   
        ## BEADS ##                                                                          #                 
        ## BONDS ##                                                                          #                 
        ## PAIRS ##                                                                          #
        self.bbpk     = {
                      'F':  FUNC.spl("0.40 0.30  0.90  0.90  0.90  0.90  0.30  0.30  0.20"), # Default beads   #@#
                      'E':  FUNC.spl("0.40 0.30  0.90  0.90  0.90  0.90  0.30  0.30  0.20"), # Default beads   #@#
                      'H':  FUNC.spl("0.90 0.90  0.90  0.90  0.90  0.90  0.90  0.90  0.90"), # Default beads   #@#
                      '1':  FUNC.spl("0.90 0.90  0.90  0.90  0.90  0.90  0.90  0.90  0.90"), # Default beads   #@#
                      '2':  FUNC.spl("0.90 0.90  0.90  0.90  0.90  0.90  0.90  0.90  0.90"), # Default beads   #@#
                      '3':  FUNC.spl("0.90 0.90  0.90  0.90  0.90  0.90  0.90  0.90  0.90"), # Default beads   #@#
                      'T':  FUNC.spl("0.40 0.30  0.90  0.90  0.90  0.90  0.30  0.30  0.20"), # Default beads   #@#
                      'S':  FUNC.spl("0.40 0.30  0.90  0.90  0.90  0.90  0.30  0.30  0.20"), # Default beads   #@#
                      'C':  FUNC.spl("0.40 0.30  0.90  0.90  0.90  0.90  0.30  0.30  0.20"), # Default beads   #@#
                        }
        ## ANGLES ##                                                                         #                 
        ## DIHEDRALS ##                                                                      #                 
        ###############################################################################################               
        
        # In a polarizable backbone there are 4 beads in the backbone. 
        # Parameters are dependent on which bead it is.
        # Only the pair interactions are dependent on secondary structure
        self.pol_bb = {
        #                        BAS DIM DBB DIP   
            'atom'   : FUNC.spl("BAS DIM DBB DIP"),
            'bond'   : [(0.215,50000),(0.236,50000),(0.249,50000),(0.206,50000),(0.249,50000),(0.380,None),],
            'angle'  : [(8,0,2.0),(1,95,1.5),],
            'dih'    : [(8,1,0.2),(8,5,0.8),],
            'imp'    : [(2,0,100),],
            'vsite'  : [(2,4,0.5),],
            'excl'   : [(6,8),(6,8)],
            'pair'   : [(-0.40,-0.40,0.0,1e-7),(-0.40,0.40,0.0,1e-7),(0.40,-0.40,0.0,1e-7),(0.40,0.40,0.0,1e-7)]
        }
        # connectivity using atoms number in first residue.
        # First bead is zero, to make 'remainder' calculation possible.
        self.pol_con  = {
            'bond'   : [(0,1),(0,3),(1,4),(3,4),(1,3),(0,4)],
            'angle'  : [(0,4,8),(3,4,8)],
            'dih'    : [(1,3,4,8),(1,3,5,7)],
            'imp'    : [(0,4,1,3)],
            'vsite'  : [(2,)],
            'excl'   : [(1,5,7),(3,5,7)],
            'pair'   : [(1,5),(1,7),(3,5),(3,7)],
        }

        # Some Forcefields use the Ca position to position the BB-bead (me like!)
        # Setting puts extra beads in the coarse grain
        self.ca2bb = 'polBB' 
        
        # BBS angle, equal for all ss types                                                         
        # Connects BB(i-1),BB(i),SC(i), except for first residue: BB(i+1),BB(i),SC(i)               
        #                      ANGLE   Ka                                                                
        self.bbsangle =      [   100,  25]                                                          #@#
        
        # Bonds for extended structures (more stable than using dihedrals)                          
        #               LENGTH FORCE                                                                
        self.ebonds   = {                                                                           #@#
               'short': [ .640, 2500],                                                              #@#
               'long' : [ .970, 2500]                                                               #@#
        }                                                                                           #@#
        #----+-----------------------+
        ## B | SIDE CHAIN PARAMETERS |
        #----+-----------------------+
        
        # Sidechain parameters for Elnedyn. (read from cg-2.1.dat). 
        # For HIS the order of bonds is changed and a bond with fc=0 is added.
        # In the elnedyn2, TRP has an extra, cross-ring constraint
        self.sidechains = {
        #RES#   BEADS                      BONDS                                                                    ANGLES                          DIHEDRALS
        'TRP': [FUNC.spl("SC4 SNd SC5 SC5"), [(0.255,73000), (0.220,None), (0.250,None), (0.280,None), (0.255,None), (0.35454,None)], [(142,30), (143,20), (104,50)], [(180,200)]],
        'TYR': [FUNC.spl("SC4 SC4 SP1"),     [(0.335, 6000), (0.335,6000), (0.240,None), (0.310,None), (0.310,None)], [(70,100), (130, 50)]],
        'PHE': [FUNC.spl("SC5 SC5 SC5"),     [(0.340, 7500), (0.340,7500), (0.240,None), (0.240,None), (0.240,None)], [(70,100), (125,100)]],
        'HIS': [FUNC.spl("SC4 SP1 SP1"),     [(0.195,94000), (0.193,None), (0.295,None), (0.216,None)],               [(135,100),(115, 50)]],
        'HIH': [FUNC.spl("SC4 SP1 SQd"),     [(0.195,94000), (0.193,None), (0.295,None), (0.216,None)],               [(135,100),(115, 50)]],
        'ARG': [FUNC.spl("N0 Qd"),           [(0.250,12500), (0.350,6200)],                                           [(150,15)]],
        'LYS': [FUNC.spl("C3 Qd"),           [(0.250,12500), (0.300,9700)],                                           [(150,20)]],
        'CYS': [FUNC.spl("C5"),              [(0.240, None)]],
        'ASP': [FUNC.spl("Qa"),              [(0.255, None)]],
        'GLU': [FUNC.spl("Qa"),              [(0.310, 2500)]],
        'ILE': [FUNC.spl("C1"),              [(0.225,13250)]],
        'LEU': [FUNC.spl("C1"),              [(0.265, None)]],
        'MET': [FUNC.spl("C5"),              [(0.310, 2800)]],
        'ASN': [FUNC.spl("P5"),              [(0.250, None)]],
        'PRO': [FUNC.spl("C3"),              [(0.190, None)]],
        'GLN': [FUNC.spl("P4"),              [(0.300, 2400)]],
        'SER': [FUNC.spl("P1"),              [(0.195, None)]],
        'THR': [FUNC.spl("P1"),              [(0.195, None)]],
        'VAL': [FUNC.spl("C2"),              [(0.200, None)]],
        'GLY': [],
        'ALA': [],
        }
        
        # Not all (eg Elnedyn) forcefields use backbone-backbone-sidechain angles and BBBB-dihedrals.
        self.UseBBSAngles        = False 
        self.UseBBBBDihedrals    = False

        # Martini 2.2p has polar and charged residues with seperate charges.
        self.polar   = []
        self.charged = []

        # If masses or charged diverge from standard (45/72 and -/+1) they are defined here.
        self.mass_charge = {
        #RES   MASS               CHARGE
        }

        # Defines the connectivity between between beads
        # Connectivity records for Elnedyn (read from cg-2.1.dat). 
        # For HIS the order of bonds is changed and a bond with fc=0 is added.
        self.connectivity = {
        #RES       BONDS                                             ANGLES                            DIHEDRALS       V-SITE
        "TRP":     [[(0, 1), (1, 2), (2, 4), (4, 3), (3, 1), (1, 4)],[(0, 1, 2), (0, 1, 4), (0, 1, 3)],[(1, 2, 3, 4)]],
        "TYR":     [[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "PHE":     [[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "HIS":     [[(0, 1), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "HIH":     [[(0, 1), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "GLN":     [[(0,1)]],
        "ASN":     [[(0,1)]],
        "SER":     [[(0,1)]],
        "THR":     [[(0,1)]],
        "ARG":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "LYS":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "ASP":     [[(0,1)]],
        "GLU":     [[(0,1)]],
        "CYS":     [[(0,1)]],
        "ILE":     [[(0,1)]],
        "LEU":     [[(0,1)]],
        "MET":     [[(0,1)]],
        "PRO":     [[(0,1)]],
        "HYP":     [[(0,1)]],
        "VAL":     [[(0,1)]],
        "ALA":     [],
        "GLY":     [],
        }

        #----+----------------+
        ## C | DNA/RNA bases  |
        #----+----------------+

        #----+----------------+
        ## C | SPECIAL BONDS  |
        #----+----------------+
        
        self.special = {
            # Used for sulfur bridges
            # ATOM 1         ATOM 2          BOND LENGTH   FORCE CONSTANT
            (("SC1","CYS"), ("SC1","CYS")):     (0.39,         5000),
            }
       
        # By default use an elastic network
        self.ElasticNetwork = False

        # Elastic networks bond shouldn't lead to exclusions (type 6) 
        # But Elnedyn has been parametrized with type 1.
        self.EBondType = 1
        
        #----+----------------+
        ## D | INTERNAL STUFF |
        #----+----------------+

        # Dictionary to get the right bond given a connectivity for pol BB.
        # We might not neede the two Dictionaries above anymore in this FF-file.
        self.bbBondDictC = dict(zip(self.pol_con['bond'],self.pol_bb['bond']))

        # Dictionary to get the right angle given a connectivity for pol BB.
        # We might not neede the two Dictionaries above anymore in this FF-file.
        self.bbAngleDictC = dict(zip(self.pol_con['angle'],self.pol_bb['angle']))
       
        # Dictionary to get the right dictionary given a connectivity for pol BB.
        # We might not neede the two Dictionaries above anymore in this FF-file.
        self.bbDihedDictC = dict(zip(self.pol_con['dih']+self.pol_con['imp'],self.pol_bb['dih']+self.pol_bb['imp']))

        ## BB VSITE AND EXCLUSION TYPE ##
        self.bbVsiteDictC = dict(zip(self.pol_con['vsite'],self.pol_bb['vsite']))
        self.bbExclusionDictC = dict(zip(self.pol_con['excl'],self.pol_bb['excl']))

        ## BB PAIR TYPE ## 
        # Dictionary of default dihedral types (*D)                                                 
        self.bbPairDictD = dict([(i,FUNC.hash(SS.bbss,self.bbpk[i])) for i in self.bbpk.keys()])
        self.bbPairDictC = dict(zip(self.pol_con['pair'],self.pol_bb['pair']))
        # Dictionary of dictionaries for specific types (*S)                                        
        # self.bbDihedDictS = dict([(i,FUNC.hash(SS.bbss,zip(self.bbdtyp[i],self.bbkdtyp[i]))) for i in self.bbdtyp.keys()])


    # The following function returns the backbone bead for a given residue and                   
    # secondary structure type.                                                                 
    # 1. Look up the proper dictionary for the residue                                          
    # 2. Get the proper type from it for the secondary structure                                
    # If the residue is not in the dictionary of specials, use the default                      
    # If the secondary structure is not listed (in the residue specific                         
    # dictionary) revert to the default.                                                        
    def bbGetBead(self,r1,ss="C"):                                                                   
        return self.pol_bb['atom']                      
    
    def bbGetBond(self,r,ca,ss):
        # The next line put everything in ca back by the amount needed to set the first feeld to be <5
        ca = tuple([i-(ca[0]-ca[0]%4) for i in ca])
        return ca in self.bbBondDictC.keys() and self.bbBondDictC[ca] or None
    
    def bbGetAngle(self,r,ca,ss):
        # The next line put everything in ca back by the amount needed to set the first feeld to be <5
        ca = tuple([i-(ca[0]-ca[0]%4) for i in ca])
        return ca in self.bbAngleDictC.keys() and self.bbAngleDictC[ca] or None

    def bbGetDihedral(self,r,ca,ss):
        # The next line put everything in ca back by the amount needed to set the first feeld to be <5
        ca = tuple([i-(ca[0]-ca[0]%4) for i in ca])
        return ca in self.bbDihedDictC.keys() and self.bbDihedDictC[ca] or None

    def bbGetVsite(self,r,ca,ss):
        shift = (ca[0]-ca[0]%4) # Used to put return values at right value.
        ca = tuple([i-(ca[0]-ca[0]%4) for i in ca])
        return ca in self.bbVsiteDictC.keys() and tuple([i+shift for i in self.bbVsiteDictC[ca][:-1]]+list(self.bbVsiteDictC[ca][-1:])) or None 

    def bbGetExclusion(self,r,ca,ss):
        shift = (ca[0]-ca[0]%4) # Used to put return values at right value.
        ca = tuple([i-(ca[0]-ca[0]%4) for i in ca])
        return ca in self.bbExclusionDictC.keys() and tuple([i+shift for i in self.bbExclusionDictC[ca]]) or None 

    def bbGetPair(self,r,ca,ss):
        ca = tuple([i-(ca[0]-ca[0]%4) for i in ca])
        # This we have to fancy up, using ss to determine the pair parameters!
        return ca in self.pol_con['pair'] and ((self.bbPairDictD[ss[0]][ss[1]],)+self.bbPairDictC[ca]) or None 

    def getCharge(self,atype,aname):
        return self.charges.get(atype,self.dummycharges.get(aname,0))

    def messages(self):
        '''Prints any force-field specific logging messages.'''
        import logging
        logging.warning('There is no such thing as a Martini polarizable backbone!!!')
        logging.warning('IF YOU ARE READING THIS, YOU SHOULD BE VERY AFFRAID!!!')
        pass
