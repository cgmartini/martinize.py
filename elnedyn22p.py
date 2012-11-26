################################
## 6 # FORCE FIELD PARAMETERS ##  -> @FF <-
################################

class elnedyn22p:
    def __init__(self):
        import SS,FUNC,IO 

        # parameters are defined here for the following (protein) forcefields:
        self.name = 'elnedyn22p'
        
        # Charged types:
        self.charges = {"Qd":1, "Qa":-1, "SQd":1, "SQa":-1, "RQd":1, "AQa":-1}                                                           #@#
        
        
        #----+---------------------+
        ## A | BACKBONE PARAMETERS |
        #----+---------------------+
        #
        # bbss  lists the one letter secondary structure code
        # bbdef lists the corresponding default backbone beads
        # bbtyp lists the corresponding residue specific backbone beads
        #
        # bbd   lists the structure specific backbone bond lengths
        # bbkb  lists the corresponding bond force constants
        #
        # bba   lists the structure specific angles
        # bbka  lists the corresponding angle force constants
        #
        # bbd   lists the structure specific dihedral angles
        # bbkd  lists the corresponding force constants
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
        ## BEADS ##                                                                          #                 
        #                               F     E     H     1     2     3     T     S     C    # SS one letter   
        self.bbdef    =     FUNC.spl(" N0   Nda    N0    Nd    Na   Nda   Nda    P5    P5")  # Default beads   #@#
        self.bbtyp    = {                                                                    #                 #@#
                     "ALA": FUNC.spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4"), # ALA specific    #@#
                     "PRO": FUNC.spl(" C5    N0    C5    N0    Na    N0    N0    P4    P4"), # PRO specific    #@#
                     "HYP": FUNC.spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4")  # HYP specific    #@#
        }                                                                                    #                 #@#
        ## BONDS ##                                                                          #                 
        self.bbldef   =             (.365, .350, .350, .350, .350, .350, .350, .350, .350)   # BB bond lengths #@#
        self.bbkb     =             (1250, 1250, 1250, 1250, 1250, 1250,  500,  400,  400)   # BB bond kB      #@#
        self.bbltyp   = {}                                                                   #                 #@#
        self.bbkbtyp  = {}                                                                   #                 #@#
        ## ANGLES ##                                                                         #                 
        self.bbadef   =             (119.2, 134,   96,   96,   96,   96,  100,  130,  127)   # BBB angles      #@#
        self.bbka     =             ( 150,   25,  700,  700,  700,  700,   25,   25,   25)   # BBB angle kB    #@#
        self.bbatyp   = {                                                                    #                 #@#
                    "PRO":          ( 119.2,134,   98,   98,   98,   98,  100,  130,  127),  # PRO specific    #@#
                    "HYP":          ( 119.2,134,   98,   98,   98,   98,  100,  130,  127)   # PRO specific    #@#
        }                                                                                    #                 #@#
        self.bbkatyp  = {                                                                    #                 #@#
                    "PRO":          ( 150,   25,  100,  100,  100,  100,   25,   25,   25),  # PRO specific    #@#
                    "HYP":          ( 150,   25,  100,  100,  100,  100,   25,   25,   25)   # PRO specific    #@#
        }                                                                                    #                 #@#
        ## DIHEDRALS ##                                                                      #                 
        self.bbddef   =             (90.7,    0, -120, -120, -120, -120)                     # BBBB dihedrals  #@#
        self.bbkd     =             ( 100,   10,  400,  400,  400,  400)                     # BBBB kB         #@#
        self.bbdmul   =             (   1,    1,    1,    1,    1,    1)                     # BBBB mltplcty   #@#
        self.bbdtyp   = {}                                                                   #                 #@#
        self.bbkdtyp  = {}                                                                   #                 #@#
                                                                                             #                 
        ###############################################################################################               
        
        # Some Forcefields use the Ca position to position the BB-bead (me like!)
        self.ca2bb = True 
        
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
        #RES#   BEADS                      BONDS                                                                                      ANGLES                          DIHEDRALS  V-SITES
        'TRP': [FUNC.spl("SC4 SNd SC5 SC5"), [(0.255,73000), (0.220,None), (0.250,None), (0.280,None), (0.255,None), (0.35454,None)], [(142,30), (143,20), (104,50)], [(180,200)]],
        'TYR': [FUNC.spl("SC4 SC4 SP1"),     [(0.335, 6000), (0.335,6000), (0.240,None), (0.310,None), (0.310,None)],                 [(70,100), (130, 50)]],
        'PHE': [FUNC.spl("SC5 SC5 SC5"),     [(0.340, 7500), (0.340,7500), (0.240,None), (0.240,None), (0.240,None)],                 [(70,100), (125,100)]],
        'HIS': [FUNC.spl("SC4 SP1 SP1"),     [(0.195,94000), (0.193,None), (0.295,None), (0.216,None)],                               [(135,100),(115, 50)]],
        'HIH': [FUNC.spl("SC4 SP1 SQd"),     [(0.195,94000), (0.193,None), (0.295,None), (0.216,None), (0.11,None)],                  [(135,100),(115, 50)]],
        'GLN': [FUNC.spl("Nda D D"),         [(0.300, 2400), (0.280,None)],                                                           [],                             [],         [(0.5,)]],
        'ASN': [FUNC.spl("Nda D D"),         [(0.250,61000), (0.280,None)],                                                           [],                             [],         [(0.5,)]],
        'SER': [FUNC.spl("N0 D D"),          [(0.195,94000), (0.280,None)],                                                           [],                             [],         [(0.5,)]],
        'THR': [FUNC.spl("N0 D D"),          [(0.195,94000), (0.280,None)],                                                           [],                             [],         [(0.5,)]],
        'ARG': [FUNC.spl("N0 Qd D"),         [(0.250,12500), (0.350,6200), (0.110,None)],                                             [(150,15)]],
        'LYS': [FUNC.spl("C3 Qd D"),         [(0.250,12500), (0.300,9700), (0.110,None)],                                             [(150,20)]],
        'ASP': [FUNC.spl("Qa D"),            [(0.255, None), (0.110,None)]],
        'GLU': [FUNC.spl("Qa D"),            [(0.310, 2500), (0.110,None)]],
        'CYS': [FUNC.spl("C5"),              [(0.240, None)]],
        'ILE': [FUNC.spl("C1"),              [(0.225,13250)]],
        'LEU': [FUNC.spl("C1"),              [(0.265, None)]],
        'MET': [FUNC.spl("C5"),              [(0.310, 2800)]],
        'PRO': [FUNC.spl("C3"),              [(0.190, None)]],
        'HYP': [FUNC.spl("P1"),              [(0.190, None)]],
        'VAL': [FUNC.spl("C2"),              [(0.200, None)]],
        'GLY': [],
        'ALA': [],
        }
        
        # Not all (eg Elnedyn) forcefields use backbone-backbone-sidechain angles and BBBB-dihedrals.
        self.UseBBSAngles        = False 
        self.UseBBBBDihedrals    = False

        # Martini 2.2p has polar and charged residues with seperate charges.
        self.polar   = ["GLN","ASN","SER","THR"]
        self.charged = ["ARG","LYS","ASP","GLU","HIH"]

        # If masses or charged diverge from standard (45/72 and -/+1) they are defined here.
        self.mass_charge = {
        #RES   MASS               CHARGE
        "GLN":[[0,36,36],         [0,0.42,-0.42]],
        "ASN":[[0,36,36],         [0,0.46,-0.46]],
        "SER":[[0,36,36],         [0,0.40,-0.40]],
        "THR":[[0,36,36],         [0,0.36,-0.36]],
        "ARG":[[72,36,36],        [0,0,1]],
        "LYS":[[72,36,36],        [0,0,1]],
        "HIH":[[72,72,36,36],    [0,0,0,1]],
        "ASP":[[36,36],           [0,-1]],
        "GLU":[[36,36],           [0,-1]],
        }

        # Defines the connectivity between between beads
        # The polar sidechains have charged dummy beads, connected with a constraint
        # The charged sidechains have a charged dummy bead.
        self.connectivity = {
        #RES       BONDS                                              ANGLES                            DIHEDRALS       V-SITE
        "TRP":     [[(0, 1), (1, 2), (2, 4), (4, 3), (3, 1), (1, 4)],[(0, 1, 2), (0, 1, 4), (0, 1, 3)],[(1, 2, 3, 4)]],
        "TYR":     [[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "PHE":     [[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)],        [(0, 1, 2), (0, 1, 3)]],
        "HIS":     [[(0, 1), (1, 2), (1, 3), (2, 3)],                [(0, 1, 2), (0, 1, 3)]],
        "HIH":     [[(0, 1), (1, 2), (1, 3), (2, 3), (3, 4)],        [(0, 1, 2), (0, 1, 3)],           [(0, 2, 3, 1)]],
        "GLN":     [[(0, 1), (2, 3)],                                [],                               [],              [(1,2,3)]],
        "ASN":     [[(0, 1), (2, 3)],                                [],                               [],              [(1,2,3)]],
        "SER":     [[(0, 1), (2, 3)],                                [],                               [],              [(1,2,3)]],
        "THR":     [[(0, 1), (2, 3)],                                [],                               [],              [(1,2,3)]],
        "ARG":     [[(0, 1), (1, 2), (2, 3)],                        [(0,1,2)]],
        "LYS":     [[(0, 1), (1, 2), (2, 3)],                        [(0,1,2)]],
        "ASP":     [[(0, 1), (1, 2)]],
        "GLU":     [[(0, 1), (1, 2)]],
        "CYS":     [[(0, 1)]],
        "ILE":     [[(0, 1)]],
        "LEU":     [[(0, 1)]],
        "MET":     [[(0, 1)]],
        "PRO":     [[(0, 1)]],
        "HYP":     [[(0, 1)]],
        "VAL":     [[(0, 1)]],
        "ALA":     [],
        "GLY":     [],
        }
       
        #----+----------------+
        ## C | SPECIAL BONDS  |
        #----+----------------+
        
        self.special = {
            # Used for sulfur bridges
            # ATOM 1         ATOM 2          BOND LENGTH   FORCE CONSTANT
            (("SC1","CYS"), ("SC1","CYS")):     (0.39,         5000),
            }
       
        # By default use an elastic network
        self.ElasticNetwork = True 

        # Elastic networks bond shouldn't lead to exclusions (type 6) 
        # But Elnedyn has been parametrized with type 1.
        self.EBondType = 1
        
        #----+----------------+
        ## D | INTERNAL STUFF |
        #----+----------------+
        
        
        ## BACKBONE BEAD TYPE ##                                                                    
        # Dictionary of default bead types (*D)                                                     
        self.bbBeadDictD  = FUNC.hash(SS.bbss,self.bbdef)                                                             
        # Dictionary of dictionaries of types for specific residues (*S)                            
        self.bbBeadDictS  = dict([(i,FUNC.hash(SS.bbss,self.bbtyp[i])) for i in self.bbtyp.keys()])                        
         
        ## BB BOND TYPE ##                                                                          
        # Dictionary of default abond types (*D)                                                    
        self.bbBondDictD = FUNC.hash(SS.bbss,zip(self.bbldef,self.bbkb))                                                   
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbBondDictS = dict([(i,FUNC.hash(SS.bbss,zip(self.bbltyp[i],self.bbkbtyp[i]))) for i in self.bbltyp.keys()])       
        # This is tricky to read, but it gives the right bondlength/force constant

        ## BBB ANGLE TYPE ##                                                                        
        # Dictionary of default angle types (*D)                                                    
        self.bbAngleDictD = FUNC.hash(SS.bbss,zip(self.bbadef,self.bbka))                                                  
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbAngleDictS = dict([(i,FUNC.hash(SS.bbss,zip(self.bbatyp[i],self.bbkatyp[i]))) for i in self.bbatyp.keys()])      
       
        ## BBBB DIHEDRAL TYPE ##                                                                    
        # Dictionary of default dihedral types (*D)                                                 
        self.bbDihedDictD = FUNC.hash(SS.bbss,zip(self.bbddef,self.bbkd,self.bbdmul))                                           
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbDihedDictS = dict([(i,FUNC.hash(SS.bbss,zip(self.bbdtyp[i],self.bbkdtyp[i]))) for i in self.bbdtyp.keys()])      

    # The following function returns the backbone bead for a given residue and                   
    # secondary structure type.                                                                 
    # 1. Look up the proper dictionary for the residue                                          
    # 2. Get the proper type from it for the secondary structure                                
    # If the residue is not in the dictionary of specials, use the default                      
    # If the secondary structure is not listed (in the residue specific                         
    # dictionary) revert to the default.                                                        
    def bbGetBead(self,r1,ss="C"):                                                                   
        return self.bbBeadDictS.get(r1,self.bbBeadDictD).get(ss,self.bbBeadDictD.get(ss))                      
    
    # For Elnedyn we need something else to get the bond length (much simpler due to Ca position BB's)
    def bbGetBond(self,r,ca,ss):
        import FUNC 
        import math
        # The 150000 forceconstant gave an error message, turning to constraints would be better.
        return ( math.sqrt(FUNC.distance2(ca[0],ca[1]))/10., None   )
    
    def bbGetAngle(self,r,ca,ss):
        import FUNC,IO 
        import math
        # Elnedyn takes angles from structure, with fc=40
        return (math.acos(FUNC.cos_angle([i-j for i,j in zip(ca[0],ca[1])],[i-j for i,j in zip(ca[2],ca[1])]))/IO.d2r, 40)

    def messages(self):
        '''Prints any force-field specific logging messages.'''
        import logging
        logging.info('The elnedyn forcefield has been implemented with some changes compared to the published parameters:')
        #logging.info('- Backbone-Backbone bonds are constraints in stead of high force constant bonds.')
        logging.info('- Backbone-Backbone bonds use high force constant bonds instead of constraints.')
        logging.info('- Trp has an extra constrain added to the sidechain.')
        logging.info('- The Backbone sidechain bonds with high force constants are replaced by constraints except for Trp and His and the polar sidechains.')
        logging.warning('Elnedyn topologies might not give numerical stable simulations with a 20fs timestep.')
        logging.warning('This can be solved by setting all S-type bead masses to 72amu.')
        logging.warning('Martini version 2.2 is in beta release. It has not been extensively tested and problems might occur.')
        pass

