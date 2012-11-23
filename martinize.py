#!/usr/bin/env python


# EDITABLE SECTIONS ARE MARKED WITH #@# 


version="2.2"
authors=["Djurre de Jong", "Tsjerk A. Wassenaar"]

# Parameters are defined for the following (protein) forcefields:
forcefields = ['martini21','martini21p','martini22','martini22p','elnedyn','elnedyn22','elnedyn22p']


notes = [
    ("TAW110929","Modified Option class to handle boolean options nicely"),
    ("TAW110929","Modified position restraints to allow arbitrary atom name selections"),
    ("TAW110929","Enabled merging chains (option -merge)"),
    ("TAW110930","Changed order of fields in CG atom representation (name,resn,chain,resi)->(name,resn,resi,chain)"),
    ("TAW110930","Enabled cystine bridges, including automatic detection"),
    ("TAW111007","Added stuff to support multiscaling"),
    ("TAW111007","Added writing of an index file for multiscaling"),
    ("TAW120215","Started rearranging things and added documentation to options"),    
    ("TAW120314","Enabled elastic networks"),
    ("TAW120314","Moved mapping stuff to class CoarseGrained - Class needs to be made into generic converter"),
    ("TAW120329","Added support for residue insertion codes in PDB files"),
    ("DdJ120401","v1.1"),
    ("DdJ120512","Fixed bug with counter in multi chain topologies"),
    ("DdJ120512","Corrected wrong collagen parameters"),
    ("DdJ120521","Corrected wrong collagen parameters"),
    ("DdJ120521","Fixed bug involving BBBB dihedrals in extended regions."),
    ("DdJ120522","v1.2"),
    ("DdJ120406","Splitted script in to 10+1 seperate modules."),
    ("DdJ120419","Added a WEB module."),
    ("DdJ240419","Moved all the forcefield specific settings to FF or CMD, to allow easy adding of forcefields."),
    ("DdJ120512","Every forcefield now has it's own file."),
    ("DdJ140612","The program is martini 2.2(p) ready (dummy beads in pdb, dummies,vsites and p-charged in itp."),
    ("DdJ140612","Added new forcefield files (martini22, martini22p, elnedyn2)"),
    ("DdJ260612","Adopted forcefield class to better fit in polar and charged particles"),
    ("DdJ260612","Added new forcefield files (elnedyn, elnedyn22, elnedyn22p)"),
    ("DdJ120712","Added basic DNA capabilities"),
    ("DdJ240712","Fixed bug with manually defined Cys-bridges due to insertion code (TAW120329)."), 
    ("DdJ240712","Added FF specific message definition."),
    ("DdJ240712","Added exclusions for p-FFs. Changed naming of dummy beads."),
    ("DdJ220812","Bug fixing bond lengths, ss for multiple chains and print message in itp."),
    ("DdJ140912","Fixed wrong combining of position restraints when merging proteins."),
    ("DdJ140912","Fixed bug in TYR/HIS in elnedyn forcefields."),
    ("DdJ140912","Added warning about differing beadnames and printing commandline arguments to itp."),
    ("DdJ270912","The Ca position is now determined based on atom name (CA) iso second atom of residue."),
    ("DdJ280912","Fixed more bugs in the elnedyn definition."),
    ("DdJ221112","Clean help text and some code."),
    ("DdJ221112","Fixed that crashed the break checking code if water chain name == protein chain name."),
    ("DdJ231112","Fixed bug when helix was starting at first residue."),
    ]

# 
# This program has grown to be pretty complete and complex. 
# The routines have been organized in sections, which are 
# tagged to make jumping to a particular section easy.
# For working versions, the sections are in different modules
#
# Index of this file:
#
#   1. Options and documentation                             @DOC
#   2. Description, options and command line parsing         @CMD
#   3. Helper functions and macros                           @FUNC
#   4. Finegrained to coarsegrained mapping                  @MAP
#   5. Secondary structure determination and interpretation  @SS
#   6. Force field parameters (MARTINI/ELNEDYN)              @FF
#   7. Elastic network                                       @ELN
#   8. Structure I/O                                         @IO
#   9. Topology generation                                   @TOP
#  10. Main                                                  @MAIN
#  11. Web-interface		                			     @WEB
#  

def cat(file_out):
    '''Function to 'compile' the martinize script into one file.'''
    import re
    files_in = 'martinize.py DOC.py CMD.py FUNC.py MAP.py SS.py '+'.py '.join(forcefields)+'.py ELN.py IO.py TOP.py MAIN.py '
    pattern1 = re.compile(files_in.replace('.py ','|')[:-1])
    pattern2 = re.compile(files_in.replace('.py ','\.|')[:-1])
    file_out = open(file_out,'w')
    tail = ''; head = True
    for f in files_in.split():
        for line in open(f).readlines():
            # Split the string to avoid the function finding itself
            if '__na'+'me__' in line:
                head = False
            if head:
                file_out.write(line)
            elif (f == 'martinize.py' and not head) and not ('import' in line and pattern1.search(line)):
                tail += pattern2.sub('',line)
            elif line[0] == '#':
                file_out.write(line)
            elif not ('import' in line and pattern1.search(line)):
                file_out.write(pattern2.sub('',line))
    file_out.write(tail)

if __name__ == '__main__':
    import sys,logging
    import DOC,CMD,MAIN
    args = sys.argv[1:]
    if '-cat' in args:
        cat('martinize-'+version+'.py')
        sys.exit()
    options,lists = DOC.options,DOC.lists
    options = CMD.option_parser(args,options,lists,version)

    MAIN.main(options)
