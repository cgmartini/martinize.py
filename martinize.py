#!/usr/bin/env python


# EDITABLE SECTIONS ARE MARKED WITH #@# 


version="2.5a"
authors=["Djurre de Jong", "Jaakko J. Uusitalo", "Tsjerk A. Wassenaar"]

# Parameters are defined for the following (protein) forcefields:
forcefields = ['martini21','martini21p','martini22','martini22p','elnedyn','elnedyn22','elnedyn22p','martini22dna']


notes = [
    ("DdJ130213","V2.3"),
    ("DdJ200613","Fixes in cysteine bridge detection and help text."),
    ("DdJ200820","Fixes in cysteine bridge length and revert elnedyn BB-bonds to bonds."),
    ("DdJ200826","Inverted 'define NO_RUBBER_BANDS', fixed writing posres when merging and added few comments."),
    ("DdJ200831","Shortened in-file changelog and fixed some comments."),
    ("DdJ181013","V2.4"),
    ]

# 
# This program has grown to be pretty complex. 
# The routines have been organized in different files.
# For working versions, all files can be incorporated by using the option -cat. 
#
# Index of the program files:
#
#   1. Options and documentation                             @DOC.py
#   2. Description, options and command line parsing         @CMD.py
#   3. Helper functions and macros                           @FUNC.py
#   4. Finegrained to coarsegrained mapping                  @MAP.py
#   5. Secondary structure determination and interpretation  @SS.py
#   6. Force field parameters (MARTINI/ELNEDYN)              @FF.py
#   7. Elastic network                                       @ELN.py
#   8. Structure I/O                                         @IO.py
#   9. Topology generation                                   @TOP.py
#  10. Main                                                  @MAIN.py
#  11. Web-interface		                			     @WEB.py
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
    # The argument cat is only given once: when concatenating to on exportable script.
    if '-cat' in args:
        cat('martinize-'+version+'.py')
        sys.exit()
    # Get the possible commandline arguments arguments and help text. 
    options,lists = DOC.options,DOC.lists
    # Parse commandline options.
    options = CMD.option_parser(args,options,lists,version)

    MAIN.main(options)
