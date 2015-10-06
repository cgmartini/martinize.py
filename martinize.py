#!/usr/bin/env python


# EDITABLE SECTIONS ARE MARKED WITH #@#


version = "2.5"
authors = ["Djurre de Jong", "Jaakko J. Uusitalo", "Tsjerk A. Wassenaar"]

notes = [
    ("DdJ181013", "V2.4"),
    ("DdJ041213", "Fixed bug where Cys-Cys constraints were not recognized as such."),
    ("DdJ110815", "Removed warnings about beta status of Martini 2.2."),
    ("DdJ110815", "V2.5"),
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
#


if __name__ == '__main__':
    import sys, logging
    import DOC, CMD, MAIN
    args = sys.argv[1:]
    # Get the possible commandline arguments arguments and help text. 
    options, lists = DOC.options, DOC.lists
    # Parse commandline options.
    options = CMD.option_parser(args, options, lists, version)

    MAIN.main(options)
