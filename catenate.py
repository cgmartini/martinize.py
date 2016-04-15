def cat():
    '''Function to 'compile' (catenate) the martinize script modules into one file.'''
    import re,os
    from martinize import version
    file_name = "martinize-%s.py"%version

    # Parameters are defined for the following (protein) forcefields:
    forcefields = [ff[:-6] for ff in os.listdir(".") if ff[-6:] == "_ff.py"]
    if os.path.isdir("ForceFields"):
        forcefields += ["ForceFields/"+ff[:-6] for ff in os.listdir("ForceFields") if ff[-6:] == "_ff.py"]

    # Not all forcefields should be included.
    print "Found %s forcefields:"%len(forcefields)
    for i,ff in enumerate(forcefields):
        print '%s. %s'%(i,ff)
    forcefields = [forcefields[i] for i in input("Which ones should be included? (enter comma seperate string):")]
 
    files_in = 'martinize.py '+'_ff.py '.join(forcefields)+'_ff.py DOC.py CMD.py FUNC.py MAP.py SS.py ELN.py IO.py TOP.py MAIN.py '
    pattern1 = re.compile(files_in.replace('.py ','|')[:-1])
    pattern2 = re.compile(files_in.replace('.py ','\.|')[:-1])
    file_out = open(file_name,'w')
    tail = ''; head = True
    for f in files_in.split():
        for line in open(f).readlines():
            # split the variable to avoid finding it self.
            if "__na"+"me__" in line:
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
    os.chmod(file_name,0755)

if __name__ == '__main__':
    cat()
