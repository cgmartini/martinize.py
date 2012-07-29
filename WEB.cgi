#!/usr/bin/python
import cgi,cgitb
# cgitb should only be enable in a debug version
cgitb.enable()
print "Content-Type: text/html"
print

########################
## 11 # WEB-INTERFACE ##  -> @WEB <-
########################

# File needs: - A working installation of DSSP
#             - Directory called 'out' at same level, to save output.
#             - Input from a html form.

dssp = '/home/djurre/public_html/martinize/dssp-2-linux-i386'
version = '2.0web'

import os,sys
import urllib2,gzip,StringIO
import DOC,CMD,MAIN,martinize

def restore_std(stdout,stdout_old,stderr,stderr_old):
    # Handle logging output
    sys.stdout = stdout_old        # Restore the old stdout
    sys.stderr = stderr_old        # Restore the old stdout
    
    for line in stderr.getvalue().split('\n'):
        print line,'<br>'
    for line in stdout.getvalue().split('\n'):
        print line,'<br>'



# This script has been called by a html form. Get the data from there:
form = cgi.FieldStorage()
# We use the current time as name for the directory where we store 
# things. This makes it pretty unique, hard to guess by hand (writing 
# a script is easy) and easy to remove old data.
import time
timestamp = time.time()
odir = 'out/%s/'%timestamp
os.mkdir(odir)
os.chdir(odir)

# We would like to print the log to a file and the webpage.
# For the latter we need specific make-up, which we'll do ourselfs.
stderr_old = sys.stderr     # To later restore stderr
stdout_old = sys.stdout     # To later restore stdout
sys.stdout = stdout = StringIO.StringIO()
sys.stderr = stderr = StringIO.StringIO()

# Based upn the form data we have to generate the arguments that would
# be given on the command line.
args = []
# On the form we can have either a uploaded file, or...
if form.has_key('inputfile') and form['inputfile'].value.strip() != '':
    base,ext = form['inputfile'].filename.rsplit('.',1)
    open(base+'.'+ext,'w').writelines(form['inputfile'].value)
    args+= ['-f',base+'.'+ext,'-x',base+'_CG.'+ext,'-o',base+'_CG.top']
# A pdb id we have to get from the database
elif form.has_key('pdbid') and len(form['pdbid'].value.strip())==4:
    # For the portability of the web code to newer martinize versions, it is better to have 
    # this here in stead of when reading in the pdb file
    try:
        url = 'http://www.rcsb.org/pdb/files/%s.pdb.gz'%form['pdbid'].value
        response = urllib2.urlopen(url)
        data = StringIO.StringIO(response.read())
        s = gzip.GzipFile(fileobj=data)
        open('%s.pdb'%(form['pdbid'].value),'w').writelines(s)
        args+= ['-f','%s.pdb'%(form['pdbid'].value),'-o','%s_CG.top'%(form['pdbid'].value),'-x','%s_CG.pdb'%(form['pdbid'].value)]
    except:
        print "Can't find %s structure in protein database.\n"%form['pdbid'].value
        restore_std(stdout,stdout_old,stderr,stderr_old)
        sys.exit()
else:
    print "You didn't specify a valid input file or PDB-code"
    restore_std(stdout,stdout_old,stderr,stderr_old)
    sys.exit()

# If a secondary structure string is given, we use that, else...
if form.has_key('ss') and form['ss'].value.strip()!='':
    args += ['-ss',form['ss'].value]
# use the a local installation of dssp to determine the secondary structure.
else:
    args += ['-dssp',dssp]

# Check if it is a valid forcefield.
if form['forcefield'].value in martinize.forcefields:
    args += '-ff %s'%form['forcefield'].value
else:
    print 'The %s forcefield has not yet been inplemented.'%form['forcefield'].value

# Cysteine bridges are defined slightly different in the webform then online.
if form.has_key('cys') and form['cys'].value!='':
    for ssbridge in form['cys'].value.split('-'):
        args += ['-cys',ssbridge]

if form.has_key('collagen'): args.append('-collagen')
if form.has_key('neutral_term'): args.append('-nt')
if form.has_key('charged_breaks'): args.append('-cb')

# Here the actual stuff is happening!
options,lists = DOC.options,DOC.lists
open('tmp.tmp','w').write(args)
options = CMD.option_parser(args,options,lists,version)
MAIN.main(options)
 
restore_std(stdout,stdout_old,stderr,stderr_old)
print " <br> The files can be downloaded <a href= 'http://md.chem.rug.nl/~djurre/martinize/%s'>here.</a>"%odir
print "<br> For example using the command: 'wget -r -np -nH --cut-dirs=3 -R index.htm* http://md.chem.rug.nl/~djurre/martinize/%s'<br><br>"%odir


# Clean up the 'out' directory: everything older then a week can be removed.
os.chdir('../..')
dirlist = os.listdir('out/')
dirlist.remove('itps')
dirlist.remove('index.htm')
for dir in dirlist:
    if time.time() - float(dir) >= 604800:    # That's a week in seconds...
        [os.remove('out/%s/%s'%(dir,f)) for f in os.listdir('out/'+dir)]   # Python can't recursively delete directories
        os.rmdir('out/'+dir)
