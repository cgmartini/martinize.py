# Testset for martinize script. May 2012, DdJ
# Updated: October 2015
# Run using: ./test.sh 2>&1 >&log&
# There are two test: 
#   none of the calls should give errors/break
#   The output is getting diff'ed with an older run to check for consistancy.
# The latter ofcourse breaks if you do not have any previous runs...
# It will give small difference for the pdb-files of polarizable FFs: particles use random placement.
# Maybe we need to add 'failing' tests?

SCRIPT='../../../GIT/martinize-2.6a.py'
DIFF='V2.5/'

function pdbgett {
    wget http://www.rcsb.org/pdb/files/$1.pdb.gz
    gzip -d $1.pdb.gz
}

VERSION=V`echo ${SCRIPT:0:${#SCRIPT}-3} |cut -d '-' -f 2`
mkdir ${VERSION}
cd ${VERSION}

# Most standard function, of course with ubiquitin
block=standard21
echo '***************'${block}'*****************'
mkdir $block 
cd $block
pdbgett 1UBQ
echo 12|pdb2gmx -f 1UBQ.pdb -o 1UBQ.gro -water spc
$SCRIPT -f 1UBQ.gro -o 1UBQ_cg.top -x 1UBQ_cg.pdb -dssp $DSSP -ff martini21
diff *_cg.pdb ../../${DIFF}/${block}/*_cg.pdb >> ../${block}.diff
diff *_cg.top ../../${DIFF}/${block}/*_cg.top >> ../${block}.diff
cd ..

# Most standard function for M2.2, of course with ubiquitin
block=standard22
echo '***************'${block}'*****************'
mkdir $block 
cd $block
pdbgett 1UBQ
$SCRIPT -f 1UBQ.pdb -o 1UBQ_cg.top -x 1UBQ_cg.pdb -dssp $DSSP -ff martini22
diff *_cg.pdb ../../${DIFF}/${block}/*_cg.pdb >> ../${block}.diff
diff *_cg.top ../../${DIFF}/${block}/*_cg.top >> ../${block}.diff
cd ..

# Most standard functions for polarizable martini, of course with ubiquitin
block=standard22p
echo '***************'${block}'*****************'
mkdir $block 
cd $block
pdbgett 1UBQ
$SCRIPT -f 1UBQ.pdb -o 1UBQ_cg.top -x 1UBQ_cg.pdb -dssp $DSSP -ff martini22p
diff *_cg.pdb ../../${DIFF}/${block}/*_cg.pdb >> ../${block}.diff
diff *_cg.top ../../${DIFF}/${block}/*_cg.top >> ../${block}.diff
cd ..

# Elnedyn (1) with a protein that does contain TRP's 
block=elnedyn_TRP
echo '***************'${block}'*****************'
mkdir $block 
cd $block
pdbgett 1L35
$SCRIPT -f 1L35.pdb -o 1L35_cg.top -x 1L35_cg.pdb -dssp $DSSP -ff elnedyn
diff *_cg.pdb ../../${DIFF}/${block}/*_cg.pdb >> ../${block}.diff
diff *_cg.top ../../${DIFF}/${block}/*_cg.top >> ../${block}.diff
cd ..

# Elnedyn (1) and giving secondary sturcture as string. Use custom elastic force.
block=elnedyn1
echo '***************'${block}'*****************'
mkdir $block 
cd $block
pdbgett 1UBQ
$SCRIPT -f 1UBQ.pdb -o 1UBQ_cg.top -x 1UBQ_cg.pdb -ss ~EEEEEETTS~EEEEE~~TTSBHHHHHHHHHHHH~~~GGGEEEEETTEE~~TTSBTGGGT~~TT~EEEEEE~~S~~ -ff elnedyn -ef 700
diff *_cg.pdb ../../${DIFF}/${block}/*_cg.pdb >> ../${block}.diff
diff *_cg.top ../../${DIFF}/${block}/*_cg.top >> ../${block}.diff
cd ..

# Elnedyn (2) and giving secondary sturcture as string. Use custom elastic force.
block=elnedyn2
echo '***************'${block}'*****************'
mkdir $block 
cd $block
pdbgett 1UBQ
$SCRIPT -f 1UBQ.pdb -o 1UBQ_cg.top -x 1UBQ_cg.pdb -ss ~EEEEEETTS~EEEEE~~TTSBHHHHHHHHHHHH~~~GGGEEEEETTEE~~TTSBTGGGT~~TT~EEEEEE~~S~~ -ff elnedyn22p
diff *_cg.pdb ../../${DIFF}/${block}/*_cg.pdb >> ../${block}.diff
diff *_cg.top ../../${DIFF}/${block}/*_cg.top >> ../${block}.diff
cd ..

# Test the seperate writing of chains, position restraints, neutral termini. Use Mscl.
block=chain21
echo '***************'${block}'*****************'
mkdir $block 
cd $block
pdbgett 2OAR
$SCRIPT -f 2OAR.pdb -o 2OAR_cg.top -x 2OAR_cg.pdb -sep -nt -p All -pf 500 -dssp $DSSP -ff martini21
diff *_cg.pdb ../../${DIFF}/${block}/*_cg.pdb >> ../${block}.diff
diff *_cg.top ../../${DIFF}/${block}/*_cg.top >> ../${block}.diff
cd ..

# Test the seperate writing of chains, position restraints, neutral termini. Use Mscl.
block=chains22
echo '***************'${block}'*****************'
mkdir $block 
cd $block
pdbgett 2OAR
$SCRIPT -f 2OAR.pdb -o 2OAR_cg.top -x 2OAR_cg.pdb -sep -nt -p All -pf 500 -dssp $DSSP -ff martini22
diff *_cg.pdb ../../${DIFF}/${block}/*_cg.pdb >> ../${block}.diff
diff *_cg.top ../../${DIFF}/${block}/*_cg.top >> ../${block}.diff
cd ..

# Test cysteine bridge (between cys 9 and cys 164), naming protein and use dihedrals for extended regions. Use mutated lysozyme.
block=Sbonds21
echo '***************'${block}'*****************'
mkdir $block 
cd $block
pdbgett 1L35 
$SCRIPT -f 1L35.pdb -o 1L35_cg.top -x 1l35_cg.pdb -cys auto -name lysozyme -dssp $DSSP -ed -ff martini21
diff *_cg.pdb ../../${DIFF}/${block}/*_cg.pdb >> ../${block}.diff
diff *_cg.top ../../${DIFF}/${block}/*_cg.top >> ../${block}.diff
cd ..

# Test cysteine bridge (between cys 9 and cys 164), naming protein and use dihedrals for extended regions. Use mutated lysozyme.
block=Sbonds22
echo '***************'${block}'*****************'
mkdir $block 
cd $block
pdbgett 1L35 
$SCRIPT -f 1L35.pdb -o 1L35_cg.top -x 1l35_cg.pdb -cys auto -name lysozyme -dssp $DSSP -ed -ff martini22
diff *_cg.pdb ../../${DIFF}/${block}/*_cg.pdb >> ../${block}.diff
diff *_cg.top ../../${DIFF}/${block}/*_cg.top >> ../${block}.diff
cd ..

# Test collagen paramters on of the six collagen parameters.
block=collagen21
echo '***************'${block}'*****************'
mkdir $block 
cd $block
pdbgett 1CAG
$SCRIPT -f 1CAG.pdb -o 1CAG_cg.top -x 1CAG_cg.pdb -collagen -ff martini21
diff *_cg.pdb ../../${DIFF}/${block}/*_cg.pdb >> ../${block}.diff
diff *_cg.top ../../${DIFF}/${block}/*_cg.top >> ../${block}.diff
cd ..

# Test collagen paramters on of the six collagen parameters.
block=collagen22
echo '***************'${block}'*****************'
mkdir $block 
cd $block
pdbgett 1CAG
$SCRIPT -f 1CAG.pdb -o 1CAG_cg.top -x 1CAG_cg.pdb -collagen -ff martini22
diff *_cg.pdb ../../${DIFF}/${block}/*_cg.pdb >> ../${block}.diff
diff *_cg.top ../../${DIFF}/${block}/*_cg.top >> ../${block}.diff
cd ..

# Test mixed DNA-protein.
block=DNAprotein
echo '***************'${block}'*****************'
mkdir $block 
cd $block
pdbgett 3SJM
$SCRIPT -f 3SJM.pdb -o 3SJM_cg.top -x 3SJM_cg.pdb -collagen -ff martini22dna
diff *_cg.pdb ../../${DIFF}/${block}/*_cg.pdb >> ../${block}.diff
diff *_cg.top ../../${DIFF}/${block}/*_cg.top >> ../${block}.diff
cd ..

cd ..
