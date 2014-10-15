# conda create --yes -n cheminf -c https://conda.binstar.org/rmcgibbo rmcgibbo-cheminf
set -e
source activate cheminf

THISDIR=`python2 -c "import os; print os.path.abspath(os.path.dirname('$0'))"`

while true; do
    TEMPDIR=`mktemp -d`
    cd $TEMPDIR
    fab --fabfile=$THISDIR/fabfile.py generate_qchemfile
    fab --fabfile=$THISDIR/fabfile.py run_qchem
    fab --fabfile=$THISDIR/fabfile.py upload
    cd $THISDIR
    rm -rf $TEMPDIR
done
