conda create --yes -n cheminf -c https://conda.binstar.org/rmcgibbo rmcgibbo-cheminf
source activate cheminf

THISDIR=`pwd`

while true; do
    TEMPDIR=`mktemp -d`
    cd $TEMPDIR
    fab --fabfile=$THISDIR/fabfile.py generate_qchemfile
    fab --fabfile=$THISDIR/fabfile.py run_qchem
    fab --fabfile=$THISDIR/fabfile.py upload
    cd $THISDIR
    rm -rf $TEMPDIR
done
