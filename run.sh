conda create --yes -n cheminf -c https://conda.binstar.org/rmcgibbo rmcgibbo-cheminf
source activate cheminf

while true; do
    fab generate_qchemfile
    fab run_qchem
    fab upload
done
