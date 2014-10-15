while true; do
    fab generate_qchemfile
    fab run_qchem
    fab upload
done
