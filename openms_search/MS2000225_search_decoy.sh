
spec_dir="../spec/MS2000225/"

if ls "$spec_dir"/*.mzML ;
then
	echo "Skip converting"
else
	./convert_raw_to_mzML.sh "$spec_dir"
fi

#./run_search_dir.sh "$spec_dir/*.mzML" "../db/MS2000225/uniprot_Human_20190401.fasta"
./run_decoy_search_dir.sh "$spec_dir/*.mzML" "../db/MS2000225/subjects_decoy.fasta"

