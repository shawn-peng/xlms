
spec_dir="../spec/RPA/"

if ls "$spec_dir"/*.mzML;
then
	echo "Skip converting"
else
	./convert_raw_to_mzML.sh "$spec_dir"
fi

./run_decoy_search_dir.sh "$spec_dir/*.mzML" "../db/RPA/subjects_decoy.fasta" BS3


