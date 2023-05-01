
spec_dir="../spec/peplib/"

if ls "$spec_dir"/*.mzML;
then
	echo "Skip converting"
else
	./convert_raw_to_mzML.sh "$spec_dir"
fi

./run_decoy_search_dir.sh "$spec_dir/*.mzML" "../db/peplib/Cas9_plus10_decoy.fasta"


