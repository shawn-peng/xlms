
spec_dir="../spec/alban/"

if ls "$spec_dir"/*.mzML;
then
	echo "Skip converting"
else
	./convert_raw_to_mzML.sh "$spec_dir"
fi

./run_search_dir.sh "$spec_dir/*.mzML" "../db/alban/Fibrin-FT-Top20.fasta"


