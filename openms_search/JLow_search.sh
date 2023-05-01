
spec_dir="../spec/JLow/"

if ls "$spec_dir"/*.mzML;
then
	echo "Skip converting"
else
	./convert_raw_to_mzML.sh "$spec_dir"
fi

./run_search_dir.sh "$spec_dir/*.mzML" "../db/JLow/UP000005640_9606.fasta" BS3


