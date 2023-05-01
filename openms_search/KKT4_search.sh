
spec_dir="../spec/KKT4/"

if ls "$spec_dir"/*.mzML;
then
	echo "Skip converting"
else
	./convert_raw_to_mzML.sh "$spec_dir"
fi

./run_search_dir.sh "$spec_dir/*.mzML" "../db/KKT4/KKT4_BS3_search_database.fasta" BS3



