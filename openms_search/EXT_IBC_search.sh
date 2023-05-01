
spec_dir="../spec/EXT_IBC/"

if ls "$spec_dir"/*.mzML;
then
	echo "Skip converting"
else
	./convert_raw_to_mzML.sh "$spec_dir"
fi

./run_search_dir.sh "$spec_dir/*.mzML" "../db/EXT_IBC/9606_human_nodecoy_060421_complete_iRT_RPS27A_split.fasta"


