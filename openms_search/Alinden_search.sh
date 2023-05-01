
spec_dir="../spec/Alinden/"

if ls "$spec_dir"/*.mzML;
then
	echo "Skip converting"
else
	./convert_raw_to_mzML.sh "$spec_dir"
fi

./run_search_dir.sh "$spec_dir/*.mzML" "../db/Alinden/eIF1-1A-2-3-4A-4B-40S-4G-DAP-4A2-4E_withIsoforms+Niels.fasta" BS3


