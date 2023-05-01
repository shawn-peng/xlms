
spec_dir="../spec/KKT_CPC/"

if ls "$spec_dir"/*.mzML;
then
	echo "Skip converting"
else
	./convert_raw_to_mzML.sh "$spec_dir"
fi

./run_decoy_search_dir.sh "$spec_dir/*.mzML" "../db/KKT_CPC/Tbrucei927_v4_2_with_small_proteomes_Akiyoshi_20190124_decoy.fasta" BS3



