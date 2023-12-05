
# datasets=$(cat <<- DOC_DATASETS
#     alban
#     Alinden
#     ALott
#     CPSF
#     D1810
#     ecoli_xl
#     MS2000225
#     peplib
#     QE
#     RPA
# DOC_DATASETS
# )

conda activate cuda11.8

declare -a datasets
datasets=(alban Alinden ALott CPSF D1810 ecoli_xl MS2000225 peplib QE RPA)

# 	unweighted_pdf_weighted_pdf_mode
# 	unweighted_cdf_weighted_pdf_mode
# 	unweighted_pdf_cdf_weighted_pdf_mode

# skew=$(sed -ne 's/^\s*init_skewness = \([0-9]*\).*/\1/p' run.py)
# echo ${skew}x
# if [[ ! -z "$2" ]]; then
# 	conf=-c $2
# fi
i=$1
echo ${datasets[i]}

dataset=${datasets[i]}

conf=$2
part=$3
mu_strategy=$4

python run.py -c $conf -d $dataset -q $part --suffix _${mu_strategy} --mu_strategy ${mu_strategy} -r 40 -j 10 -p -i # &> ../log/${conf}_${dataset}.log

# tail -f ../log/${conf}_${dataset}.log



