
datasets=$(cat <<- DOC_DATASETS
    alban
    Alinden
    ALott
    CPSF
    D1810
    ecoli_xl
    MS2000225
    peplib
    QE
    RPA
DOC_DATASETS
)
# 	unweighted_pdf_weighted_pdf_mode
# 	unweighted_cdf_weighted_pdf_mode
# 	unweighted_pdf_cdf_weighted_pdf_mode

# skew=$(sed -ne 's/^\s*init_skewness = \([0-9]*\).*/\1/p' run.py)
# echo ${skew}x
# if [[ ! -z "$2" ]]; then
# 	conf=-c $2
# fi
conf=$2
for dataset in $datasets;
do
	CUDA_VISIBLE_DEVICES=$1 python run.py -c $conf -d $dataset &> ../log/${conf}_1S_${dataset}.log &
	echo $! >> running_pids.txt
done

tail -f ../log/${conf}_1S_*.log



