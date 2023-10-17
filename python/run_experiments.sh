
experiments=$(cat <<- DOC_EXPERIMENTS
	unweighted_pdf_mode
DOC_EXPERIMENTS
)
# 	unweighted_pdf_weighted_pdf_mode
# 	unweighted_cdf_weighted_pdf_mode
# 	unweighted_pdf_cdf_weighted_pdf_mode

# skew=$(sed -ne 's/^\s*init_skewness = \([0-9]*\).*/\1/p' run.py)
# echo ${skew}x
for conf in $experiments;
do
	CUDA_VISIBLE_DEVICES=$1 python run.py $conf &> ../log/${conf}.log &
done

tail -f ../log/${conf}.log



