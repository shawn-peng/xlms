#OpenPepXLLF -in /home/yisupeng/workspace/ms/spec/xl/ecoli_xl/DESTRibo_092010exp2_LEF10.mzML -database /home/yisupeng/workspace/xlms/XLSearch/testDatabase/ecoli_riboprot_fwd.fasta -algorithm:number_top_hits 10 -out_mzIdentML ../results/openpepxllf/knime4.6/DESTRibo_092010exp2_LEF10.mzid

SPEC_PATH_PATTERN="$1"
echo "$SPEC_PATH_PATTERN"
DB="$2"
CROSSLINKER="$3"
if [ -z "$3" ];
then
	CROSSLINKER="DSS"
fi
echo CROSSLINKER $CROSSLINKER
SPECDIR=$(dirname "$SPEC_PATH_PATTERN")
echo SPECDIR $SPECDIR
RESDIR=../results/openpepxllf/knime4.6/$(basename $SPECDIR)_decoy
mkdir -p $RESDIR
echo RESDIR $RESDIR
EXT=$(echo "$SPEC_PATH_PATTERN" | sed -e "s/.*\.\(\w*\)$/\1/")
echo EXT $EXT
LOGDIR=log
echo Logging into $LOGDIR/$(basename $SPECDIR)_decoy
mkdir -p $LOGDIR/$(basename $SPECDIR)_decoy

for f in $SPEC_PATH_PATTERN;
do
	echo $(basename $f)
	echo OpenPepXLLF -in $f -database $DB -threads 20 -algorithm:number_top_hits 10 -out_idXML $RESDIR/$(echo $(basename $f) | sed -e "s/$EXT/idXML/") &> $LOGDIR/$(basename $SPECDIR)_decoy/$(echo $(basename $f) | sed -e "s/$EXT/log/")
	OpenPepXLLF -in $f -database $DB -threads 20 -algorithm:number_top_hits 10 -decoy_string reverse_ -out_idXML $RESDIR/$(echo $(basename $f) | sed -e "s/$EXT/idXML/") &>> $LOGDIR/$(basename $SPECDIR)_decoy/$(echo $(basename $f) | sed -e "s/$EXT/log/") &
done

