
raw_dir="$1"
raw_dir=$(realpath "$raw_dir")

pushd /home/yisupeng/softs/ThermoRawFileParser/

mono ThermoRawFileParser.exe -d $raw_dir

popd #/home/yisupeng/softs/ThermoRawFileParser/

