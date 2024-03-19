# Pre-processing of all kinetics schemes
source ~/.zshrc
set-opensmoke-0.21

BASEDIR=${PWD}

run_preprocessor_in_folders() {

  echo "Preprocessing $1"
  for folder in $1/*; do
    cd "$folder"
    OpenSMOKEpp_CHEMKIN_PreProcessor.sh
    cd $BASEDIR
  done

}

run_preprocessor_in_folders "one-step"
run_preprocessor_in_folders "two-step"
run_preprocessor_in_folders "materials"
run_preprocessor_in_folders "skeletal"
run_preprocessor_in_folders "evaporation"
run_preprocessor_in_folders "biomass"

