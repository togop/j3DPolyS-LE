#! /bin/bash

PARAMS=""
# default values
OUT_DIR=./out
ANALYSE=""
BOUNDARY=""
BOUNDARY_FACTOR=""
BOUNDARY_DIRECTION=""
BOUNDARY_SCORE=""
Z_LOOP=""
UNIDIRECTIONAL=""
RADIUS_CONTACT=""
CONTACT_PROBABILITY=""
NLEF=""
KM=""
INIT_MODE=""
LOG=""
HIC3D=""
PCA=""
INPUT_DAT_FILE=./input.dat

while (( "$#" )); do
  case $1 in
    -o:*)
        OUT_DIR=${1:3}
        echo "argument -o: $OUT_DIR"
        shift
        ;;
    --output_folder:*)
        OUT_DIR=${1:16}
        echo "argument -o: $OUT_DIR"
        shift
        ;;
    -a:*|--analyse:*)
        ANALYSE=$1
        echo "argument ANALYSE: $1"
        shift
        ;;
    -b:*|--boundary:*)
        BOUNDARY=$1
        echo "argument BOUNDARY: $1"
        shift
        ;;
    -bf:*|--boundary_factor:*)
        BOUNDARY_FACTOR=$1
        echo "argument BOUNDARY_FACTOR: $1"
        shift
        ;;
    -bd:*|--boundary_direction:*)
        BOUNDARY_DIRECTION=$1
        echo "argument BOUNDARY_DIRECTION: $1"
        shift
        ;;
    -bs|--boundary_score)
        BOUNDARY_SCORE=$1
        echo "argument BOUNDARY_SCORE: $1"
        shift
        ;;
    -z|--z_loop)
        Z_LOOP=$1
        echo "argument Z_LOOP: $1"
        shift
        ;;
    -u|--unidirectional)
        UNIDIRECTIONAL=$1
        echo "argument UNIDIRECTIONAL: $1"
        shift
        ;;
    -r:*|--radius_contact:*)
        RADIUS_CONTACT=$1
        echo "argument RADIUS_CONTACT: $1"
        shift
        ;;
    -cp|--contact_probability)
        CONTACT_PROBABILITY=$1
        echo "argument CONTACT_PROBABILITY: $1"
        shift
        ;;
    -l:*|--nlef:*)
        NLEF=$1
        echo "argument NLEF: $1"
        shift
        ;;
    -m:*|--km:*)
        KM=$1
        echo "argument KM: $1"
        shift
        ;;
    -im:*|--init_mode:*)
        INIT_MODE=$1
        echo "argument init_mode: $1"
        shift
        ;;
    --hic3d:*)
        HIC3D=$1
        echo "argument hic3d: $1"
        shift
        ;;
    --log:*)
        LOG=$1
        echo "argument LOG: $1"
        shift
        ;;
    --pca)
        PCA=$1
        echo "argument PCA: $1"
        shift
        ;;
    --) # end argument parsing
        shift
        break
        ;;
    -*|--*=) # unsupported flags
        echo "Error in mpirun_3dpolys-le.sh: Unsupported flag $1" >&2
        exit 1
        ;;
    *) # preserve positional arguments
        #INPUT_DAT_FILE="$1"
        PARAMS="$PARAMS $1"
        shift
        ;;
  esac
done
echo "OUT_DIR: $OUT_DIR"

# set positional arguments in their proper place
eval set -- "$PARAMS"

if [[ -z "$1" ]]
then
    INPUT_DAT_FILE=./input.dat
else
    INPUT_DAT_FILE=$1
fi

#if [[ -z "$2" ]]
#then
#    OUT_DIR=./out
#else
#    OUT_DIR=$2
#fi

if [[ "$(ls -A ${OUT_DIR})" && -z "$ANALYSE" ]]
then
    echo "WARN: Clear not empty output folder: ${OUT_DIR}"
    rm -r ${OUT_DIR}
fi

echo "CALL: mpirun 3dpolys-le -o:${OUT_DIR} ${ANALYSE} ${BOUNDARY} ${BOUNDARY_FACTOR} ${BOUNDARY_DIRECTION} ${BOUNDARY_SCORE} ${Z_LOOP} ${UNIDIRECTIONAL} ${RADIUS_CONTACT} ${CONTACT_PROBABILITY} ${NLEF} ${KM} ${INIT_MODE} ${LOG} ${HIC3D} ${PCA} ${INPUT_DAT_FILE}" # --mca orte_base_help_aggregate 0
mpirun             3dpolys-le -o:${OUT_DIR} ${ANALYSE} ${BOUNDARY} ${BOUNDARY_FACTOR} ${BOUNDARY_DIRECTION} ${BOUNDARY_SCORE} ${Z_LOOP} ${UNIDIRECTIONAL} ${RADIUS_CONTACT} ${CONTACT_PROBABILITY} ${NLEF} ${KM} ${INIT_MODE} ${LOG} ${HIC3D} ${PCA} ${INPUT_DAT_FILE}
