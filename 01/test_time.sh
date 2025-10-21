#!/usr/bin/env bash

if [ $# -lt 1 ]; then
    echo "Usage: ./test_time.sh <program> [INDEX] [LOGS] [MULT]"
    exit 1
fi

EXEC="$1"
IDXS="${2:-"154978,155868,155913,155967,155826,155891,155954,155975,155851,155892,155958,156010,155866,155905,155963,156029"}"
USE_LOGS="${3:-}"
MULT="${4:-10}"

OUTFILE=$(mktemp)
for idx in $(echo "$IDXS" | tr ',' '\n'); do
    size=50
    while (( size <= 500 )); do
        LOG="logs/in_${idx}_${size}.log"
        if [ ! -z $USE_LOGS ]; then
            LOG=$USE_LOGS
        fi
        DATA="instances/in_${idx}_${size}.txt"
        touch $DATA
        echo ">>>>>>>>> checking $DATA \w $EXEC ..."

        time=$((size/MULT))

        echo "========= solver:" 1>> $LOG
        ./timed.sh ${time} ${EXEC} ${DATA} ${OUTFILE} ${time} 2> $LOG
        echo "========= checker:" 1>> $LOG
        CMAX=`./checker.py $DATA ${OUTFILE} 2>> $LOG`
        EXIT_CODE=$?
        if [ $EXIT_CODE -eq 0 ]; then
            echo "CMax: $CMAX"
        else
            echo "Incorrect solution! (check logs)"
        fi
        echo "<<<<<<<<<"
        size=$((size + 50))
    done
done

rm -f $OUTFILE
