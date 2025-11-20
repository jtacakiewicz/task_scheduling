#!/usr/bin/env bash

for idx in 155975; do
    echo "=============== Solver ${idx} ==============="
    for x in 50 100 150 200 250 300 350 400 450 500; do
        echo ">>>>>> Instance ${x}"
        ./${idx}* in/in_155963_${x}.txt  out/out.txt
        python valid.py in/in_155963_${x}.txt  out/out.txt
        result=${?}
        if [[ $result == 1 ]]; then
            echo "ERRORS DETECTED!"
        else
            val=$(head -n 1 out/out.txt)
            echo "U for ${idx};${x}:    $val"
        fi

        echo "<<<<<<"
    done
done

