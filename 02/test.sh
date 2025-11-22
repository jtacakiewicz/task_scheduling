#!/usr/bin/env bash

instances=155975
solvers="154978 155826 155851 155866 155868 155891 155892 155905 155913 155954 155958 155963 155967 155975 156010 156029"
solvers="155975"
for idx in $solvers; do
    echo "=============== Solver ${idx} ==============="
    for x in 50 100 150 200 250 300 350 400 450 500; do
        echo ">>>>>> Instance ${x}"
        in_file=in/in_${instances}_${x}.txt
        out_file=/tmp/out.txt
        if [[ -f "$out_file" ]]; then
            rm "$out_file"
        fi

        { time ./${idx}* ${in_file} ${out_file} $((x / 10)); } 2> /tmp/result${idx}_${x}.txt
        result_solve=$?
        if [[ $result_solve == 0 ]]; then
            time=`head /tmp/result${idx}_${x}.txt -n 2 | tail -n 1 | sed 's/.*real.*.m//g'`
            echo -e "\e[34mTime taken: ${time}\e[0m"
            python valid.py ${in_file} ${out_file}
        fi

        result=${?} || $result_solve
        if [[ $result == 1 ]]; then
            echo -e "\e[31mERRORS DETECTED!\e[0m"
        else
            val=$(head -n 1 $out_file)
            echo -e "\e[32mU for ${idx};${x}:    $val\e[0m"
        fi

        echo "<<<<<<"
    done
done

