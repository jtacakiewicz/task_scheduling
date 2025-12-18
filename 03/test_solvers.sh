#!/usr/bin/env bash

CHECK=weryfikator
solvers="154978 155826 155851 155866 155868 155891 155892 155905 155913 155954 155958 155963 155967 155975 156010 156029"
instances="$solvers"
instances="155892"

solvers="155975"
RESULTS="result.txt"

echo "=============== Weryfikator ${idx} ==============="
for instance in ${instances}; do
    for solver in ${solvers}; do
        echo ">>>>>>>>>>>>> Solver ${solver}"
        echo ">> $solver" >> "$RESULTS"
        for x in 50 100 150 200 250 300 350 400 450 500; do
            echo ">>>>>> Instance ${x} ${instance}"
            in_file=in/in_${instance}_${x}.txt
            out_file=out/out_${x}.txt
            time ./solvers/${solver} "$in_file" "$out_file" $((x / 10));
            ./${CHECK} "$in_file" "$out_file"

            result=$?
            if [[ $result != 0 ]]; then
                echo -e "\e[31mERRORS DETECTED! code: $result\e[0m"
                echo "-1" >> "$RESULTS"
            else
                val=$(head -n 1 $out_file)
                echo "$val" >> "$RESULTS"
                echo -e "\e[32mU for ${idx};${x}:    $val\e[0m"
            fi

            echo "<<<<<<"
        done
    done
done
