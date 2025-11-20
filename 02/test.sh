#!/usr/bin/env bash



for idx in 154978 155826 155851 155866 155868 155891 155892 155905 155913 155954 155958 155963 155967 155975 156010 156029; do; do
    echo "=============== Solver ${idx} ==============="
    for x in 50 100 150 200 250 300 350 400 450 500; do
        echo ">>>>>> Instance ${x}"
        { time ./solvers/${idx}* in/in_155963_${x}.txt  out/out.txt $((x / 10)); } 2> /tmp/result${idx}_${x}.txt
        time=`head /tmp/result${idx}_${x}.txt -n 2 | tail -n 1 | sed 's/.*real.*.m//g'`
        echo -e "\e[34mTime taken: ${time}\e[0m"
        python valid.py in/in_155963_${x}.txt  out/out.txt
        result=${?}
        if [[ $result == 1 ]]; then
            echo -e "\e[31mERRORS DETECTED!\e[0m"
        else
            val=$(head -n 1 out/out.txt)
            echo -e "\e[32mU for ${idx};${x}:    $val\e[0m"
        fi

        echo "<<<<<<"
    done
done

