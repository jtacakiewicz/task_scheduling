#!/usr/bin/env bash

# Konfiguracja
VERIFIER_EXEC="weryfikator"
RESULTS_FILE="wyniki.txt"
TIMES_FILE="czasy.txt"

# Instancje i rozmiary
solvers="154978 155826 155851 155866 155868 155891 155892 155905 155913 155954 155958 155963 155967 155975 156010 156029"
sizes="50 100 150 200 250 300 350 400 450 500"
instance="155975"

# Czyszczenie plików wyjściowych przed startem
> "$RESULTS_FILE"
> "$TIMES_FILE"

for solver in ${solvers}; do
    echo "Processing Solver: ${solver}"
    
    # Określenie komendy uruchamiającej
    # Jeśli plik kończy się na .py, użyj python3, w przeciwnym razie uruchom bezpośrednio
    if [[ "$solver" == *.py ]]; then
        cmd="python3 $solver"
    else
        cmd="./$solver"
    fi

    for x in ${sizes}; do
        in_file="in/in_${instance}_${x}.txt"
        tmp_out="/tmp/out_${instance}_${x}.txt"
        time_limit=$((x / 10))
        
        start_time=$(date +%s%N)
        
        # Uruchomienie solvera (dynamiczna komenda $cmd)
        $cmd "$in_file" "$tmp_out" "$time_limit" > /dev/null 2>&1
        solver_res=$?
        
        end_time=$(date +%s%N)
        
        duration_ns=$((end_time - start_time))
        duration_s=$(awk "BEGIN {printf \"%.3f\", $duration_ns/1000000000}")

        if [[ $solver_res != 0 ]]; then
            echo "ERR" >> "$RESULTS_FILE"
            echo "$duration_s" >> "$TIMES_FILE"
            echo -e "Size ${x}: \e[31mSOLVER ERROR (code: $solver_res)\e[0m"
        else
            ./${VERIFIER_EXEC} "$in_file" "$tmp_out" > /dev/null 2>&1
            verifier_res=$?

            if [[ $verifier_res != 0 ]]; then
                echo "ERR" >> "$RESULTS_FILE"
                echo -e "Size ${x}: \e[31mVERIFIER ERROR\e[0m"
            else
                val=$(head -n 1 "$tmp_out")
                echo "$val" >> "$RESULTS_FILE"
                echo -e "Size ${x}: \e[32mOK ($val)\e[0m"
            fi
            echo "$duration_s" >> "$TIMES_FILE"
        fi
    done
done

echo "-------------------------------------------------------"
echo "Gotowe! Wyniki zapisano w $RESULTS_FILE, a czasy w $TIMES_FILE"
