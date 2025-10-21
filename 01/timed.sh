#!/usr/bin/env bash
# Usage: ./timed.sh <timeout_in_seconds> <command> [args...]
# Example: ./timed.sh 5 sleep 10
#          ./timed.sh 5 ls -la

TIMEOUT=$1
TIMEOUT_LONG=$(($1 * 2))
shift

# Use a temp file to capture the time output
TIMEFILE=$(mktemp)

# Run the command with /usr/bin/time, capturing its output
# & put it in the background to monitor its PID
/usr/bin/env timeout $TIMEOUT_LONG /usr/bin/env time -f "%e" -o "$TIMEFILE" "$@" > /dev/stderr 2>&1
EXIT_CODE=$?

# Read the elapsed time from the temp file
if [[ -f "$TIMEFILE" ]]; then
    ELAPSED_TIME=$(cat "$TIMEFILE")
    rm -f "$TIMEFILE"
fi

# Show results
if [[ "$exit_code" -eq 0 && $(python3 -c "print(1 if float('$ELAPSED_TIME') <= float('$TIMEOUT') else 0)") -eq 1 ]]; then
    echo "Completed successfully in: ${ELAPSED_TIME}s."
else
    if [ -z $ELAPSED_TIME ]; then
        echo "Command was killed after: ${TIMEOUT_LONG}s."
    else
        echo "Command took too long: ${ELAPSED_TIME}s."
    fi
fi

exit $EXIT_CODE
