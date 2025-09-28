#!/bin/bash

cd src/
make clean
make
cd ..

dir="./instances/"
find "$dir" -type f -iname "*.txt" | while read file; do
    echo "$file"

    count=0
    for i in {1..1}; do
        echo "Starting thread with seed $i for file $file" &

        src/is -f "$file" -s "$i" &

        ((count++))
        if (( count % 3 == 0 )); then
            wait
        fi
    done

    wait 
done
