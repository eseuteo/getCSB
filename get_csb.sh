#! /bin/bash

head -n 16 $1 > ./header_file.txt
tail -n +17 $1 > ./frags_file.csv
python3 main.py ./header_file.txt ./frags_file.csv
rm header_file.txt
rm frags_file.csv
rm temp_csb.csv
