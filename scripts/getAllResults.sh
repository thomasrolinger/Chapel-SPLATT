#!/bin/sh

# YELP
> YELP.csv

# NELL-2
> NELL-2.csv

# NELL-1
> NELL-1.csv

echo""
echo "##############################################"
echo "+ Processing output for C......"
echo""

python parseTimings.py ../benchmark/output_data_C/YELP/ output_C.csv C
cat output_C.csv > YELP.csv

python parseTimings.py ../benchmark/output_data_C/NELL-2/ output_C.csv C
cat output_C.csv > NELL-2.csv

python parseTimings.py ../benchmark/output_data_C/NELL-1/ output_C.csv C
cat output_C.csv > NELL-1.csv

echo""
echo "##############################################"
echo "+ Processing output for Chapel......"
echo""

python parseTimings.py ../benchmark/output_data_CHAPEL/YELP/ output_CHAPEL.csv CHAPEL
tail --lines=+2 output_CHAPEL.csv >> YELP.csv

python parseTimings.py ../benchmark/output_data_CHAPEL/NELL-2/ output_CHAPEL.csv CHAPEL
tail --lines=+2 output_CHAPEL.csv >> NELL-2.csv

python parseTimings.py ../benchmark/output_data_CHAPEL/NELL-1/ output_CHAPEL.csv CHAPEL
tail --lines=+2 output_CHAPEL.csv >> NELL-1.csv

echo""
echo "##############################################"
echo "+ DONE"
echo""

