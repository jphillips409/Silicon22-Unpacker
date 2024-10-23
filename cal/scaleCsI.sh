#!/bin/bash

# Script for calibrating CsI using proton pedestal and punchthrough data
#
# How to run:
# ./scaleCsI.sh [pedestal/punchthrough input.txt] [calibration output.txt]

#
# Format expected in pedestal/punchthrough input file:
# ----------------------------------------------------------------------------
# Det(#) Strip(#) Ped(channels) Ped(MeV) Punchthr(channels) Punchthr(MeV)
# ----------------------------------------------------------------------------
#
# Format produced in pedestal/punchthrough input file:
# ----------------------------------------------------------------------------
# Det(#) Strip(#) Slope Intercept
# ----------------------------------------------------------------------------

if [ -s $2 ]
then
    rm $2
fi

while read l
do
    e=($l)
    echo -n "${e[0]} ${e[1]} " >> $2

    # calculate linear calibration as: y = mx + b
    m=$(bc <<< "scale=4;${e[5]}/(${e[4]}-${e[2]})")
    b=$(bc <<< "scale=4;${e[3]}-$m*${e[2]}")

    echo "$m $b" >> $2
done < $1
