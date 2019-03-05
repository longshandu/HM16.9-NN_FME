#!/usr/bin/env bash

# HM-Tool Directory
HM=~/git-repos/HM16.9

# Quantization Parameters
QPs=(22 27 32 37)

# Video Title
VID="BlowingBubbles"

# Video title is saved in VID variable
# "echo |" is put before the command in order to simulate "ENTER" keypress after HM finishes running
# Rename the output files
# Rinse and Repeat

cd $HM/bin
for qp in "${QPs[@]}";
do
    echo | ./TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp
    mv SSE.csv SSE_$qp.csv
done
