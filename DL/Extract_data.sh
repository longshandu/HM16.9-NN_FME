#!/usr/bin/env bash

# Quantization Parameters
QPs=(22 27 32 37)

# Video Title
VID="BlowingBubbles"

# Video title is saved in VID variable
# "echo |" is put before the command in order to simulate "ENTER" keypress after HM finishes running
# Rename the output files
# Rinse and Repeat

for qp in "${QPs[@]}";
do
    echo | ../bin/TAppEncoderStatic -c ../cfg/encoder_lowdelay_P_main.cfg -c ../cfg/per-sequence/$VID.cfg -q $qp
    mv SSE.csv SSE_$qp.csv
done

rm rec.yuv str.bin
