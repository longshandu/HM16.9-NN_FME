#!/usr/bin/env bash

# HM-Tool Directory
HM=~/git-repos/HM16.9

# Useful Parameters
fs=0
f=600

# Video title is saved in VID variable
# "echo |" is put before the command in order to simulate "ENTER" keypress after HM finishes running
# Rename the output files
# Rinse and Repeat

qp=22
VID="BlowingBubbles"
echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp

mv SSE.csv SSE_$qp.csv

qp=27
echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp

mv SSE.csv SSE_$qp.csv

qp=32
echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp

mv SSE.csv SSE_$qp.csv

qp=37
echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qps

mv SSE.csv SSE_$qp.csv


################################

#qp=27
#VID="Kimono"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 5 -fs $fs

#VID="BQTerrace"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 5 -fs $fs

#VID="PartyScene"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 20 -fs $fs

#VID="BlowingBubbles"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 30 -fs $fs

#VID="BQSquare"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 30 -fs $fs

#VID="RaceHorses"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 30 -fs $fs

#VID="Johnny"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 10 -fs $fs

#VID="SlideShow"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 10 -fs $fs

#VID="SlideEditing"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 10 -fs $fs

#mv mv_nn.csv ../positions/mv_$qp.csv
#mv SSE_errors.csv ../positions/SSE_$qp.csv

##################################

#qp=32
#VID="Kimono"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 5 -fs $fs

#VID="BQTerrace"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 5 -fs $fs

#VID="PartyScene"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 20 -fs $fs

#VID="BlowingBubbles"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 30 -fs $fs

#VID="BQSquare"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 30 -fs $fs

#VID="RaceHorses"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 30 -fs $fs

#VID="Johnny"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 10 -fs $fs

#VID="SlideShow"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 10 -fs $fs

#VID="SlideEditing"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 10 -fs $fs

#mv mv_nn.csv ../positions/mv_$qp.csv
#mv SSE_errors.csv ../positions/SSE_$qp.csv

##################################

#qp=37
#VID="Kimono"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 5 -fs $fs

#VID="BQTerrace"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 5 -fs $fs

#VID="PartyScene"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 20 -fs $fs

#VID="BlowingBubbles"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 30 -fs $fs

#VID="BQSquare"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 30 -fs $fs

#VID="RaceHorses"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 30 -fs $fs

#VID="Johnny"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 10 -fs $fs

#VID="SlideShow"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 10 -fs $fs

#VID="SlideEditing"
#echo | $HM/bin/TAppEncoderStatic -c $HM/cfg/encoder_lowdelay_P_main.cfg -c $HM/cfg/per-sequence/$VID.cfg -q $qp -f 10 -fs $fs

#mv mv_nn.csv ../positions/mv_$qp.csv
#mv SSE_errors.csv ../positions/SSE_$qp.csv

#################################

