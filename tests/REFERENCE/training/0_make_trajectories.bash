#!/bin/bash
IS=0.1
IS2=0.2

m4 -DM4VAR_INITIALSTRAIN="[$IS,   0.0,   0.0,   0.0,   0.0,   0.0]" ProjectParameters.m4.json > ProjectParameters_00.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   $IS,   0.0,   0.0,   0.0,   0.0]" ProjectParameters.m4.json > ProjectParameters_01.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   0.0,   $IS,   0.0,   0.0,   0.0]" ProjectParameters.m4.json > ProjectParameters_02.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   0.0,   0.0,  $IS2,   0.0,   0.0]" ProjectParameters.m4.json > ProjectParameters_03.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   0.0,   0.0,   0.0,  $IS2,   0.0]" ProjectParameters.m4.json > ProjectParameters_04.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   0.0,   0.0,   0.0,   0.0,  $IS2]" ProjectParameters.m4.json > ProjectParameters_05.json
m4 -DM4VAR_INITIALSTRAIN="[$IS,   $IS,   0.0,   0.0,   0.0,   0.0]" ProjectParameters.m4.json > ProjectParameters_06.json
m4 -DM4VAR_INITIALSTRAIN="[$IS,   0.0,   $IS,   0.0,   0.0,   0.0]" ProjectParameters.m4.json > ProjectParameters_07.json
m4 -DM4VAR_INITIALSTRAIN="[$IS,   0.0,   0.0,  $IS2,   0.0,   0.0]" ProjectParameters.m4.json > ProjectParameters_08.json
m4 -DM4VAR_INITIALSTRAIN="[$IS,   0.0,   0.0,   0.0,  $IS2,   0.0]" ProjectParameters.m4.json > ProjectParameters_09.json
m4 -DM4VAR_INITIALSTRAIN="[$IS,   0.0,   0.0,   0.0,   0.0,  $IS2]" ProjectParameters.m4.json > ProjectParameters_10.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   $IS,   $IS,   0.0,   0.0,   0.0]" ProjectParameters.m4.json > ProjectParameters_11.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   $IS,   0.0,  $IS2,   0.0,   0.0]" ProjectParameters.m4.json > ProjectParameters_12.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   $IS,   0.0,   0.0,  $IS2,   0.0]" ProjectParameters.m4.json > ProjectParameters_13.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   $IS,   0.0,   0.0,   0.0,  $IS2]" ProjectParameters.m4.json > ProjectParameters_14.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   0.0,   $IS,  $IS2,   0.0,   0.0]" ProjectParameters.m4.json > ProjectParameters_15.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   0.0,   $IS,   0.0,  $IS2,   0.0]" ProjectParameters.m4.json > ProjectParameters_16.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   0.0,   $IS,   0.0,   0.0,  $IS2]" ProjectParameters.m4.json > ProjectParameters_17.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   0.0,   0.0,  $IS2,  $IS2,   0.0]" ProjectParameters.m4.json > ProjectParameters_18.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   0.0,   0.0,  $IS2,   0.0,  $IS2]" ProjectParameters.m4.json > ProjectParameters_19.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   0.0,   0.0,   0.0,  $IS2,  $IS2]" ProjectParameters.m4.json > ProjectParameters_20.json
m4 -DM4VAR_INITIALSTRAIN="[$IS,   $IS,   $IS,   0.0,   0.0,   0.0]" ProjectParameters.m4.json > ProjectParameters_21.json
m4 -DM4VAR_INITIALSTRAIN="[$IS,   $IS,   0.0,  $IS2,   0.0,   0.0]" ProjectParameters.m4.json > ProjectParameters_22.json
m4 -DM4VAR_INITIALSTRAIN="[$IS,   $IS,   0.0,   0.0,  $IS2,   0.0]" ProjectParameters.m4.json > ProjectParameters_23.json
m4 -DM4VAR_INITIALSTRAIN="[$IS,   $IS,   0.0,   0.0,   0.0,  $IS2]" ProjectParameters.m4.json > ProjectParameters_24.json
m4 -DM4VAR_INITIALSTRAIN="[$IS,   0.0,   $IS,  $IS2,   0.0,   0.0]" ProjectParameters.m4.json > ProjectParameters_25.json
m4 -DM4VAR_INITIALSTRAIN="[$IS,   0.0,   $IS,   0.0,  $IS2,   0.0]" ProjectParameters.m4.json > ProjectParameters_26.json
m4 -DM4VAR_INITIALSTRAIN="[$IS,   0.0,   $IS,   0.0,   0.0,  $IS2]" ProjectParameters.m4.json > ProjectParameters_27.json
m4 -DM4VAR_INITIALSTRAIN="[$IS,   0.0,   0.0,  $IS2,  $IS2,   0.0]" ProjectParameters.m4.json > ProjectParameters_28.json
m4 -DM4VAR_INITIALSTRAIN="[$IS,   0.0,   0.0,  $IS2,   0.0,  $IS2]" ProjectParameters.m4.json > ProjectParameters_29.json
m4 -DM4VAR_INITIALSTRAIN="[$IS,   0.0,   0.0,   0.0,  $IS2,  $IS2]" ProjectParameters.m4.json > ProjectParameters_30.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   $IS,   $IS,  $IS2,   0.0,   0.0]" ProjectParameters.m4.json > ProjectParameters_31.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   $IS,   $IS,   0.0,  $IS2,   0.0]" ProjectParameters.m4.json > ProjectParameters_32.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   $IS,   $IS,   0.0,   0.0,  $IS2]" ProjectParameters.m4.json > ProjectParameters_33.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   $IS,   0.0,  $IS2,  $IS2,   0.0]" ProjectParameters.m4.json > ProjectParameters_34.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   $IS,   0.0,  $IS2,   0.0,  $IS2]" ProjectParameters.m4.json > ProjectParameters_35.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   $IS,   0.0,   0.0,  $IS2,  $IS2]" ProjectParameters.m4.json > ProjectParameters_36.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   0.0,   $IS,  $IS2,  $IS2,   0.0]" ProjectParameters.m4.json > ProjectParameters_37.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   0.0,   $IS,  $IS2,   0.0,  $IS2]" ProjectParameters.m4.json > ProjectParameters_38.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   0.0,   $IS,   0.0,  $IS2,  $IS2]" ProjectParameters.m4.json > ProjectParameters_39.json
m4 -DM4VAR_INITIALSTRAIN="[0.0,   0.0,   0.0,  $IS2,  $IS2,  $IS2]" ProjectParameters.m4.json > ProjectParameters_40.json

for i in {00..40}
do
	TRAJ="trajectory_$i"
	echo $TRAJ 
	mkdir -p $TRAJ
	cp materials.json $TRAJ/materials.json
	cp model.mdpa $TRAJ/model.mdpa
	cp MainKratos.py $TRAJ/MainKratos.py
	mv ProjectParameters_$i.json $TRAJ/ProjectParameters.json
done
