#!/bin/bash
IS=0.0005
m4 -DM4VAR_INITIALSTRAIN="[$IS, 0.0, 0.0, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_00.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, $IS, 0.0, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_01.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, $IS, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_02.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, $IS, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_03.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, 0.0, $IS, 0.0]" t_ProjectParameters.m4 > ProjectParameters_04.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, 0.0, 0.0, $IS]" t_ProjectParameters.m4 > ProjectParameters_05.json
m4 -DM4VAR_INITIALSTRAIN="[$IS, $IS, 0.0, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_06.json
m4 -DM4VAR_INITIALSTRAIN="[$IS, 0.0, $IS, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_07.json
m4 -DM4VAR_INITIALSTRAIN="[$IS, 0.0, 0.0, $IS, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_08.json
m4 -DM4VAR_INITIALSTRAIN="[$IS, 0.0, 0.0, 0.0, $IS, 0.0]" t_ProjectParameters.m4 > ProjectParameters_09.json
m4 -DM4VAR_INITIALSTRAIN="[$IS, 0.0, 0.0, 0.0, 0.0, $IS]" t_ProjectParameters.m4 > ProjectParameters_10.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, $IS, $IS, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_11.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, $IS, 0.0, $IS, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_12.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, $IS, 0.0, 0.0, $IS, 0.0]" t_ProjectParameters.m4 > ProjectParameters_13.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, $IS, 0.0, 0.0, 0.0, $IS]" t_ProjectParameters.m4 > ProjectParameters_14.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, $IS, $IS, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_15.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, $IS, 0.0, $IS, 0.0]" t_ProjectParameters.m4 > ProjectParameters_16.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, $IS, 0.0, 0.0, $IS]" t_ProjectParameters.m4 > ProjectParameters_17.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, $IS, $IS, 0.0]" t_ProjectParameters.m4 > ProjectParameters_18.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, $IS, 0.0, $IS]" t_ProjectParameters.m4 > ProjectParameters_19.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, 0.0, $IS, $IS]" t_ProjectParameters.m4 > ProjectParameters_20.json
#m4 -DM4VAR_INITIALSTRAIN="[$IS, $IS, $IS, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_21.json
#m4 -DM4VAR_INITIALSTRAIN="[$IS, $IS, 0.0, $IS, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_22.json
#m4 -DM4VAR_INITIALSTRAIN="[$IS, $IS, 0.0, 0.0, $IS, 0.0]" t_ProjectParameters.m4 > ProjectParameters_23.json
#m4 -DM4VAR_INITIALSTRAIN="[$IS, $IS, 0.0, 0.0, 0.0, $IS]" t_ProjectParameters.m4 > ProjectParameters_24.json
#m4 -DM4VAR_INITIALSTRAIN="[$IS, 0.0, $IS, $IS, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_25.json
#m4 -DM4VAR_INITIALSTRAIN="[$IS, 0.0, $IS, 0.0, $IS, 0.0]" t_ProjectParameters.m4 > ProjectParameters_26.json
#m4 -DM4VAR_INITIALSTRAIN="[$IS, 0.0, $IS, 0.0, 0.0, $IS]" t_ProjectParameters.m4 > ProjectParameters_27.json
#m4 -DM4VAR_INITIALSTRAIN="[$IS, 0.0, 0.0, $IS, $IS, 0.0]" t_ProjectParameters.m4 > ProjectParameters_28.json
#m4 -DM4VAR_INITIALSTRAIN="[$IS, 0.0, 0.0, $IS, 0.0, $IS]" t_ProjectParameters.m4 > ProjectParameters_29.json
#m4 -DM4VAR_INITIALSTRAIN="[$IS, 0.0, 0.0, 0.0, $IS, $IS]" t_ProjectParameters.m4 > ProjectParameters_30.json
#m4 -DM4VAR_INITIALSTRAIN="[0.0, $IS, $IS, $IS, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_31.json
#m4 -DM4VAR_INITIALSTRAIN="[0.0, $IS, $IS, 0.0, $IS, 0.0]" t_ProjectParameters.m4 > ProjectParameters_32.json
#m4 -DM4VAR_INITIALSTRAIN="[0.0, $IS, $IS, 0.0, 0.0, $IS]" t_ProjectParameters.m4 > ProjectParameters_33.json
#m4 -DM4VAR_INITIALSTRAIN="[0.0, $IS, 0.0, $IS, $IS, 0.0]" t_ProjectParameters.m4 > ProjectParameters_34.json
#m4 -DM4VAR_INITIALSTRAIN="[0.0, $IS, 0.0, $IS, 0.0, $IS]" t_ProjectParameters.m4 > ProjectParameters_35.json
#m4 -DM4VAR_INITIALSTRAIN="[0.0, $IS, 0.0, 0.0, $IS, $IS]" t_ProjectParameters.m4 > ProjectParameters_36.json
#m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, $IS, $IS, $IS, 0.0]" t_ProjectParameters.m4 > ProjectParameters_37.json
#m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, $IS, $IS, 0.0, $IS]" t_ProjectParameters.m4 > ProjectParameters_38.json
#m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, $IS, 0.0, $IS, $IS]" t_ProjectParameters.m4 > ProjectParameters_39.json
#m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, $IS, $IS, $IS]" t_ProjectParameters.m4 > ProjectParameters_40.json

for i in {00..20}
do
	TRAJ="trajectory_$i"
	echo $TRAJ 
	mkdir -p $TRAJ
	cp t_materials.json $TRAJ/materials.json
	cp t_model.mdpa $TRAJ/model.mdpa
	cp t_MainKratos.py $TRAJ/MainKratos.py
	mv ProjectParameters_$i.json $TRAJ/ProjectParameters.json
done
