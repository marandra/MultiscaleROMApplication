IS=0.001
m4 -DM4VAR_INITIALSTRAIN="[$IS, 0.0, 0.0, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_00.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, $IS, 0.0, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_01.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, $IS, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_02.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, $IS, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_03.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, 0.0, $IS, 0.0]" t_ProjectParameters.m4 > ProjectParameters_04.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, 0.0, 0.0, $IS]" t_ProjectParameters.m4 > ProjectParameters_05.json
m4 -DM4VAR_INITIALSTRAIN="[$IS, $IS, 0.0, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_06.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, $IS, $IS, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_07.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, $IS, $IS, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_08.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, $IS, $IS, 0.0]" t_ProjectParameters.m4 > ProjectParameters_09.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, 0.0, $IS, $IS]" t_ProjectParameters.m4 > ProjectParameters_10.json
m4 -DM4VAR_INITIALSTRAIN="[$IS, $IS, $IS, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_11.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, $IS, $IS, $IS]" t_ProjectParameters.m4 > ProjectParameters_12.json

for i in 00 01 02 03 04 05 06 07 08 09 10 11 12
do
	TRAJ="trajectory_$i"
	echo $TRAJ 
	mkdir -p $TRAJ
	cp t_materials.json $TRAJ/materials.json
	cp t_model.mdpa $TRAJ/model.mdpa
	cp t_MainKratos.py $TRAJ/MainKratos.py
	mv ProjectParameters_$i.json $TRAJ/ProjectParameters.json
done
