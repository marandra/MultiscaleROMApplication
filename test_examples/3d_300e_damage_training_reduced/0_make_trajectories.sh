m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, 0.0, 0.0, 0.001]" t_ProjectParameters.m4 > ProjectParameters_00.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, 0.0, 0.001, 0.0]" t_ProjectParameters.m4 > ProjectParameters_01.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, 0.0, 0.001, 0.001]" t_ProjectParameters.m4 > ProjectParameters_02.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, 0.001, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_03.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, 0.001, 0.0, 0.001]" t_ProjectParameters.m4 > ProjectParameters_04.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, 0.001, 0.001, 0.0]" t_ProjectParameters.m4 > ProjectParameters_05.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.001, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_07.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.001, 0.0, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_15.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.001, 0.001, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_23.json
m4 -DM4VAR_INITIALSTRAIN="[0.001, 0.0, 0.0, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_31.json
m4 -DM4VAR_INITIALSTRAIN="[0.001, 0.0, 0.001, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_39.json
m4 -DM4VAR_INITIALSTRAIN="[0.001, 0.001, 0.0, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_47.json

for i in 00 01 02 03 04 05 07 15 23 31 39 47
do
	TRAJ="trajectory_$i"
	echo $TRAJ 
	mkdir -p $TRAJ
	cp t_materials.json $TRAJ/materials.json
	cp t_model.mdpa $TRAJ/model.mdpa
	cp t_MainKratos.py $TRAJ/MainKratos.py
	mv ProjectParameters_$i.json $TRAJ/ProjectParameters.json
done

