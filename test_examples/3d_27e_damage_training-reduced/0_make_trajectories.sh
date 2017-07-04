IS=0.001
m4 -DM4VAR_INITIALSTRAIN="[$IS, 0.0, 0.0, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_0.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, $IS, 0.0, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_1.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, $IS, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_2.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, $IS, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_3.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, 0.0, $IS, 0.0]" t_ProjectParameters.m4 > ProjectParameters_4.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, 0.0, 0.0, $IS]" t_ProjectParameters.m4 > ProjectParameters_5.json
m4 -DM4VAR_INITIALSTRAIN="[$IS, $IS, $IS, 0.0, 0.0, 0.0]" t_ProjectParameters.m4 > ProjectParameters_6.json
m4 -DM4VAR_INITIALSTRAIN="[0.0, 0.0, 0.0, $IS, $IS, $IS]" t_ProjectParameters.m4 > ProjectParameters_7.json
m4 -DM4VAR_INITIALSTRAIN="[$IS, $IS, $IS, $IS, $IS, $IS]" t_ProjectParameters.m4 > ProjectParameters_8.json

for i in 0 1 2 3 4 5 6 7 8
do
	TRAJ="trajectory_$i"
	echo $TRAJ 
	mkdir -p $TRAJ
	cp materials.py $TRAJ
	cp t_model.mdpa $TRAJ/model.mdpa
	cp t_MainKratos.py $TRAJ/MainKratos.py
	mv ProjectParameters_$i.json $TRAJ/ProjectParameters.json
done
