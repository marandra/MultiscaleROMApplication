m4 -DM4VAR_INITIALSTRAIN="[0.02, 0.00, 0.00, 0.00, 0.00, 0.00]" t_ProjectParameters.m4 > ProjectParameters_0.json
m4 -DM4VAR_INITIALSTRAIN="[0.00, 0.02, 0.00, 0.00, 0.00, 0.00]" t_ProjectParameters.m4 > ProjectParameters_1.json
m4 -DM4VAR_INITIALSTRAIN="[0.00, 0.00, 0.02, 0.00, 0.00, 0.00]" t_ProjectParameters.m4 > ProjectParameters_2.json
m4 -DM4VAR_INITIALSTRAIN="[0.00, 0.00, 0.00, 0.02, 0.00, 0.00]" t_ProjectParameters.m4 > ProjectParameters_3.json
m4 -DM4VAR_INITIALSTRAIN="[0.00, 0.00, 0.00, 0.00, 0.02, 0.00]" t_ProjectParameters.m4 > ProjectParameters_4.json
m4 -DM4VAR_INITIALSTRAIN="[0.00, 0.00, 0.00, 0.00, 0.00, 0.02]" t_ProjectParameters.m4 > ProjectParameters_5.json
m4 -DM4VAR_INITIALSTRAIN="[0.02, 0.02, 0.02, 0.00, 0.00, 0.00]" t_ProjectParameters.m4 > ProjectParameters_6.json
m4 -DM4VAR_INITIALSTRAIN="[0.00, 0.00, 0.00, 0.02, 0.02, 0.02]" t_ProjectParameters.m4 > ProjectParameters_7.json
m4 -DM4VAR_INITIALSTRAIN="[0.02, 0.02, 0.02, 0.02, 0.02, 0.02]" t_ProjectParameters.m4 > ProjectParameters_8.json

for i in 0 1 2 3 4 5 6 7 8
do
	TRAJ="trajectory_$i"
	echo $TRAJ 
	mkdir -p $TRAJ
	cp materials.py $TRAJ
	cp t_3d_27elem.mdpa $TRAJ/3d_27elem.mdpa
	cp t_MainKratos.py $TRAJ/MainKratos.py
	mv ProjectParameters_$i.json $TRAJ/ProjectParameters.json
done
