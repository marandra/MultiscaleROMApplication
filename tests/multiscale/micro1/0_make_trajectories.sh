STRA="0.000001" 
m4 -DM4VAR_INITIALSTRAIN="[0.0  , 0.0  , 0.0  , 0.0  , 0.0  , $STRA]" micro_ProjectParameters_m4.json > ProjectParameters_00.json
m4 -DM4VAR_INITIALSTRAIN="[0.0  , 0.0  , 0.0  , 0.0  , $STRA, 0.0  ]" micro_ProjectParameters_m4.json > ProjectParameters_01.json
m4 -DM4VAR_INITIALSTRAIN="[0.0  , 0.0  , 0.0  , 0.0  , $STRA, $STRA]" micro_ProjectParameters_m4.json > ProjectParameters_02.json
m4 -DM4VAR_INITIALSTRAIN="[0.0  , 0.0  , 0.0  , $STRA, 0.0  , 0.0  ]" micro_ProjectParameters_m4.json > ProjectParameters_03.json
m4 -DM4VAR_INITIALSTRAIN="[0.0  , 0.0  , 0.0  , $STRA, 0.0  , $STRA]" micro_ProjectParameters_m4.json > ProjectParameters_04.json
m4 -DM4VAR_INITIALSTRAIN="[0.0  , 0.0  , 0.0  , $STRA, $STRA, 0.0  ]" micro_ProjectParameters_m4.json > ProjectParameters_05.json
m4 -DM4VAR_INITIALSTRAIN="[0.0  , 0.0  , $STRA, 0.0  , 0.0  , 0.0  ]" micro_ProjectParameters_m4.json > ProjectParameters_06.json
m4 -DM4VAR_INITIALSTRAIN="[0.0  , $STRA, 0.0  , 0.0  , 0.0  , 0.0  ]" micro_ProjectParameters_m4.json > ProjectParameters_07.json
m4 -DM4VAR_INITIALSTRAIN="[0.0  , $STRA, $STRA, 0.0  , 0.0  , 0.0  ]" micro_ProjectParameters_m4.json > ProjectParameters_08.json
m4 -DM4VAR_INITIALSTRAIN="[$STRA, 0.0  , 0.0  , 0.0  , 0.0  , 0.0  ]" micro_ProjectParameters_m4.json > ProjectParameters_09.json
m4 -DM4VAR_INITIALSTRAIN="[$STRA, 0.0  , $STRA, 0.0  , 0.0  , 0.0  ]" micro_ProjectParameters_m4.json > ProjectParameters_10.json
m4 -DM4VAR_INITIALSTRAIN="[$STRA, $STRA, 0.0  , 0.0  , 0.0  , 0.0  ]" micro_ProjectParameters_m4.json > ProjectParameters_11.json

for i in 00 01 02 03 04 05 06 07 08 09 10 11
do
	TRAJ="trajectory_$i"
	echo $TRAJ 
	mkdir -p $TRAJ
	cp micro_materials.json $TRAJ/
	cp micro_model.mdpa $TRAJ/
	cp micro_MainKratos.py $TRAJ/MainKratos.py
	mv ProjectParameters_$i.json $TRAJ/ProjectParameters.json
done

