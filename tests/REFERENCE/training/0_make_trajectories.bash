#!/bin/bash
STRAIN_SET="training_strain_set.dat"
NR_TRAJ=`wc $STRAIN_SET | awk {'print $1'}`
LAST_TRAJ=`expr $NR_TRAJ - 1` 
LEN=${#LAST_TRAJ}
for ((i=0; i<$NR_TRAJ; i++))
do
        printf -v ID "%0${LEN}d" $i
	TRAJ="trajectory_"$ID
        LINE=`expr $ID + 1`
        STRAIN_VECTOR=`head -$LINE $STRAIN_SET | tail -1 | awk {'print "["$1",   "$2",   "$3",   "$4",   "$5",   "$6"]"'}`
        echo $ID": "$STRAIN_VECTOR
        mkdir -p $TRAJ
        m4 -DM4VAR_INITIALSTRAIN="${STRAIN_VECTOR}" ProjectParameters.m4.json > ProjectParameters.json
        mv ProjectParameters.json $TRAJ
	cp materials.json model.mdpa MainKratos.py $TRAJ
done
