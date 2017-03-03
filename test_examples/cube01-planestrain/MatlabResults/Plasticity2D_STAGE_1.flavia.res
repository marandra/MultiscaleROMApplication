GiD Post Results File 1.0 
# Datos de puntos de gauss del Set_1
GaussPoints "GP_Unico_Set_1" Elemtype Quadrilateral "Set_1"
Number of Gauss Points: 1
Nodes not included
Natural Coordinates: Internal
End GaussPoints
GaussPoints "GP_Set_1" Elemtype Quadrilateral "Set_1"
Number of Gauss Points: 4
Nodes not included
Natural Coordinates: Given
-0.577350 -0.577350
0.577350 -0.577350
-0.577350 0.577350
0.577350 0.577350
End GaussPoints
Result "Displacements//Fluctuations" "Load Analysis" 1 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0214285714285714
4 0 -0.0214285714285714
End Values
Result "Displacements//Total" "Load Analysis" 1 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0214285714285714
4 0 -0.0214285714285714
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 1 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0214285714285714
4 0 -0.0214285714285714
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 1 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 5.49450549450549 -5.00401002685617e-17 1.64835164835165 0
 5.49450549450549 -5.00401002685617e-17 1.64835164835165 0
 5.49450549450549 -1.00080200537123e-16 1.64835164835165 0
 5.49450549450549 -1.00080200537123e-16 1.64835164835165 0
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 1 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.05 -0.0214285714285714 -8.67361737988404e-19 0
 0.05 -0.0214285714285714 -8.67361737988404e-19 0
 0.05 -0.0214285714285714 -1.73472347597681e-18 0
 0.05 -0.0214285714285714 -1.73472347597681e-18 0
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 1 Scalar OnGaussPoints "GP_Set_1"
Values
1 0.137362637362637
 0.137362637362637
 0.137362637362637
 0.137362637362637
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 1 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 1 Scalar OnGaussPoints "GP_Set_1"
Values
1 4.88362330621736
 4.88362330621736
 4.88362330621736
 4.88362330621736
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 1 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 1 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 0
End Values
Result "Displacements//Fluctuations" "Load Analysis" 2 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.1 0
3 0.1 -0.0428571428571429
4 0 -0.0428571428571429
End Values
Result "Displacements//Total" "Load Analysis" 2 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.1 0
3 0.1 -0.0428571428571429
4 0 -0.0428571428571429
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 2 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0214285714285714
4 0 -0.0214285714285714
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 2 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 10.989010989011 -1.50120300805685e-16 3.2967032967033 -1.33440267382831e-16
 10.989010989011 -9.38218519968687e-16 3.2967032967033 -1.33440267382831e-16
 10.989010989011 -2.50200501342809e-16 3.2967032967033 -5.33761069531325e-16
 10.989010989011 -1.03829872050581e-15 3.2967032967033 -5.33761069531325e-16
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 2 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.1 -0.0428571428571429 -2.60208521396521e-18 -3.46944695195361e-18
 0.1 -0.0428571428571429 -8.67361737988404e-19 -3.46944695195361e-18
 0.1 -0.0428571428571429 -4.33680868994202e-18 -1.38777878078145e-17
 0.1 -0.0428571428571429 -2.60208521396521e-18 -1.38777878078145e-17
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 2 Scalar OnGaussPoints "GP_Set_1"
Values
1 0.549450549450549
 0.549450549450549
 0.549450549450549
 0.549450549450549
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 2 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 2 Scalar OnGaussPoints "GP_Set_1"
Values
1 9.76724661243471
 9.76724661243471
 9.76724661243471
 9.76724661243471
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 2 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 2 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 0
End Values
Result "Displacements//Fluctuations" "Load Analysis" 3 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.15 0
3 0.15 -0.0642857142857143
4 0 -0.0642857142857143
End Values
Result "Displacements//Total" "Load Analysis" 3 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.15 0
3 0.15 -0.0642857142857143
4 0 -0.0642857142857143
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 3 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0214285714285714
4 0 -0.0214285714285714
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 3 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 16.4835164835165 -2.00160401074247e-16 4.94505494505494 0
 16.4835164835165 -2.00160401074247e-16 4.94505494505494 0
 16.4835164835165 -3.0024060161137e-16 4.94505494505494 0
 16.4835164835165 -3.0024060161137e-16 4.94505494505494 0
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 3 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.15 -0.0642857142857143 -3.46944695195361e-18 0
 0.15 -0.0642857142857143 -3.46944695195361e-18 0
 0.15 -0.0642857142857143 -5.20417042793042e-18 0
 0.15 -0.0642857142857143 -5.20417042793042e-18 0
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 3 Scalar OnGaussPoints "GP_Set_1"
Values
1 1.23626373626374
 1.23626373626374
 1.23626373626374
 1.23626373626374
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 3 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 3 Scalar OnGaussPoints "GP_Set_1"
Values
1 14.6508699186521
 14.6508699186521
 14.6508699186521
 14.6508699186521
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 3 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 3 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 0
End Values
Result "Displacements//Fluctuations" "Load Analysis" 4 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.2 0
3 0.2 -0.0857142857142857
4 0 -0.0857142857142857
End Values
Result "Displacements//Total" "Load Analysis" 4 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.2 0
3 0.2 -0.0857142857142857
4 0 -0.0857142857142857
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 4 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0214285714285714
4 0 -0.0214285714285714
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 4 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 21.978021978022 -2.00160401074247e-16 6.59340659340659 0
 21.978021978022 -2.00160401074247e-16 6.59340659340659 0
 21.978021978022 -4.00320802148494e-16 6.59340659340659 0
 21.978021978022 -4.00320802148494e-16 6.59340659340659 0
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 4 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.2 -0.0857142857142857 -3.46944695195361e-18 0
 0.2 -0.0857142857142857 -3.46944695195361e-18 0
 0.2 -0.0857142857142857 -6.93889390390723e-18 0
 0.2 -0.0857142857142857 -6.93889390390723e-18 0
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 4 Scalar OnGaussPoints "GP_Set_1"
Values
1 2.1978021978022
 2.1978021978022
 2.1978021978022
 2.1978021978022
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 4 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 4 Scalar OnGaussPoints "GP_Set_1"
Values
1 19.5344932248694
 19.5344932248694
 19.5344932248694
 19.5344932248694
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 4 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 4 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 0
End Values
Result "Displacements//Fluctuations" "Load Analysis" 5 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.25 0
3 0.25 -0.107142857142857
4 0 -0.107142857142857
End Values
Result "Displacements//Total" "Load Analysis" 5 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.25 0
3 0.25 -0.107142857142857
4 0 -0.107142857142857
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 5 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0214285714285714
4 0 -0.0214285714285714
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 5 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 27.4725274725275 1.576196438326e-15 8.24175824175824 1.33440267382831e-16
 27.4725274725275 1.17587563617751e-15 8.24175824175824 1.33440267382831e-16
 27.4725274725275 1.576196438326e-15 8.24175824175824 1.06752213906265e-15
 27.4725274725275 1.17587563617751e-15 8.24175824175824 1.06752213906265e-15
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 5 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.25 -0.107142857142857 -3.46944695195361e-18 3.46944695195361e-18
 0.25 -0.107142857142857 -1.04083408558608e-17 3.46944695195361e-18
 0.25 -0.107142857142857 -3.46944695195361e-18 2.77555756156289e-17
 0.25 -0.107142857142857 -1.04083408558608e-17 2.77555756156289e-17
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 5 Scalar OnGaussPoints "GP_Set_1"
Values
1 3.43406593406593
 3.43406593406593
 3.43406593406593
 3.43406593406593
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 5 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 5 Scalar OnGaussPoints "GP_Set_1"
Values
1 24.4181165310868
 24.4181165310868
 24.4181165310868
 24.4181165310868
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 5 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 5 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 0
End Values
Result "Displacements//Fluctuations" "Load Analysis" 6 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.3 0
3 0.3 -0.128571428571429
4 0 -0.128571428571429
End Values
Result "Displacements//Total" "Load Analysis" 6 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.3 0
3 0.3 -0.128571428571429
4 0 -0.128571428571429
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 6 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0214285714285714
4 0 -0.0214285714285714
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 6 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 32.967032967033 -4.00320802148494e-16 9.89010989010989 0
 32.967032967033 -4.00320802148494e-16 9.89010989010989 0
 32.967032967033 -6.00481203222741e-16 9.89010989010989 0
 32.967032967033 -6.00481203222741e-16 9.89010989010989 0
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 6 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.3 -0.128571428571429 -6.93889390390723e-18 0
 0.3 -0.128571428571429 -6.93889390390723e-18 0
 0.3 -0.128571428571429 -1.04083408558608e-17 0
 0.3 -0.128571428571429 -1.04083408558608e-17 0
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 6 Scalar OnGaussPoints "GP_Set_1"
Values
1 4.94505494505494
 4.94505494505494
 4.94505494505494
 4.94505494505494
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 6 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 6 Scalar OnGaussPoints "GP_Set_1"
Values
1 29.3017398373041
 29.3017398373041
 29.3017398373041
 29.3017398373041
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 6 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 6 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 0
End Values
Result "Displacements//Fluctuations" "Load Analysis" 7 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.35 0
3 0.35 -0.15
4 0 -0.15
End Values
Result "Displacements//Total" "Load Analysis" 7 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.35 0
3 0.35 -0.15
4 0 -0.15
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 7 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0214285714285714
4 0 -0.0214285714285714
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 7 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 38.4615384615385 3.15239287665201e-15 11.5384615384615 0
 38.4615384615385 6.70510655545251e-15 11.5384615384615 0
 38.4615384615385 3.15239287665201e-15 11.5384615384615 0
 38.4615384615385 6.70510655545251e-15 11.5384615384615 0
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 7 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.35 -0.15 -6.93889390390723e-18 0
 0.35 -0.15 -6.93889390390723e-18 0
 0.35 -0.15 -6.93889390390723e-18 0
 0.35 -0.15 -6.93889390390723e-18 0
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 7 Scalar OnGaussPoints "GP_Set_1"
Values
1 6.73076923076923
 6.73076923076923
 6.73076923076923
 6.73076923076923
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 7 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 7 Scalar OnGaussPoints "GP_Set_1"
Values
1 34.1853631435215
 34.1853631435215
 34.1853631435215
 34.1853631435215
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 7 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 7 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 0
End Values
Result "Displacements//Fluctuations" "Load Analysis" 8 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.4 0
3 0.4 -0.171428571428571
4 0 -0.171428571428571
End Values
Result "Displacements//Total" "Load Analysis" 8 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.4 0
3 0.4 -0.171428571428571
4 0 -0.171428571428571
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 8 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0214285714285715
4 0 -0.0214285714285715
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 8 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 43.9560439560439 -4.00320802148494e-16 13.1868131868132 0
 43.9560439560439 -4.00320802148494e-16 13.1868131868132 0
 43.9560439560439 -8.00641604296988e-16 13.1868131868132 0
 43.9560439560439 -8.00641604296988e-16 13.1868131868132 0
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 8 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.4 -0.171428571428571 -6.93889390390723e-18 0
 0.4 -0.171428571428571 -6.93889390390723e-18 0
 0.4 -0.171428571428571 -1.38777878078145e-17 0
 0.4 -0.171428571428571 -1.38777878078145e-17 0
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 8 Scalar OnGaussPoints "GP_Set_1"
Values
1 8.79120879120879
 8.79120879120879
 8.79120879120879
 8.79120879120879
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 8 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 8 Scalar OnGaussPoints "GP_Set_1"
Values
1 39.0689864497388
 39.0689864497388
 39.0689864497388
 39.0689864497388
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 8 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 8 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 0
End Values
Result "Displacements//Fluctuations" "Load Analysis" 9 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.45 0
3 0.45 -0.192857142857143
4 0 -0.192857142857143
End Values
Result "Displacements//Total" "Load Analysis" 9 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.45 0
3 0.45 -0.192857142857143
4 0 -0.192857142857143
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 9 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0214285714285714
4 0 -0.0214285714285714
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 9 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 49.4505494505494 -8.00641604296988e-16 14.8351648351648 0
 49.4505494505494 -8.00641604296988e-16 14.8351648351648 0
 49.4505494505494 -4.00320802148494e-16 14.8351648351648 0
 49.4505494505494 3.15239287665201e-15 14.8351648351648 0
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 9 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.45 -0.192857142857143 -1.38777878078145e-17 0
 0.45 -0.192857142857143 -1.38777878078145e-17 0
 0.45 -0.192857142857143 -6.93889390390723e-18 0
 0.45 -0.192857142857143 -6.93889390390723e-18 0
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 9 Scalar OnGaussPoints "GP_Set_1"
Values
1 11.1263736263736
 11.1263736263736
 11.1263736263736
 11.1263736263736
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 9 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 9 Scalar OnGaussPoints "GP_Set_1"
Values
1 43.9526097559562
 43.9526097559562
 43.9526097559562
 43.9526097559562
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 9 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 9 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 0
End Values
Result "Displacements//Fluctuations" "Load Analysis" 10 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.5 0
3 0.5 -0.214285714285714
4 0 -0.214285714285714
End Values
Result "Displacements//Total" "Load Analysis" 10 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.5 0
3 0.5 -0.214285714285714
4 0 -0.214285714285714
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 10 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0214285714285714
4 0 -0.0214285714285714
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 10 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 54.9450549450549 -6.00481203222741e-16 16.4835164835165 2.66880534765663e-16
 54.9450549450549 2.55191167342927e-15 16.4835164835165 2.66880534765663e-16
 54.9450549450549 -6.00481203222741e-16 16.4835164835165 1.06752213906265e-15
 54.9450549450549 2.55191167342927e-15 16.4835164835165 1.06752213906265e-15
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 10 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.5 -0.214285714285714 -1.04083408558608e-17 6.93889390390723e-18
 0.5 -0.214285714285714 -1.73472347597681e-17 6.93889390390723e-18
 0.5 -0.214285714285714 -1.04083408558608e-17 2.77555756156289e-17
 0.5 -0.214285714285714 -1.73472347597681e-17 2.77555756156289e-17
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 10 Scalar OnGaussPoints "GP_Set_1"
Values
1 13.7362637362637
 13.7362637362637
 13.7362637362637
 13.7362637362637
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 10 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 10 Scalar OnGaussPoints "GP_Set_1"
Values
1 48.8362330621736
 48.8362330621736
 48.8362330621736
 48.8362330621736
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 10 Scalar OnGaussPoints "GP_Set_1"
Values
1 0
 0
 0
 0
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 10 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 0
End Values
Result "Displacements//Fluctuations" "Load Analysis" 11 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.55 0
3 0.55 -0.251258614094061
4 0 -0.251258614094061
End Values
Result "Displacements//Total" "Load Analysis" 11 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.55 0
3 0.55 -0.251258614094061
4 0 -0.251258614094061
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 11 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0369728998083472
4 0 -0.0369728998083471
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 11 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 56.8563398376856 7.50510764646606e-14 17.829006638799 0
 56.8563398376856 7.50510764646606e-14 17.829006638799 0
 56.8563398376856 7.50510764646606e-14 17.829006638799 0
 56.8563398376856 7.50510764646606e-14 17.829006638799 0
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 11 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.55 -0.251258614094061 -2.08166817117217e-17 0
 0.55 -0.251258614094061 -2.08166817117217e-17 0
 0.55 -0.251258614094061 -1.38777878078145e-17 0
 0.55 -0.251258614094061 -1.38777878078145e-17 0
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 11 Scalar OnGaussPoints "GP_Set_1"
Values
1 14.7249697992164
 14.7249697992164
 14.7249697992164
 14.7249697992164
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 11 Scalar OnGaussPoints "GP_Set_1"
Values
1 0.0366901801091756
 0.0366901801091756
 0.0366901801091756
 0.0366901801091756
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 11 Scalar OnGaussPoints "GP_Set_1"
Values
1 50.3669018010917
 50.3669018010917
 50.3669018010917
 50.3669018010917
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 11 Scalar OnGaussPoints "GP_Set_1"
Values
1 1
 1
 1
 1
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 11 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 1
End Values
Result "Displacements//Fluctuations" "Load Analysis" 12 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.6 0
3 0.6 -0.293568895241909
4 0 -0.29356889524191
End Values
Result "Displacements//Total" "Load Analysis" 12 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.6 0
3 0.6 -0.293568895241909
4 0 -0.29356889524191
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 12 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0423102811478481
4 0 -0.0423102811478481
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 12 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 57.6148746497015 1.28785870856518e-13 18.992901539821 4.80831093885203e-16
 57.6148746497015 1.29674049276218e-13 18.992901539821 4.80831093885203e-16
 57.6148746497015 1.29674049276218e-13 18.992901539821 9.61662187770406e-16
 57.6148746497015 1.29674049276218e-13 18.992901539821 9.61662187770406e-16
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 12 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.6 -0.29356889524191 -6.93889390390723e-18 1.38777878078145e-17
 0.6 -0.29356889524191 -2.08166817117217e-17 1.38777878078145e-17
 0.6 -0.29356889524191 -1.38777878078145e-17 2.77555756156289e-17
 0.6 -0.29356889524191 -2.77555756156289e-17 2.77555756156289e-17
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 12 Scalar OnGaussPoints "GP_Set_1"
Values
1 15.1729975590058
 15.1729975590058
 15.1729975590058
 15.1729975590058
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 12 Scalar OnGaussPoints "GP_Set_1"
Values
1 0.0852044682109183
 0.0852044682109184
 0.0852044682109184
 0.0852044682109184
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 12 Scalar OnGaussPoints "GP_Set_1"
Values
1 50.8520446821092
 50.8520446821092
 50.8520446821092
 50.8520446821092
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 12 Scalar OnGaussPoints "GP_Set_1"
Values
1 1
 1
 1
 1
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 12 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 1
End Values
Result "Displacements//Fluctuations" "Load Analysis" 13 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.65 0
3 0.65 -0.336311082444734
4 0 -0.336311082444734
End Values
Result "Displacements//Total" "Load Analysis" 13 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.65 0
3 0.65 -0.336311082444734
4 0 -0.336311082444734
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 13 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.0499999999999999 0
3 0.0499999999999999 -0.0427421872028249
4 0 -0.0427421872028248
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 13 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 58.3450831183217 1.19904086659517e-13 20.0771462704945 -4.76950859518665e-17
 58.3450831183217 1.26121335597418e-13 20.0771462704945 -4.76950859518666e-17
 58.3450831183217 1.19015908239817e-13 20.0771462704945 -9.53901719037332e-17
 58.3450831183217 1.26121335597418e-13 20.0771462704945 -9.53901719037332e-17
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 13 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.65 -0.336311082444734 -1.38777878078145e-17 0
 0.65 -0.336311082444734 -1.38777878078145e-17 0
 0.65 -0.336311082444734 -1.38777878078145e-17 0
 0.65 -0.336311082444734 -1.38777878078145e-17 0
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 13 Scalar OnGaussPoints "GP_Set_1"
Values
1 15.6422534444483
 15.6422534444483
 15.6422534444483
 15.6422534444483
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 13 Scalar OnGaussPoints "GP_Set_1"
Values
1 0.13404105799709
 0.13404105799709
 0.13404105799709
 0.13404105799709
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 13 Scalar OnGaussPoints "GP_Set_1"
Values
1 51.3404105799709
 51.3404105799709
 51.3404105799709
 51.3404105799709
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 13 Scalar OnGaussPoints "GP_Set_1"
Values
1 1
 1
 1
 1
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 13 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 1
End Values
Result "Displacements//Fluctuations" "Load Analysis" 14 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.7 0
3 0.7 -0.379438393853138
4 0 -0.379438393853138
End Values
Result "Displacements//Total" "Load Analysis" 14 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.7 0
3 0.7 -0.379438393853138
4 0 -0.379438393853138
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 14 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0431273114084033
4 0 -0.0431273114084033
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 14 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 59.0517125917704 1.08801856413265e-13 21.0886889449451 4.38155831563084e-16
 59.0517125917704 1.16351372980716e-13 21.0886889449451 4.38155831563084e-16
 59.0517125917704 1.08801856413265e-13 21.0886889449451 1.8386111581628e-15
 59.0517125917704 1.16351372980716e-13 21.0886889449451 1.8386111581628e-15
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 14 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.7 -0.379438393853138 -6.93889390390723e-18 1.38777878078145e-17
 0.7 -0.379438393853138 -2.08166817117217e-17 1.38777878078145e-17
 0.7 -0.379438393853138 -6.93889390390723e-18 5.55111512312578e-17
 0.7 -0.379438393853138 -2.08166817117217e-17 5.55111512312578e-17
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 14 Scalar OnGaussPoints "GP_Set_1"
Values
1 16.1334426953013
 16.1334426953013
 16.1334426953013
 16.1334426953013
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 14 Scalar OnGaussPoints "GP_Set_1"
Values
1 0.183159618344348
 0.183159618344348
 0.183159618344348
 0.183159618344348
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 14 Scalar OnGaussPoints "GP_Set_1"
Values
1 51.8315961834435
 51.8315961834435
 51.8315961834435
 51.8315961834435
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 14 Scalar OnGaussPoints "GP_Set_1"
Values
1 1
 1
 1
 1
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 14 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 1
End Values
Result "Displacements//Fluctuations" "Load Analysis" 15 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.75 0
3 0.75 -0.422908947626093
4 0 -0.422908947626093
End Values
Result "Displacements//Total" "Load Analysis" 15 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.75 0
3 0.75 -0.422908947626093
4 0 -0.422908947626093
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 15 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.0499999999999999 0
3 0.0499999999999999 -0.0434705537729555
4 0 -0.0434705537729554
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 15 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 59.7387505547889 9.1926466438963e-14 22.0340125386877 -8.62195973877955e-17
 59.738750554789 9.41469124882133e-14 22.0340125386877 -8.62195973877957e-17
 59.7387505547889 9.10382880192628e-14 22.0340125386877 -2.67332049987312e-16
 59.738750554789 9.41469124882133e-14 22.0340125386877 -2.67332049987313e-16
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 15 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.75 -0.422908947626093 -2.77555756156289e-17 0
 0.75 -0.422908947626093 -2.77555756156289e-17 0
 0.75 -0.422908947626093 -1.38777878078145e-17 0
 0.75 -0.422908947626093 -1.38777878078145e-17 0
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 15 Scalar OnGaussPoints "GP_Set_1"
Values
1 16.6472379561323
 16.6472379561323
 16.6472379561323
 16.6472379561323
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 15 Scalar OnGaussPoints "GP_Set_1"
Values
1 0.232524866288339
 0.232524866288339
 0.232524866288339
 0.232524866288339
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 15 Scalar OnGaussPoints "GP_Set_1"
Values
1 52.3252486628834
 52.3252486628834
 52.3252486628834
 52.3252486628834
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 15 Scalar OnGaussPoints "GP_Set_1"
Values
1 1
 1
 1
 1
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 15 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 1
End Values
Result "Displacements//Fluctuations" "Load Analysis" 16 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.8 0
3 0.8 -0.466685372496863
4 0 -0.466685372496863
End Values
Result "Displacements//Total" "Load Analysis" 16 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.8 0
3 0.8 -0.466685372496863
4 0 -0.466685372496863
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 16 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0437764248707698
4 0 -0.0437764248707698
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 16 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 60.4095334421279 8.70414851306123e-14 22.9191234336563 -7.77939512718408e-17
 60.4095334421279 8.70414851306123e-14 22.9191234336563 -7.7793951271841e-17
 60.4095334421279 9.45910016980633e-14 22.9191234336563 -2.4120753402008e-16
 60.4095334421279 8.39328606616618e-14 22.9191234336563 -2.41207534020081e-16
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 16 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.8 -0.466685372496863 -1.38777878078145e-17 0
 0.8 -0.466685372496863 -1.38777878078145e-17 0
 0.8 -0.466685372496863 -2.77555756156289e-17 0
 0.8 -0.466685372496863 -2.77555756156289e-17 0
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 16 Scalar OnGaussPoints "GP_Set_1"
Values
1 17.184271453372
 17.184271453372
 17.184271453372
 17.184271453372
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 16 Scalar OnGaussPoints "GP_Set_1"
Values
1 0.282106015722343
 0.282106015722343
 0.282106015722343
 0.282106015722343
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 16 Scalar OnGaussPoints "GP_Set_1"
Values
1 52.8210601572234
 52.8210601572234
 52.8210601572234
 52.8210601572234
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 16 Scalar OnGaussPoints "GP_Set_1"
Values
1 1
 1
 1
 1
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 16 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 1
End Values
Result "Displacements//Fluctuations" "Load Analysis" 17 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.85 0
3 0.85 -0.51073440494953
4 0 -0.51073440494953
End Values
Result "Displacements//Total" "Load Analysis" 17 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.85 0
3 0.85 -0.51073440494953
4 0 -0.51073440494953
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 17 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0440490324526673
4 0 -0.0440490324526673
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 17 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 61.066844072651 7.01660951563099e-14 23.7495546899664 -7.02298188606913e-17
 61.066844072651 7.01660951563099e-14 23.7495546899664 -7.02298188606915e-17
 61.066844072651 7.01660951563099e-14 23.7495546899664 -2.17754223112666e-16
 61.066844072651 7.01660951563099e-14 23.7495546899664 -2.17754223112666e-16
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 17 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.85 -0.51073440494953 -2.77555756156289e-17 0
 0.85 -0.51073440494953 -2.77555756156289e-17 0
 0.85 -0.51073440494953 -2.77555756156289e-17 0
 0.85 -0.51073440494953 -2.77555756156289e-17 0
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 17 Scalar OnGaussPoints "GP_Set_1"
Values
1 17.7451308777475
 17.7451308777475
 17.7451308777475
 17.7451308777475
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 17 Scalar OnGaussPoints "GP_Set_1"
Values
1 0.3318762550514
 0.3318762550514
 0.331876255051401
 0.331876255051401
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 17 Scalar OnGaussPoints "GP_Set_1"
Values
1 53.318762550514
 53.318762550514
 53.318762550514
 53.318762550514
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 17 Scalar OnGaussPoints "GP_Set_1"
Values
1 1
 1
 1
 1
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 17 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 1
End Values
Result "Displacements//Fluctuations" "Load Analysis" 18 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.9 0
3 0.9 -0.55502649089499
4 0 -0.55502649089499
End Values
Result "Displacements//Total" "Load Analysis" 18 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.9 0
3 0.9 -0.55502649089499
4 0 -0.55502649089499
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 18 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.0499999999999999 0
3 0.0499999999999999 -0.0442920859454595
4 0 -0.0442920859454595
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 18 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 61.712997474317 5.24025267623074e-14 24.5303798019355 -6.34380484008742e-17
 61.712997474317 5.24025267623074e-14 24.5303798019355 -6.34380484008744e-17
 61.712997474317 5.19584375524573e-14 24.5303798019355 -1.96695693786673e-16
 61.712997474317 5.50670620214078e-14 24.5303798019355 -1.96695693786673e-16
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 18 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.9 -0.55502649089499 -2.77555756156289e-17 0
 0.9 -0.55502649089499 -2.77555756156289e-17 0
 0.9 -0.55502649089499 -1.38777878078145e-17 0
 0.9 -0.55502649089499 -1.38777878078145e-17 0
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 18 Scalar OnGaussPoints "GP_Set_1"
Values
1 18.3303579275
 18.3303579275
 18.3303579275
 18.3303579275
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 18 Scalar OnGaussPoints "GP_Set_1"
Values
1 0.381812263291402
 0.381812263291402
 0.381812263291402
 0.381812263291402
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 18 Scalar OnGaussPoints "GP_Set_1"
Values
1 53.818122632914
 53.818122632914
 53.818122632914
 53.818122632914
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 18 Scalar OnGaussPoints "GP_Set_1"
Values
1 1
 1
 1
 1
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 18 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 1
End Values
Result "Displacements//Fluctuations" "Load Analysis" 19 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.95 0
3 0.95 -0.599535403868835
4 0 -0.599535403868835
End Values
Result "Displacements//Total" "Load Analysis" 19 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.95 0
3 0.95 -0.599535403868835
4 0 -0.599535403868835
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 19 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.05 0
3 0.05 -0.0445089129738455
4 0 -0.0445089129738455
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 19 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 62.349915485641 5.41788836017076e-14 25.2662335471501 -5.73382122502455e-17
 62.349915485641 5.41788836017076e-14 25.2662335471501 -5.73382122502457e-17
 62.349915485641 5.41788836017076e-14 25.2662335471501 -1.77782572499411e-16
 62.349915485641 5.41788836017076e-14 25.2662335471501 -1.77782572499411e-16
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 19 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 0.95 -0.599535403868835 -1.38777878078145e-17 0
 0.95 -0.599535403868835 -1.38777878078145e-17 0
 0.95 -0.599535403868835 -2.77555756156289e-17 0
 0.95 -0.599535403868835 -2.77555756156289e-17 0
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 19 Scalar OnGaussPoints "GP_Set_1"
Values
1 18.9404486954961
 18.9404486954961
 18.9404486954961
 18.9404486954961
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 19 Scalar OnGaussPoints "GP_Set_1"
Values
1 0.431893769598614
 0.431893769598614
 0.431893769598614
 0.431893769598614
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 19 Scalar OnGaussPoints "GP_Set_1"
Values
1 54.3189376959862
 54.3189376959862
 54.3189376959862
 54.3189376959862
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 19 Scalar OnGaussPoints "GP_Set_1"
Values
1 1
 1
 1
 1
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 19 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 1
End Values
Result "Displacements//Fluctuations" "Load Analysis" 20 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 1 0
3 1 -0.644237887225536
4 0 -0.644237887225536
End Values
Result "Displacements//Total" "Load Analysis" 20 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 1 0
3 1 -0.644237887225536
4 0 -0.644237887225536
End Values
Result "Displacements//Incremental Fluct." "Load Analysis" 20 Vector OnNodes
ComponentNames "X-DISPL" "Y-DISPL"
Values
1 0 0
2 0.0499999999999999 0
3 0.0499999999999999 -0.0447024833567012
4 0 -0.0447024833567012
End Values
Result "Stresses//On Gauss Points" "Load Analysis" 20 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Stress XX" "Stress YY" "Stress ZZ" "Stress XY"
Values
1 62.9791909026097 6.3504757008559e-14 25.9613372910061 -5.185795414862e-17
 62.9791909026097 6.3504757008559e-14 25.9613372910061 -5.18579541486201e-17
 62.9791909026097 6.3504757008559e-14 25.9613372910061 -1.60790511794491e-16
 62.9791909026097 6.3504757008559e-14 25.9613372910061 -1.60790511794492e-16
End Values
Result "StrainsFluct//On Gauss Points" "Load Analysis" 20 PlainDeformationMatrix OnGaussPoints "GP_Set_1"
ComponentNames "Strain XX" "Strain YY" "Strain ZZ" "Strain XY"
Values
1 1 -0.644237887225536 -2.77555756156289e-17 0
 1 -0.644237887225536 -2.77555756156289e-17 0
 1 -0.644237887225536 -2.77555756156289e-17 0
 1 -0.644237887225536 -2.77555756156289e-17 0
End Values
Result "Energy//GPs//On (1st-4th) Gauss Points" "Load Analysis" 20 Scalar OnGaussPoints "GP_Set_1"
Values
1 19.5758552785725
 19.5758552785725
 19.5758552785725
 19.5758552785725
End Values
Result  "Constitutive Model//Plastic Equiv. Strain//Sobre punto de Gauss" "Load Analysis" 20 Scalar OnGaussPoints "GP_Set_1"
Values
1 0.482103157865078
 0.482103157865078
 0.482103157865078
 0.482103157865078
End Values
Result  "Constitutive Model//Norm of deviatoric Stress//Sobre punto de Gauss" "Load Analysis" 20 Scalar OnGaussPoints "GP_Set_1"
Values
1 54.8210315786508
 54.8210315786508
 54.8210315786508
 54.8210315786508
End Values
Result  "Constitutive Model//Load Index//On Gauss Points" "Load Analysis" 20 Scalar OnGaussPoints "GP_Set_1"
Values
1 1
 1
 1
 1
End Values
Result  "Constitutive Model//Load Index//Sobre elemento (any PG)" "Load Analysis" 20 Scalar OnGaussPoints "GP_Unico_Set_1"
Values
1 1
End Values
