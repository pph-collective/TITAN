#initialize this as a bash script
#!/bin/bash
#define our variable ranges/sequences
for i in .1
do
  for j in .68
  do
    for k in .01 .05 1 .005
     do
         # replace sActsS with i, testingS with j, and adhS with k
         # creates a copy of the file so that you don't replace those variable names in the original, allowing you to change them with each subsequent loop
         sed "s/MAT_HF_WHITE/$i/g" test.py > OD_cal.py
         sed -i "s/MAT_HM_WHITE/$i/g" OD_cal.py
         sed -i "s/MAT_HM_BLACK/$j/g" OD_cal.py
         sed -i "s/MAT_HF_BLACK/$j/g" OD_cal.py
         sed -i "s/calParam/$k/g" OD_cal.py
         # run subTitan with your flags
         ./subTitan.sh OD_cal.py -j cal_$i\_$j\_$k -N 55000 -s 0 -t 300 -n 1
     done
  done
done

 
