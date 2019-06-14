#initialize this as a bash script
!/bin/bash
#define our variable ranges/sequences
for i in .1 .2 .3 .4 .5
do
      # replace sActsS with i, testingS with j, and adhS with k
      # creates a copy of the file so that you don't replace those variable names in the original, allowing you to change them with each subsequent loop
  sed "s/mortcalchange/$i/g" calibration_RI.py > caltest.py
      # run subTitan with your flags
      ~/./data/shared/OD_Modelling/subTitan.sh caltest.py -j calibration_mort_longerburn_$i -N 55000 -t 300 -n 10
done

 
