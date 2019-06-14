#initialize this as a bash script
#!/bin/bash
#define our variable ranges/sequences
#for i in $(seq .01 .02 .09)
for i in $(seq .008 .001 .015)
do
  for j in $(seq .1 .05 .3)
  do
    for k in .089
    #for k in $(seq .085 .001  .09)
    do
      # replace sActsS with i, testingS with j, and adhS with k
      # creates a copy of the file so that you don't replace those variable names in the original, allowing you to change them with each subsequent loop
      sed "s/MATcalChange/$i/g" trt.py > calTrt.py
      sed -i "s/discScalar/$j/g" calTrt.py
      # run subTitan with your flags
      #sed -i "s/probCal/0.1/g" calTrt.py
      #sed -i "s/oatIDU/$k/g" calTrt.py
      ./subTitan.sh calTrt.py -j trt_$i\_$j -N 55000 -s 0 -t 260 -n 5 -b 104 -m 3
    done
  done
done
 
