#initialize this as a bash script
#!/bin/bash
#define our variable ranges/sequences
for i in .3 #.008 .01 .05 .1 .5
#for i in $(seq .08 .01 .3) 
do
         # replace sActsS with i, testingS with j, and adhS with k
         # creates a copy of the file so that you don't replace those variable names in the original, allowing you to change them with each subsequent loop
    sed "s/HFscalar/$i/g" incarcerate.py > RI_OD_sanity.py
         # run subTitan with your flags
    ./subTitan.sh RI_OD_sanity.py -j sanity_060719 -N 55000 -s 0 -t 260 -n 100 -b 104 -r 1 -m 3
done

 
