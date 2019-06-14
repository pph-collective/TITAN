#initialize this as a bash script
#!/bin/bash
#define our variable ranges/sequences
for i in $(seq .6 .2 2) 
do
         # replace sActsS with i, testingS with j, and adhS with k
         # creates a copy of the file so that you don't replace those variable names in the original, allowing you to change them with each subsequent loop
    sed "s/mortScalar/$i/g" calDeath.py > calDeath_copy.py
         # run subTitan with your flags
    ./subTitan.sh calDeath_copy.py -j death_$i -N 55000 -s 0 -t 260 -n 5 -b 104 -r 1 -m 3
done

 
