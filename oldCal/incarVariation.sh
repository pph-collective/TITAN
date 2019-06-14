#initialize this as a bash script
#!/bin/bash
#define our variable ranges/sequences
for i in $(seq .003 .0001 .005) 
do
         # replace sActsS with i, testingS with j, and adhS with k
         # creates a copy of the file so that you don't replace those variable names in the original, allowing you to change them with each subsequent loop
    sed "s/incarParam/$i/g" incarceration.py > calIncar.py
         # run subTitan with your flags
    ./subTitan.sh calIncar.py -j incarceration20191705_$i -N 55000 -s 0 -t 260 -n 5 -b 104 -r 20 -m 3
done

 
