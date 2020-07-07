#!/bin/bash
for i in 0.00 0.15 0.90 1.00 #0.30 0.45 0.60 0.75 0.90 1.00 #$(seq [start] [step] [stop]) # second loop
  do # second loop, can repeat as needed depending on number of variables
      # sed finds and replaces any string in a file. Make sure that the string you replace is not repeated where it shouldn't be replaced!
      # the first uses your blank file with keywords and creates a second .py file
      sed "s/PrEPlevelsubstitute/$i/g" ~/data/TITAN/TITAN/settings/Atlanta/Atlanta_baseCase_Aly.py > ~/data/TITAN/TITAN/settings/Atlanta/Atlanta_baseCase_Aly_PrEPscenarios_A2.py  # [hold.py] > [submit.py]
      #sed -i "s/[keyword]/$j/g" [submit.py] # -i replaces the keyword in-place. use this for any subsequent loops (k,l,m,etc)
      ~/data/TITAN/TITAN/subTitan.sh ~/data/TITAN/TITAN/settings/Atlanta/Atlanta_baseCase_Aly_PrEPscenarios_A2.py -r 20 -N 17440 -n 5 -j $i\PrEP_BasART_AdhSensScenario_TESTphily #[pathToSubtitan] [submit.py] [Titan arguments]
  done # closes inner loop

 # closes outer loop. one "done" per "do"
 
