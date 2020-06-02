for i in $(seq 0.01 0.01 0.15)
do # do the following each time
      # sed finds and replaces any string in a file. Make sure that the string you replace is not repeated where it shouldn't be replaced!
      # the first uses your blank file with keywords and creates a second .py file
      sed "s/dx_black/$i/g" demographics.yml > for_calibration/Atlanta/demographics.yml
      sed "s/dx_white/$i/g" -i for_calibration/Atlanta/demographics.yml
      ./subTitan.sh ~/data/sbessey/TITAN/for_calibration/Atlanta -j cal_black$i\white$i -n 20 -r 5 -f atl_dx_0407
      sleep 2m
done # closes outer loop. one "done" per "do"
