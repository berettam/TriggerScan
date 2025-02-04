#!/bin/bash
# rm TrigParForScan_Pileup_15nA/*
AllChannelsOfInterest=(1 2 3 4 5 6 7 8 9 11)

ScanChannels=(1 2 3 4 5 6 7 8 9 11)
ScanRuns=(800019) # 500939 500941 500942 500943)
declare -A mintmap
mintmap[800019]=600
declare -A maxtmap
maxtmap[800019]=20000

AveragesScan=(2 4 6 8 10) #(2 4) #
DebounceScan=(2 4 6 8 10) #(2 4) # 

#check directories
BASEDIR='/user/gr2/freddo/berettam/Documents/'
Path_Base=$BASEDIR'/TrigScanOutput'; [ ! -d "$Path_Base" ] && mkdir $Path_Base && echo "$Path_Base directory created."
TrigParDir=$BASEDIR'/TrigScanOutput/TrigScanNoise_Parameters'; [ ! -d "$TrigParDir" ] && mkdir $TrigParDir && echo "$TrigParDir directory created."

 
for Average in ${AveragesScan[@]}; do
  echo "Average = $Average"
  AvgDir=$Path_Base/Avg${Average}; [ ! -d "$AvgDir" ] && mkdir $AvgDir && echo "$AvgDir directory created."
  for Debounce in ${DebounceScan[@]}; do
    echo "  Debounce = $Debounce"
    DebDir=$AvgDir/Deb$Debounce; [ ! -d "$DebDir" ] && mkdir $DebDir && echo "  $DebDir directory created."
    #create and fill the file of parameters
    trigparfile=$TrigParDir/TrigParamsScan_Avg${Average}_Debounce${Debounce}.txt 
    touch $trigparfile
    for t in ${AllChannelsOfInterest[@]}; do
      echo $t $Average 2 $Debounce >> $trigparfile
    done

    for i in "${!ScanRuns[@]}"; do
      for j in "${!ScanChannels[@]}"; do
        echo ${ScanRuns[i]} ${ScanChannels[j]} $Average $debounce ${mintmap[${ScanRuns[i]}]} ${maxtmap[${ScanRuns[i]}]}
        $PWD/LaunchTriggerScan.sh -r ${ScanRuns[i]} -p $DebDir -c ${ScanChannels[j]} -f $trigparfile -s 0.1 -t ${mintmap[${ScanRuns[i]}]}  -T ${maxtmap[${ScanRuns[i]}]} -e 0.005 -w 1
      done
    done
  done
  

done

