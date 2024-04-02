#!/bin/bash
# rm TrigParForScan_Pileup_15nA/*
AllChannelsOfInterest=(11 12 14 15 16 17 18)

# ScanChannels=(62 64 78 84 84 74 76)
# ScanRuns=(500904 500909 500911 500913 500915 500918)
ScanChannels=(11 12 14 15 16 17 18)
ScanRuns=(500904 500911 500915 500918) # 500939 500941 500942 500943)
declare -A mintmap
mintmap[500904]=8000
mintmap[500911]=34000
mintmap[500915]=8000
mintmap[500918]=50000
declare -A maxtmap
maxtmap[500904]=16000
maxtmap[500911]=46000
maxtmap[500915]=20000
maxtmap[500918]=80000

AveragesScan=(4) #(2 4) #
DebounceScan=(4) #(2 4) # 

#check directories
Path_Base='/data/users/berettam/DianaHallC/ThresholdScan_NoiseStudies'; [ ! -d "$Path_Base" ] && mkdir $Path_Base && echo "$Path_Base directory created."
TrigParDir='/data/users/berettam/DianaHallC/ThresholdScan_NoiseStudies/TrigScanNoise_Parameters'; [ ! -d "$TrigParDir" ] && mkdir $TrigParDir && echo "$TrigParDir directory created."

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
        /data/users/berettam/DianaHallC/TriggerScan/LaunchTriggerScan.sh -r ${ScanRuns[i]} -p $DebDir -c ${ScanChannels[j]} -f $trigparfile -s 0.05 -t ${mintmap[${ScanRuns[i]}]}  -T ${maxtmap[${ScanRuns[i]}]} -e 0.005 -w 10
      done
    done
  done
  

done

