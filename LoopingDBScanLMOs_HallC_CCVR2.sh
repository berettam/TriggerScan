#!/bin/bash
# rm TrigParForScan_Pileup_15nA/*
AllChannelsOfInterest=(1 2 4 5 6)

# ScanChannels=(62 64 78 84 84 74 76)
# ScanRuns=(500904 500909 500911 500913 500915 500918)
# ScanChannels=(1 2 4 5 6)
# ScanRuns=(500936 500939) # 500939 500941 500942 500943)
# AveragesScan=(5) #(2 4) #
# DebounceScan=(5) #(2 4) # 


ScanChannels=(1 2 6)
ScanRuns=(500915 500915) # 500939 500941 500942 500943)
AveragesScan=(6 7) #(2 4) #
DebounceScan=(6 7) #(2 4) # 

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
        echo ${ScanRuns[i]} ${ScanChannels[j]} $Average $debounce
        /data/users/berettam/DianaHallC/TriggerScan/LaunchTriggerScan.sh -r ${ScanRuns[i]} -p $DebDir -c ${ScanChannels[j]} -f $trigparfile -s 0.03 -t 100 -T 1100 -e 0.005
      done
    done
  done
  

done

