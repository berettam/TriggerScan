#!/bin/bash
# rm TrigParForScan_Pileup_15nA/*
AllChannelsOfInterest=(62 64 74 76 78 84)

# ScanChannels=(62 64 78 84 84 74 76)
# ScanRuns=(900964 900965 900966 900967 900968 900968)
ScanChannels=(76)
ScanRuns=(900968)
AveragesScan=(2 4 6 8 10 20)
DebounceScan=(2 4 6 8 10)

#check directories
Path_Base='TrigScanPileup_15nA_MultiSpace'; [ ! -d "$Path_Base" ] && mkdir $Path_Base && echo "$Path_Base directory created."
TrigParDir='TrigScanPileup_15nA_Parameters'; [ ! -d "$TrigParDir" ] && mkdir $TrigParDir && echo "$TrigParDir directory created."

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

    for i in "${!ScanChannels[@]}"; do
      ./LaunchTriggerScan.sh -r ${ScanRuns[i]} -p $DebDir -c ${ScanChannels[i]} -f $trigparfile -s 0.01 -t 3000 -T 8000 -e 0.01 
    done
  done
  

done

