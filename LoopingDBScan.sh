#!/bin/bash
rm TrigParForScan/*

AllChannelsOfInterest=(2 4 6 8 10 12 14 16 18 20)

ScanChannels=(12 14 16 18 20)
AveragesScan=(2 6 10 20 40 60 80)
DebounceScan=(2 4 6 8 10)

#check directories
Path_Base='TrigScan_MultiSpace'; [ ! -d "$Path_Base" ] && mkdir $Path_Base && echo "$Path_Base directory created."
TrigParDir='TrigScan_Parameters'; [ ! -d "$TrigParDir" ] && mkdir $TrigParDir && echo "$TrigParDir directory created."

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

    for ch in ${ScanChannels[@]}; do
      ./LaunchTriggerScan.sh -r 950035 -p $DebDir -c $ch -f $trigparfile -s 0.01
    done

  done
  

done

