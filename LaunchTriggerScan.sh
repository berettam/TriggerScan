#!/bin/bash

############################################################
# Help                                                     #
############################################################
Help()
{
	# Display Help
	echo "Script to perform the threshold scan on multiple channels/runs."
	echo "Calls the code TriggerScan.cpp (compiled to TrigScan) recursively on a specified run/channel range"
	echo
	echo "Syntax: ./WriteAverages.sh [-h|r|R|p|c|e|w|s]"
	echo "options:"
	echo "r	The list of runs (COMMA SEPARATED) to be used. "
	echo "R	Last run. If provided, check for all the runs between first -r and -R"
	echo "p	path where to save the outputs. In this folder the scripts creates run<runnumber> directories."
	echo "c	COMMA SEPARATED list of channels."
	echo "e	minimum rate to be analyzed [0.001]"
	echo "w	window to be fetched from the RDCF [10 s]"
	echo "s	step of the threshold scan [0.025]"
	echo "t	start time - when I start looking at the run in sec [100]"
	echo "T	stop time - when I stop looking at the run in sec [2000]"
	echo "f	Trigger parameter file [ default is empty]"
	echo "	- space separated txt"
	echo "	- ch avg th debounce"
	echo
	echo "-r, -p , and -c are mandatory!!!!"
	echo
	echo "Example: "
	echo ".LaunchTriggerScan.sh -r 750055 -c 1 -p /raid1-pcsingle/users/shared_data_meno3/PROCESSED/macros/triggerscanoutput/"
}



#checking if there are all the arguments
if [ $# -eq 0 ]
  then
    echo "No arguments supplied, printing help."
	Help
	exit 1
	
fi

echo "###############################"
echo "SCRIPT TO SCAN TRIGGER THs"
echo "###############################"

#check the args
while getopts "hr:R:p:c:e:w:s:t:T:f:" option; do
	case $option in
		h) #display help
			Help
			exit;;
		r) #first or only run to be converted
			set -f
			IFS=',' 
			runList=(${OPTARG});;
		R) #Last run to be converted
			LASTRUN=${OPTARG};;
		p) #path to save the output
			outdir=${OPTARG};;
		c) #list of channels
			set -f
			IFS=',' 
			chArray=(${OPTARG});;
		e) #minumum rate
			minrate=${OPTARG};;
		w) #window length
			winlength=${OPTARG};;
		s) #th step
			thstep=${OPTARG};;
		t) #min time
			tmin=${OPTARG};;
		T) #max time
			tmax=${OPTARG};;
		f) #trigger parameters file
			trigparfile=${OPTARG};;
		?/) #default
			echo "Invalid option!, printing help"
			echo
			Help
			exit;;
	esac
done

#checking mandatory args
if [ -z "$runList" ] ; then
	echo 'Missing -r, first run or run list' >&2
	echo
	echo 'Printing Help'
	Help
	exit 1
fi
if [ -z "$outdir" ]; then
	echo 'Missing -p, output parent directory' >&2
	echo
	echo 'Printing Help'
	Help
	exit 1
fi
if [ -z "$chArray" ] ; then
	echo 'Missing -c, Channel list' >&2
	echo
	echo 'Printing Help'
	Help
	exit 1
fi

#Defaults for non-mandatory arguments
if [ -z "${minrate}" ] ; then
		echo
		minrate=0.001
        echo 'Default minimum rate: '${minrate}
		echo
fi

if [ -z "${winlenght}" ] ; then
		echo
        winlength=10;
		echo 'Default window lengh: '${winlength}
		echo
fi

if [ -z "${thstep}" ] ; then
		echo
        thstep=0.025;
		echo 'Default threshold step: '${thstep}
		echo
fi

if [ -z "${tmin}" ] ; then
		echo
        tmin=100;
		echo 'Default tmin: '${tmin}
		echo
fi

if [ -z "${tmax}" ] ; then
		echo
        tmax=2000;
		echo 'Default tmax: '${tmax}
		echo
fi

if [ -z "${trigparfile}" ]; then
	echo
	echo "No trigger file provided, getting params from DB"
	trigparfile=""
fi
#composing the array list in case of -R option
echo
if [ -z "$LASTRUN" ]
then 
	#vopy the list
	runArray=(${runList[@]})
else
	#all the run, printing the list
	runArray=(${runList[0]})
	for (( rr=$((${runList[0]}+1)); rr<=${LASTRUN};++rr )) #rr in $(seq $((${runList[0]}+1)) 1 ${LASTRUN})
	do
		runArray+=($rr) 
	done
fi


#################################################
#Main body of the script
#################################################

echo
echo "Scanning the runs: "
for run in ${runArray[@]}
do
	echo $run
done
echo
echo "outputs in " ${outdir}
echo
echo "Scanning channels: "
for ch in ${chArray[@]}
do
	echo $ch
done
echo
echo "Fetching "${winlength}" s windows from "${tmin}" to "${tmax}
echo

#run cycle
for f in "${runArray[@]}" 
do 
	echo "Run " $f;
	#creating the directory -- if does not exists
	Path_output=${outdir}'/run'${f}'/'
	if ! test -d "${Path_output}";
		then 
			echo "--> Creating the output directory";
			mkdir -p ${Path_output}
	fi
	#cycle over channels	
	for i in "${!chArray[@]}"
		do
			echo "--> CH " ${chArray[$i]} 
			echo
			./TrigScan -r $f -c ${chArray[$i]} -t${tmin} -T${tmax} -w${winlength} -o ${Path_output} -m${minrate} -s${thstep} -f${trigparfile} -I0.000001 -A10
	done
done
