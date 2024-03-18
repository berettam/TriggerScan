/*

g++ -o TrigScan -std=c++0x -D_HAVE_FFTW3_ `root-config --cflags --glibs` -lRooFit `diana-config --cflags --libs` TriggerScan.cpp

*/
#include <cstdlib>
#include <complex>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <getopt.h>

#include <TPaveText.h>
#include <TString.h>
#include <TMath.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TRandom.h>
#include <TMarker.h>
#include <TH2F.h>
#include <TSystem.h>
#include <TStyle.h>

#include "QCuore.hh"
#include "QComplex.hh"
#include "QVector.hh"
#include "QVectorC.hh"
#include "QMatrix.hh"
#include "QMatrixC.hh"
#include "QRdcfRootFileReader.hh"
#include "QRdcfRootFileReaderHandler.hh"
#include "QRealComplexFFTW3.hh"
#include "QFFT.hh"
#include "QGlobalHandle.hh"
#include "QChain.hh"
#include "QRawWaveformGetter.hh"
#include "QDetTrgParamsDerivativeHandle.hh"
#include "QDetTrgParamsDerivative.hh"
#include "QGlobalDataManager.hh"
#include "QDetChannelHandle.hh"

using namespace Cuore;
//apollo like trigger
int SearchForTrigger(QVector timestream,double Threshold, double AvgWin, double Debounce, double adc2mV) {

	// threshold in derParams is in mV per ms
  // threshold is Multiplied by AvgWin to speed up calculation
  // Getting the SAME th normalization as in QApolloDerivativeTrigger.cc
  Threshold *= (AvgWin/adc2mV); //Threshold avg slope for determining a pulse (in units of ADC * Sampling Frequency)
	Debounce *= 2; 	
	AvgWin *= 2; 	
	
  // Number of samples above threshold
  int fCount = 0;

  // Loop from start to end in derivative
  unsigned int bLow, bUp;
  int k=0;
	int TrigNumber = 0;
	int DeadTime = 100*2;
	while(k < (int) timestream.Size()-AvgWin-DeadTime){
    bLow = k;
    bUp = k + AvgWin;

    if(bLow >= timestream.Size() || bUp >= timestream.Size()) return 0;
    // Evaluate derivative jump in a wider window (width=AvgWin)
    int low = timestream[bLow];
    int up = timestream[bUp];
    int delta = up - low;

    // If delta exceeds threshold, check for triggers
    if (delta > Threshold) {
      fCount++;
      // If above threshold for enough points (Debounce), fill trigger
      if (fCount == Debounce) {
        TrigNumber++;
        fCount = 0;
        k += DeadTime;
      }
    } else {
      // Not above threshold, reset count
      fCount = 0;
    }
		k++;
  }
  return TrigNumber;
}

void GetSamples(const QRdcfRootFileReader* reader, const int ch, const Long64_t &startns, const Long64_t &endns, QVector& wave){
    QError err = QERR_SUCCESS;
    err = reader->GetSamples(wave, ch, startns, endns);

    if (err == QERR_OUT_OF_RANGE){
        return;
    }
    return;
}

void GetTimestream(const QRdcfRootFileReader* reader, const int &ch, const double &tstart, const double &tend, QVector& wave){
    Long64_t startns = (Long64_t) (std::round(tstart*1.0e9));
    Long64_t endns = (Long64_t) (std::round(tend*1.0e9));

	QVector temp;
    GetSamples(reader, ch, startns, endns, temp);
	wave = temp;

}


void parseListToValues(const std::string& listOfValues, std::vector<std::string>& parsedList) {
    if (listOfValues != "") {
        std::stringstream ss(listOfValues);
        while (ss.good()) {
            std::string substr;
            getline(ss, substr, ' ');
            parsedList.push_back(substr);
        }
    }
}

double CheckForTriggers(std::map<double,QVector> VecTCH,double threshmV, double average, double debounce, double adc2mV){
	double TotalWindows = 0.;
	double TrigWindows = 0.;

	for(std::map<double,QVector>::iterator iter = VecTCH.begin(); iter != VecTCH.end(); ++iter)
	{
		QVector ts = iter->second;
	
		TotalWindows++;
		
		//count how many pulses there are 	
		int TriggerNumber = SearchForTrigger(ts,threshmV,average,debounce, adc2mV);
		if(TriggerNumber>0) TrigWindows+=TriggerNumber;
	}//cycle over times is over --> from now working on the avgs

	return TrigWindows;///TotalWindows*100.;
}


void Usage()
{

  std::cout << std::endl << std::endl;
  std::cout << "---------------------------------------------------------------------------------" << std::endl;
  if( std::getenv("USER") != NULL )
    std::cout << "Hi " << std::getenv("USER") <<"! The usage of this wonderful program is: " << std::endl;
    else
      std::cout << "Hi there! The usage of this wonderful program is:" << std::endl;
  std::cout<<"./TrigScan"<<std::endl;
  std::cout<<"Run with the following options "<<std::endl;
  std::cout<<"options "<<std::endl;
  std::cout<<"------------------------------------------------------------"<<std::endl;
  std::cout<<"-r (or --run)             run to be processed                                   [REQUIRED] "<<std::endl;
  std::cout<<"-c (or --cahnnel)         channel to be analyzed                                [REQUIRED] "<<std::endl;
  std::cout<<"-t (or --tstart)          start time for scan                                   [100 s]"<<std::endl;
  std::cout<<"-T (or --tstop)           stop time for the scan                                [2000s]"<<std::endl;
  std::cout<<"-w (or --twindow)         window length to do the scan                          [10 s]"<<std::endl;
  std::cout<<"-o (or --outputpath)      path where to save the output                         [REQUIRED]"<<std::endl;
  std::cout<<"-m (or --minrate)         rate at wich to stop the scan                         [0.001 Hz]"<<std::endl;
  std::cout<<"-s (or --step)            step in threshold used for the scan                   [0.025]"<<std::endl;
  std::cout<<"-I (or --minthreshold)    threshold value as start for the scan                 [0.001]"<<std::endl; 
  std::cout<<"-A (or --maxthreshold)    threshold value as stop for the scan                  [1e8]"<<std::endl; 
  std::cout<<"-f (or --trigparamsfile)  file for the trigger scans. If not provided, using DB [empty]"<<std::endl; 
  std::cout<<"-h this help"<<std::endl;
}

/*
 *
 * MAIN
 *
 */

#ifndef __CINT__
int main(int argc, char** argv){
  /*
  int run = atoi(argv[1]);
  int ch = atoi(argv[2]);
  int tstart = atoi(argv[3]);
  int tend = atoi(argv[4]);
  double twindow = atof(argv[5]);
	double MinRate = atof(argv[7]);
	double thstep_set = atof(argv[8]);
	std::cout<<"Min Rate = "<<MinRate<<std::endl;
  TString output = Form("%s/run%d", argv[6],run);
  //check if there is the trigparams file
  //if not, the trig parames will be read from the DB
  bool InputTrigFile = false; 
	if(argc==10)	InputTrigFile = true;

  if(argc<9 ){
    std::cout<< "Scan trigger parameters for a specific channel" << std::endl;
    std::cout<< "Run as "<< argv[0] << " <run> <ch> <tstart> <tend> <twindow> <output_path> <min rate> <th step> <triggerparamsfile>" << std::endl;
    std::cout<<std::endl;
    std::cout<< "The triggerparamsfile assumes this format for each line, without header: "<<std::endl;
    std::cout<< "ch average threshold debounce"<<std::endl;
    std::cout<<std::endl;
    return 1;
  }
*/
  if(argc<2 ){
    std::cout<< "No arguments provided!!!" << std::endl;
    Usage();
    return 1;
  }

  int run = 0;
  int ch = 0;
  int tstart = 0;
  int tend = 0;
  double twindow = 0;
	double MinRate = 0;
	double thstep_set = 0;
  double thmin = 0;
  double thmax = 0;
  TString output = "";//Form("%s/run%d", argv[6],run);
  //check if there is the trigparams file
  //if not, the trig parames will be read from the DB
  bool InputTrigFile = false; 
	//if(argc==10)	InputTrigFile = true;
  std::string trigfile = "";

  
  // Parse the options                                                                                                                                                                                                                                                                                                                              
  
  static struct option long_options[] = {
                                        { "run",        required_argument, 0, 'r' },
                                        { "channel",required_argument,0,'c'},
                                        { "tstart",optional_argument,0,'t'},
                                        { "tstop",optional_argument,0,'T'},
                                        { "twindow",optional_argument,0,'w'},
                                        { "outpath",required_argument,0,'o'},
                                        { "minrate", optional_argument, 0, 'm' },
                                        { "step",optional_argument,0,'s'},
                                        { "minthreshold",optional_argument,0,'I'},
                                        { "maxthreshold",optional_argument,0,'A'},
                                        { "trigparamsfile",optional_argument,0,'f'},
                                        { "help",no_argument, 0, 'h' },
                                        {0, 0, 0, 0}
  };
  const char* const short_options = "r:c:t::T::w::o:m::s::I::A::f::h";
  int c;

  while ((c = getopt_long(argc, argv, short_options, long_options, nullptr)) != -1 ) {
    
    switch (c) {
      case 'r': {
        run = atoi(optarg);
        break;
      }
      case 'c':{
        ch=atoi(optarg);
        break;
      }
      case 'o':{
        output=optarg;
        std::cout<<"\n\n\noutdir = "<<output<<"\n\n\n";
        break;
      }
      case 'h':{
        Usage();
        break;
      }
      case 't':{
        tstart= (optarg==NULL) ? 100 : atof(optarg);
        break;
      }
      case 'T':{
        tend=(optarg==NULL) ? 2000: atof(optarg);
        break;
      }
      case 'w':{
        std::cout<<"QUAAAAAAA"<<std::endl;
        twindow= (optarg==NULL) ? 10 : atof(optarg);
        break;
      }
      case 'm':{
        MinRate=(optarg==NULL) ? 0.001 : atoi(optarg);
        break;
      }
      case 's':{
        thstep_set=(optarg==NULL) ? 0.025 : atof(optarg);
        break;
      }
      
      case 'I':{
        thmin=(optarg==NULL) ? 0.0001 : atof(optarg);
        break;
      }
      case 'A':{
        thmax= (optarg==NULL) ? 1e8: atof(optarg);
        break;
      }
      case 'f':{
        if(optarg == NULL){
          InputTrigFile = false;
        }
        else{
          std::cout<<"QUA!!"<<std::endl;
          InputTrigFile = true;
          trigfile = optarg;
        }
        break;
      }    
      default: {
        std::cout<<"Unknown argument!!!!"<<std::endl;
        exit(1);
      }
    }
  }
    


  
	


    /********************************THESE MAY NEED TO BE RE-DEFINED********************************/
	// std::vector<int> auxchannels = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
	std::vector<int> auxchannels = {62,64,74,76,78,84};
	/********************************************************************************************************************/
    
	double adc2mV = 8.0108642578125e-02/64.0;
  double fSampleFreq = 2000;//berettam: putting a fixed value, should change according to channel FIXME
	
	//READING the DAQ parameters form DBi
	
	std::cout<<"Getting DAQ parameters from the DB"<<std::endl;
	QGlobalDataManager dmdaq;
	//std::map<int,double> adc2mV_byCH;
	//std::map<int,double> SF_byCH;
	QDetChannelHandle handle(ch);
	handle.SetRun(run);
	dmdaq.Get(&handle,"DB");
	QDetChannel detChannel = handle.Get();
	adc2mV = 1000.*( detChannel.fDaqSet.fVAdcLimitMax - detChannel.fDaqSet.fVAdcLimitMin ) / pow(2,(double)( detChannel.fDig.fDaqNBits ));
	fSampleFreq = (double) detChannel.fDaqSet.fSamplingFrequency; 
	
	std::cout<<"SF = "<<fSampleFreq<<" -- adc2mV = "<< adc2mV<<std::endl;

	//declaring the reading structure
  QRdcfRootFileReaderHandler& handler  = QRdcfRootFileReaderHandler::GetInstance();
  QError err = QERR_SUCCESS;
	//map to keep the readers
	std::map< int , const QRdcfRootFileReader* > readersCH;
	readersCH[ch] = handler.GetReader(run, ch, err);
	if(!readersCH[ch] || err != QERR_SUCCESS){
		std::cout << err << std::endl;
		std::cout << "No file found for ch: "<< ch << std::endl;
		return 231192;
	}

    //Create the QVector(C) instances we will need
    QVector ts(0);

    //Get the trigger thresholds for this channel-runs from the database
	std::map< int,double > threshmVCH;
	std::map< int,double > averageCH;
	std::map< int,double > debounceCH;
  if(InputTrigFile){
		//open the file
		std::cout<<"Reading trigger from file: "<<trigfile<<std::endl;
		std::cout<<"assuming this format for each line: "<<std::endl;
		std::cout<<"ch average threshold debounce"<<std::endl;
		std::string line;
		ifstream trigfilename(trigfile.c_str()) ;
		while (getline(trigfilename,line))
		{
			std::vector<std::string> vec;
			//std::cout<<line<<std::endl;
			parseListToValues(line,vec);
			averageCH[stoi(vec[0])] = stof(vec[1]);
			threshmVCH[stoi(vec[0])] = stof(vec[2]);
			debounceCH[stoi(vec[0])] = stof(vec[3]);
		}
		trigfilename.close();
	}
	else{
		//trigger params from the DB	
		std::cout<<"Getting the trigger parameters from the DB"<<std::endl;
		QGlobalDataManager dm;
		QDetTrgParamsDerivative *result = NULL;
		
		for(std::vector<int>::iterator it = auxchannels.begin(); it!=auxchannels.end(); ++it )
		{
				QDetChannelHandle handle(*it);
				handle.SetRun(run);
				dm.Get(&handle,"DB");
				QDetChannel detChannel = handle.Get();
				result = (QDetTrgParamsDerivative *)(detChannel.fTrgSet.fTriggers[0].GetTrgParams());
				threshmVCH[*it] =result->fThreshold;
				averageCH[*it] =result->fAverage;
				debounceCH[*it] =result->fDebounce;
		}	
	}	
	//printing the channels parameters
	std::cout<<"Using the following parameters for triggering:"<<std::endl;
	std::cout<<"CH\tAvg\tTH\tDeb"<<std::endl;
	//if the params were read from the DB, write them in a txt file
	// so it is easier to customize them in case
	std::ofstream outtrigp ;
	if(!InputTrigFile){
		outtrigp.open(Form("TrigParams_Run%d.txt",run),std::ios::out);
	}
	for(size_t i =0; i<auxchannels.size(); i++){
		std::cout<<auxchannels[i]<<" "<<averageCH[auxchannels[i]]<<" "<<threshmVCH[auxchannels[i]]<<" "<<debounceCH[auxchannels[i]]<<std::endl;
		if(!InputTrigFile){
			outtrigp<<auxchannels[i]<<" "<<averageCH[auxchannels[i]]<<" "<<threshmVCH[auxchannels[i]]<<" "<<debounceCH[auxchannels[i]]<<std::endl;
		}
	}
	if(!InputTrigFile) outtrigp.close();
	std::cout<<"--> Read trigger Done"<<std::endl<<std::endl;


	// take the number of points to predeclare the vector
    int N = twindow*fSampleFreq; 


    //finally, some counters we need
    double t0 = tstart;
    

	//loop over the noise events
	//for the channel of interest
	//saving the ok windows w and their t0 in a map
	std::cout<<"Looping over time for ch "<<ch<<std::endl;
	std::cout<<"twindow = "<<twindow<<std::endl;
	std::map< double, QVector> VecTCH;
    while(t0 < tend){ 
        //std::cout << "Reached time " << t0 << " s" << std::endl;

        //Get timestream from the channel of interest
        GetTimestream(readersCH[ch], ch,t0,t0+twindow, ts);
		
		if((int) ts.Size() != N){
            std::cout << "Size did not match at t0 = " << t0 << std::endl;
            std::cout << "Exiting the program " << std::endl;
            break; //EXIT ONCE THE TIMESTREAM IS NOT SUFFICIENTLY LONG. DO NOT DELETE
        }
		//saving in the map the event
		VecTCH[t0] = ts;
        //increasing the time
		t0+= twindow;
	}

    std::cout <<std::endl<< "Getting " << VecTCH.size() <<" windows " << std::endl;
	std::cout<<"-->Done Reading"<<std::endl<<std::endl;

	//now doing the trigger analisys

	//HERE THE CYCLE TO SCAN THE TRIGGER TH
	std::cout<<"Scan the trigger Threshold"<<ch<<std::endl;

	//non scanned parameters	
	double avg = averageCH[ch];
	double deb = debounceCH[ch];	

	//limits	
	//variables
	double rate;
	double TriggeredWindows =0;	
	double tottime = tend-tstart;
	//saving the steps
	std::vector<double> thresholds;
	std::vector<double> rates;
	//starting parameters
	rate = 10000;
	double thstep = thstep_set;//0.025;
	//memory of the previous iteration
	double PrevTrigWindows =  0;

	//SCAN ROUTINE
	for(double thscan = thmin; thscan<thmax && rate>=MinRate; thscan+=thstep)
	{
		//get the trigger over the windows
		TriggeredWindows = CheckForTriggers(VecTCH,thscan,avg,deb,adc2mV);
		//calculate the rate at this threshold
		rate = TriggeredWindows/tottime;
		//save the current point
		thresholds.push_back(thscan);
		rates.push_back(rate);
		std::cout<<Form("Th = %f - TWindow = %.0f - rate = %.3f",thscan,TriggeredWindows,rate)<<std::endl;
		//If the number of triggered samples is equal to the previous step
		//move double, to avoid getting stuck
		if(TriggeredWindows == PrevTrigWindows) 
		{
			thstep*=2;
		}
		else{
			thstep = thstep_set;
		}
		//update the memory
		PrevTrigWindows = TriggeredWindows;
		
	}
	
	std::cout<<"-->Done Scanning"<<std::endl<<std::endl;

	//define the output file	
	TFile* newFile = new TFile(Form("%s/ScanTHs_run%d_CH%d.root",output.Data(),run,ch),"RECREATE");
	TGraph *gWave = new TGraph(thresholds.size());
    TCanvas *can_ts = new TCanvas(Form("RatevsThreshold_run%d_CH%d",run,ch), Form("RateVsThreshold_run%d_CH%d",run,ch), 1000, 750);
	
	std::ofstream newFileTXT(Form("%s/ScanTHs_run%d_CH%d.txt",output.Data(),run,ch));
  
	newFileTXT <<  "Threshold,Rate"<<std::endl;
	
	for(size_t kk = 0; kk < thresholds.size(); ++kk)
	{
		gWave->SetPoint(kk,thresholds[kk],rates[kk]);
		newFileTXT <<  thresholds[kk]<<","<<rates[kk]<<std::endl;
	}
	newFileTXT.close();
    gWave->SetTitle(Form("Th Scan for CH %d; threshold; Rate [Hz] ",ch));
    gWave->SetLineColor(kBlack);
    gWave->SetLineWidth(2);
    gWave->SetMarkerStyle(20);
    gWave->Draw("AP");

	//Drawing a line with the chosen threshold from SB
	TLine line (threshmVCH[ch],MinRate,threshmVCH[ch],*max_element(std::begin(rates), std::end(rates)));
	line.SetLineColor(kRed);
	line.SetLineWidth(2);
	line.Draw("SAME");


    gPad->SetGridx();
    gPad->SetGridy();
    can_ts->Update();
    can_ts->Write("can_ts");
	newFile->Close();
    std::cout << "Parameters saved" << std::endl;
  return 0;
}
#endif
