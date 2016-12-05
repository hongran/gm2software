/********************************************************************\

Name:         TestMagnetFe.cxx
Created by:   Matteo Bartolini
Modified by:  Joe Grange, Ran Hong

Contents:     readout code to talk to Galil motion control

$Id$

\********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "midas.h"
#include "mcstd.h"
#include "gclib.h"
#include "gclibo.h"
#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <unistd.h>
#include <sys/timeb.h>
#include <fstream>
#include <sstream>
#include <thread>
#include <mutex>

#define GALIL_EXAMPLE_OK G_NO_ERROR //return code for correct code execution
#define GALIL_EXAMPLE_ERROR -100
using namespace std;

/* make frontend functions callable from the C framework */
#ifdef __cplusplus
extern "C" {
#endif

  // i am defining some Galil libraries variables

  //----------------------------------------------------------
  /*-- Globals -------------------------------------------------------*/

  /* The frontend name (client name) as seen by other MIDAS clients   */
  char *frontend_name = "TestMagnetGalilFe";
  /* The frontend file name, don't change it */
  char *frontend_file_name = __FILE__;

  /* frontend_loop is called periodically if this variable is TRUE    */
  BOOL frontend_call_loop = FALSE;

  /* a frontend status page is displayed with this frequency in ms */
  INT display_period = 3000;

  /* maximum event size produced by this frontend */
  INT max_event_size = 100000;

  /* maximum event size for fragmented events (EQ_FRAGMENTED) */
  INT max_event_size_frag = 5 * 1024 * 1024;

  /* buffer size to hold events */
  INT event_buffer_size = 100 * 10000;


  /*-- Function declarations -----------------------------------------*/
  INT frontend_init();
  INT frontend_exit();
  INT begin_of_run(INT run_number, char *error);
  INT end_of_run(INT run_number, char *error);
  INT pause_run(INT run_number, char *error);
  INT resume_run(INT run_number, char *error);
  INT frontend_loop();

  INT read_event(char *pevent, INT off);

  INT poll_event(INT source, INT count, BOOL test);
  INT interrupt_configure(INT cmd, INT source, POINTER_T adr);


  /*-- Equipment list ------------------------------------------------*/


  EQUIPMENT equipment[] = {

    {"Galil",                /* equipment name */
      {1, 0,                   /* event ID, trigger mask */
	"SYSTEM",               /* event buffer */
	EQ_POLLED,            /* equipment type */
	0,                      /* event source */
	"MIDAS",                /* format */
	TRUE,                   /* enabled */
	RO_RUNNING   /* read when running*/
	//|RO_ODB
	,
	500,                  /* poll every 0.5 sec */
	0,                      /* stop run after this event limit */
	0,                      /* number of sub events */
	0,                      /* log history, logged once per minute */
	"", "", "",},
      read_event,       /* readout routine */
    },

    {""}
  };


#ifdef __cplusplus
}
#endif

HNDLE hDB, hkeyclient;

typedef struct GalilDataStruct{
  INT TimeStamp;
  INT PositionArray[4];
  INT VelocityArray[4];
  INT OutputVArray[4];
}GalilDataStruct;

typedef struct GalilDataStructD{
  double PositionArray[4];
  double VelocityArray[4];
  double OutputVArray[4];
}GalilDataStructD;

INT StepSize[4];
INT StepNumber[4];
int IX,IY,IZ,IS;

thread monitor_thread;

//Flags
bool ReadyToMove = false;
bool ReadyToRead = false;

mutex mlock;

void GalilMonitor(const GCon &g);

bool MonitorActive;
int ReadGroupSize = 50;

GCon g = 0; //var used to refer to a unique connection. A valid connection is nonzero.

/********************************************************************\
  Callback routines for system transitions

  These routines are called whenever a system transition like start/
  stop of a run occurs. The routines are called on the following
occations:

frontend_init:  When the frontend program is started. This routine
should initialize the hardware.

frontend_exit:  When the frontend program is shut down. Can be used
to releas any locked resources like memory, commu-
nications ports etc.

begin_of_run:   When a new run is started. Clear scalers, open
rungates, etc.

end_of_run:     Called on a request to stop a run. Can send
end-of-run event and close run gates.

pause_run:      When a run is paused. Should disable trigger events.

resume_run:     When a run is resumed. Should enable trigger events.
\********************************************************************/

/*-- Frontend Init -------------------------------------------------*/

INT frontend_init()
{ 
  
  char buf[1023]; //traffic buffer
  GReturn b = G_NO_ERROR;
  b=GOpen("192.168.1.13 -s ALL -t 1000 -d",&g);
  GInfo(g, buf, sizeof(buf)); //grab connection string
//  cout << "connection string is" << " "<<  buf << "\n";
  cm_msg(MINFO,"init","connection string is %s",buf);
  if (b==G_NO_ERROR){
    //cout << "connection successful\n";
    cm_msg(MINFO,"init","Galil connection successful");
  }
  else {
 //   cout << "connection failed \n";
    cm_msg(MERROR,"init","Galil connection failed");
  }

  //Load script
  char ScriptName[500];
  INT ScriptName_size = sizeof(ScriptName);
  db_get_value(hDB,0,"/Equipment/Galil/Settings/CmdScript",ScriptName,&ScriptName_size,TID_STRING,0);
  //string FullScriptName = string("/home/rhong/gm2/g2-field-team/field-daq/resources/GalilMotionScripts/")+string(ScriptName)+string(".dmc");
  string FullScriptName = string("/Users/rhong/gm2software/experiments/TestMagnet/GalilMotionScripts/")+string(ScriptName)+string(".dmc");
//  cout <<"Galil Script to load: " << FullScriptName<<endl;
  cm_msg(MINFO,"init","Galil Script to load: %s",FullScriptName.c_str());

  GProgramDownload(g,"",0); //to erase prevoius programs
  //dump the buffer
  int rc = GALIL_EXAMPLE_OK; //return code
  char buffer[5000000];
  rc = GMessage(g, buffer, sizeof(buffer));
  b=GProgramDownloadFile(g,FullScriptName.c_str(),0);
  GCmd(g, "XQ #TH1,0");
  MonitorActive=true;
  //Start thread
  monitor_thread = thread(GalilMonitor,g);
//  GTimeout(g,2000);//adjust timeout
  //-------------end code to communicate with Galil------------------

  return SUCCESS;
}

/*-- Frontend Exit -------------------------------------------------*/

INT frontend_exit()
{
  mlock.lock();
  GCmd(g,"AB");
  mlock.unlock();
  cm_msg(MINFO,"end_of_run","Motion aborted.");

  mlock.lock();
  GCmd(g,"HX");
  cm_msg(MINFO,"end_of_run","Galil execution stopped.");
  mlock.unlock();

  mlock.lock();
  MonitorActive=false;
  mlock.unlock();
//  cm_msg(MINFO,"end_of_run","Trying to join threads.");
  monitor_thread.join();
  cm_msg(MINFO,"end_of_run","All threads joined.");
  return SUCCESS;
}

/*-- Begin of Run --------------------------------------------------*/

INT begin_of_run(INT run_number, char *error)
{
  //Get run number
  INT RunNumber;
  INT RunNumber_size = sizeof(RunNumber);
  cm_get_experiment_database(&hDB, NULL);
  db_get_value(hDB,0,"/Runinfo/Run number",&RunNumber,&RunNumber_size,TID_INT, 0);

  INT StepSize_size = sizeof(StepSize);
  INT StepNumber_size = sizeof(StepNumber);
  db_get_value(hDB,0,"/Equipment/Galil/AutoControl/RelPos",StepSize,&StepSize_size,TID_INT,0);
  db_get_value(hDB,0,"/Equipment/Galil/AutoControl/StepNumber",StepNumber,&StepNumber_size,TID_INT,0);
 
  IX=IY=IZ=IS=0;

  ReadyToMove = true;
  
  return SUCCESS;
}

/*-- End of Run ----------------------------------------------------*/

INT end_of_run(INT run_number, char *error)

{
  return SUCCESS;
}

/*-- Pause Run -----------------------------------------------------*/

INT pause_run(INT run_number, char *error)
{
//  GCmd(g,"STA");
  return SUCCESS;
}

/*-- Resuem Run ----------------------------------------------------*/

INT resume_run(INT run_number, char *error)
{ 
//  GCmd(g,"SHA");
//  GCmd(g,"BGA"); 
  return SUCCESS;
}

/*-- Frontend Loop -------------------------------------------------*/

INT frontend_loop()
{
  /* if frontend_call_loop is true, this routine gets called when
     the frontend is idle or once between every event */
  return SUCCESS;
}

/*------------------------------------------------------------------*/

/********************************************************************\

  Readout routines for different events

  \********************************************************************/

/*-- Trigger event routines ----------------------------------------*/

INT poll_event(INT source, INT count, BOOL test)
  /* Polling routine for events. Returns TRUE if event
     is available. If test equals TRUE, don't return. The test
     flag is used to time the polling */
{
  static unsigned int i;
  if (test) {
    for (i = 0; i < count; i++) {
      usleep(10);
    }
    return 0;
  }

  if (ReadyToMove)return 1;
  else return 0;
}

/*-- Interrupt configuration ---------------------------------------*/

INT interrupt_configure(INT cmd, INT source, POINTER_T adr)
{
  switch (cmd) {
    case CMD_INTERRUPT_ENABLE:
      break;
    case CMD_INTERRUPT_DISABLE:
      break;
    case CMD_INTERRUPT_ATTACH:
      break;
    case CMD_INTERRUPT_DETACH:
      break;
  }
  return SUCCESS;
}

/*-- Event readout -------------------------------------------------*/

INT read_event(char *pevent, INT off){
  ReadyToMove = false;
  char CmdBuffer[500];
  //Regeister current position from odb
  INT CurrentPositions[4];
  INT CurrentPositions_size = sizeof(CurrentPositions);
  db_get_value(hDB,0,"/Equipment/Galil/Monitors/Positions",CurrentPositions,&CurrentPositions_size,TID_INT,0);
  INT CurrentVelocities[4];
  INT CurrentVelocities_size = sizeof(CurrentPositions);
  db_get_value(hDB,0,"/Equipment/Galil/Monitors/Velocities",CurrentVelocities,&CurrentVelocities_size,TID_INT,0);
  //Take a measurement
  
  //Move platform
  IX++;
  if (IX==StepNumber[0]){
    IX=0;
    //move back to X0
    sprintf(CmdBuffer,"PRA=%d",-(StepNumber[0]-1)*StepSize[0]);
    GCmd(g,CmdBuffer);
    GCmd(g,"BGA");
    IY++;
    if (IY==StepNumber[1]){
      IY=0;
      //move back to Y0
      sprintf(CmdBuffer,"PRB=%d",-(StepNumber[1]-1)*StepSize[1]);
      GCmd(g,CmdBuffer);
      GCmd(g,"BGB");
      IZ++;
      if (IZ==StepNumber[2]){
	IZ=0;
	//move back to Z0
	sprintf(CmdBuffer,"PRC=%d",-(StepNumber[2]-1)*StepSize[2]);
	GCmd(g,CmdBuffer);
	GCmd(g,"BGC");
	IS++;
	if (IS==StepNumber[3]){
	  ReadyToMove = false;
	  return bk_size(pevent);
	}else{
	  //move forward in S
	  sprintf(CmdBuffer,"PRD=%d",StepSize[3]);
	  GCmd(g,CmdBuffer);
	  GCmd(g,"BGD");
	}
      }else{
	//move forward in Z
	sprintf(CmdBuffer,"PRC=%d",StepSize[2]);
	GCmd(g,CmdBuffer);
	GCmd(g,"BGC");
      }
    }else{
      //move forward in Y
      sprintf(CmdBuffer,"PRB=%d",StepSize[1]);
      GCmd(g,CmdBuffer);
      GCmd(g,"BGB");
    }
  }else{
    //move forward in X
    sprintf(CmdBuffer,"PRA=%d",StepSize[0]);
    GCmd(g,CmdBuffer);
    GCmd(g,"BGA");
  }

  //Check if motion is stopped
  while(1){
    db_get_value(hDB,0,"/Equipment/Galil/Monitors/Velocities",CurrentVelocities,&CurrentVelocities_size,TID_INT,0);
    if (CurrentVelocities[0]==0 && CurrentVelocities[1]==0 && CurrentVelocities[2]==0 && CurrentVelocities[3]==0)break;
    sleep(1);
  }
  ReadyToRead = true;
  return bk_size(pevent);
}

//GalilMonitor
void GalilMonitor(const GCon &g){
  int MonitorThreadActive = 1;
  db_set_value(hDB,0,"/Equipment/Galil/Monitors/Active",&MonitorThreadActive,sizeof(MonitorThreadActive), 1 ,TID_BOOL); 
  char buffer[50000];
  hkeyclient=0;
  string Header;
  int  rc = GALIL_EXAMPLE_OK; //return code
  /* residule string */
  string ResidualString = string("");
  INT command = 0;

  int position =0;
//  timeb starttime,currenttime;
//  ftime(&starttime);
 
  //Readout loop
  int i=0;
  int jj=0;
  double Time,Time0;
  while (1){
    bool localMonitorActive;
    mlock.lock();
    localMonitorActive = MonitorActive;
    mlock.unlock();
    if (!localMonitorActive)break;
    //Read Message to buffer
    mlock.lock();
    rc = GMessage(g, buffer, sizeof(buffer));
    mlock.unlock();
//    ftime(&currenttime);
//    double time = (currenttime.time-starttime.time)*1000 + (currenttime.millitm - starttime.millitm);
//    cout<<buffer<<endl;

    string BufString = string(buffer);
    //Add the residual from last read
    if (ResidualString.size()!=0)BufString = ResidualString+BufString;
    ResidualString.clear();

    size_t foundnewline = BufString.find("\n");
    //  static  bool flag = false;

    GalilDataStruct GalilDataUnit;
    GalilDataStructD GalilDataUnitD;

    jj=0;
    while (foundnewline!=string::npos){
      stringstream iss (BufString.substr(0,foundnewline-1));
      // output returned by Galil is stored in the following variables

      iss >> Header;
      if(Header.compare("Galil")==0){
	//iss >> GalilDataUnit.TimeStamp;
	iss >> Time;
	if (i==0 && jj==0)Time0=Time;
	Time-=Time0;
	for (int j=0;j<4;j++){
	  iss >> GalilDataUnitD.PositionArray[j];
	}
	for (int j=0;j<4;j++){
	  iss >> GalilDataUnitD.VelocityArray[j];
	}
	for (int j=0;j<4;j++){
	  iss >> GalilDataUnitD.OutputVArray[j];
	}
	//Convert to INT
	GalilDataUnit.TimeStamp = INT(Time);
	for (int j=0;j<4;j++){
	  GalilDataUnit.PositionArray[j] = INT(GalilDataUnitD.PositionArray[j]);
	}
	for (int j=0;j<4;j++){
	  GalilDataUnit.VelocityArray[j] = INT(GalilDataUnitD.VelocityArray[j]);
	}
	for (int j=0;j<4;j++){
	  GalilDataUnit.OutputVArray[j] = INT(GalilDataUnitD.OutputVArray[j]);
	}
      }

      BufString = BufString.substr(foundnewline+1,string::npos);
      foundnewline = BufString.find("\n");
      jj++;
    }
    if (BufString.size()!=0){
      ResidualString = BufString;
    }

    //Update odb for monitoring
    db_set_value(hDB,0,"/Equipment/Galil/Monitors/Positions",&GalilDataUnit.PositionArray,sizeof(GalilDataUnit.PositionArray), 4 ,TID_INT); 
    db_set_value(hDB,0,"/Equipment/Galil/Monitors/Velocities",&GalilDataUnit.VelocityArray,sizeof(GalilDataUnit.VelocityArray), 4 ,TID_INT); 
    db_set_value(hDB,0,"/Equipment/Galil/Monitors/ControlVoltages",&GalilDataUnit.OutputVArray,sizeof(GalilDataUnit.OutputVArray), 4 ,TID_INT); 

    //Check commandand execute;
    INT command_size = sizeof(command);
    db_get_value(hDB,0,"/Equipment/Galil/ManualControl/Cmd",&command,&command_size,TID_INT,0);
    char CmdBuffer[500];
    if (command == 1){
      INT AbsPos[4];
      INT P_size = sizeof(AbsPos);
      db_get_value(hDB,0,"/Equipment/Galil/ManualControl/AbsPos",AbsPos,&P_size,TID_INT,0);
      sprintf(CmdBuffer,"PA %d,%d,%d,%d",AbsPos[0],AbsPos[1],AbsPos[2],AbsPos[3]);
      GCmd(g,CmdBuffer);
      GCmd(g,"BG");
      command=0;
      db_set_value(hDB,0,"/Equipment/Galil/ManualControl/Cmd",&command,sizeof(command), 1 ,TID_INT); 
    }else if (command == 2){
      INT RelPos[4];
      INT P_size = sizeof(RelPos);
      db_get_value(hDB,0,"/Equipment/Galil/ManualControl/RelPos",RelPos,&P_size,TID_INT,0);
      sprintf(CmdBuffer,"PR %d,%d,%d,%d",RelPos[0],RelPos[1],RelPos[2],RelPos[3]);
      GCmd(g,CmdBuffer);
      GCmd(g,"BG");
      command=0;
      db_set_value(hDB,0,"/Equipment/Galil/ManualControl/Cmd",&command,sizeof(command), 1 ,TID_INT); 
    }else if (command == 3){
      GCmd(g,"DP");
      command=0;
      db_set_value(hDB,0,"/Equipment/Galil/ManualControl/Cmd",&command,sizeof(command), 1 ,TID_INT); 
    }
    i++;
  }
  MonitorThreadActive = 0;
  db_set_value(hDB,0,"/Equipment/Galil/Monitors/Active",&MonitorThreadActive,sizeof(MonitorThreadActive), 1 ,TID_BOOL); 
}
