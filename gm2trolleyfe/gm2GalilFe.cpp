/********************************************************************\

Name:         gm2GalilFe.cpp
Created by:   Matteo Bartolini
Modified by:  Joe Grange, Ran Hong

Contents:     readout code to talk to Galil motion control

$Id$

\********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "midas.h"
#include "mcstd.h"
#include "experim.h"
#include "gclib.h"
#include "gclibo.h"
#include <iostream>
#include <string>
#include <iomanip>
//#include "/home/galil/DAQ/midas/drivers/device/nulldev.h"
//#include "/home/galil/DAQ/midas/drivers/bus/null.h"
//#include "/home/galil/DAQ/midas/drivers/class/hv.h"
//#include "/home/galil/DAQ/midas/drivers/bus/rs232.h"
#include <termios.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/ioctl.h>
#include <sys/timeb.h>
#include <fstream>
#include <sstream>


#define GALIL_EXAMPLE_OK G_NO_ERROR //return code for correct code execution
#define GALIL_EXAMPLE_ERROR -100
using namespace std;

/* make frontend functions callable from the C framework */
#ifdef __cplusplus
extern "C" {
#endif

  ofstream myfile; 

  // i am defining some Galil libraries variables
  INT level1=2;
  float axes[3];
  INT setaxes[3];
  float speed[3];
  float acceleration[3];
  float torque[3];
  INT getaxes[3];
  HNDLE hDB, hkeyclient;
  char  name[32];
  int   size; //size of axes[3]
  INT size1; // size of setaxes[3]
  float TensionArray[2];
  float VelocityArray[2];
  float PositionArray[2];
  float OutputVArray[2];
  float BarcodeArray[3];

  int i;
  string s;
  int s1;
  GReturn b = G_NO_ERROR;
  int rc = GALIL_EXAMPLE_OK; //return code
  char buf[1023]; //traffic buffer
  char buf1[1024];
  //char buf2[1024];
  GCon g = 0; //var used to refer to a unique connection. A valid connection is nonzero.

  //----------------------------------------------------------

  ofstream GalilOutFile;

  /*-- Globals -------------------------------------------------------*/

  /* The frontend name (client name) as seen by other MIDAS clients   */
  char *frontend_name = "gm2GalilFe";
  /* The frontend file name, don't change it */
  char *frontend_file_name = __FILE__;

  /* frontend_loop is called periodically if this variable is TRUE    */
  BOOL frontend_call_loop = FALSE;

  /* a frontend status page is displayed with this frequency in ms */
  INT display_period = 3000;

  /* maximum event size produced by this frontend */
  INT max_event_size = 10000;

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

  INT read_galil_event(char *pevent, INT off);
  INT read_trigger_event(char *pevent, INT off);

  INT poll_event(INT source, INT count, BOOL test);
  INT interrupt_configure(INT cmd, INT source, POINTER_T adr);
  INT db_set_value(HNDLE hDB, HNDLE hKeyRoot, const char *key_name, const void *data, INT data_size, INT num_values, DWORD type); 
  INT db_find_key(HNDLE hdB, HNDLE hKey, const char *key_name, HNDLE *subhkey);
  INT cm_get_experiment_database(HNDLE *hDB, HNDLE *hKeyClient);


  /* device driver list */
  /*DEVICE_DRIVER hv_driver[] = {
    {"Dummy Device", nulldev, 16, null},
    {""}
    };*/


  /*-- Equipment list ------------------------------------------------*/


  EQUIPMENT equipment[] = {


    {"Galil",                /* equipment name */
      {3, 0,                   /* event ID, trigger mask */
	"SYSTEM",               /* event buffer */
	EQ_PERIODIC,            /* equipment type */
	0,                      /* event source */
	"MIDAS",                /* format */
	TRUE,                   /* enabled */
	RO_RUNNING | RO_TRANSITIONS |   /* read when running and on transitions */
	  RO_ODB,                 /* and update ODB */
	10,                  /* read every 0.01 sec */
	0,                      /* stop run after this event limit */
	0,                      /* number of sub events */
	60,                      /* log history, logged once per minute */
	"", "", "",},
      read_galil_event,       /* readout routine */
    },



    {""}
  };


#ifdef __cplusplus
}
#endif

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
  //b=GOpen("/dev/ttyUSB0 -t 1000 -s MG -d", &g);
  b=GOpen("192.168.1.13 -s ALL -t 1000 -d",&g);
  //GOpen("00:50:4c:38:19:AA -s ALL -t 1000 -d", &g);
  GInfo(g, buf, sizeof(buf)); //grab connection string
  cout << "buf is" << " "<<  buf << "\n";
  if (b==G_NO_ERROR){
    cout << "connection successful\n";
  }
  else {cout << "connection failed \n";}

  GTimeout(g,2000);//adjust timeout
  //int i = 0;
  //int s;

  //-------------end code to communicate with Galil------------------



  return SUCCESS;
}

/*-- Frontend Exit -------------------------------------------------*/

INT frontend_exit()
{
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
  char filename[1000];
  sprintf(filename,"/Users/rhong/gm2software/experiments/gm2Trolley/GalilOutput%04d.txt",RunNumber);
  GalilOutFile.open(filename,ios::out);
  //Load script
  GProgramDownload(g,"",0); //to erase prevoius programs
  b=GProgramDownloadFile(g,"/Users/rhong/gm2software/experiments/gm2Trolley/feedbacktestCurrent.dmc",0);
  GCmd(g, "XQ");
  return SUCCESS;
}

/*-- End of Run ----------------------------------------------------*/

INT end_of_run(INT run_number, char *error)

{
  GCmd(g,"AB");
  GalilOutFile.close();
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


  return 0;
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






INT read_galil_event(char *pevent, INT off){
  float *pdata, a;
  float *pspid;
  float *pacc;

  char buffer[500];
  char buffer1[500];
  char buffer2[500];
  hkeyclient=0;

  double Tension1 = 0.0;
  double Tension2 = 0.0;
  double Velocity1 = 0.0;
  double Velocity2 = 0.0;
  double Position1 = 0.0;
  double Position2 = 0.0;
  double OutputV1 = 0.0;
  double OutputV2 = 0.0;

  double Barcode1 = 0.0;
  double Barcode2 = 0.0;
  double Barcode3 = 0.0;

  rc = GCmdD(g,"MG @AN[1]",&Tension1);
  rc = GCmdD(g,"MG @AN[2]",&Tension2);
  rc = GCmdD(g,"TPA",&Position1);
  rc = GCmdD(g,"TPC",&Position2);
  rc = GCmdD(g,"TVA",&Velocity1);
  rc = GCmdD(g,"TVC",&Velocity2);
  rc = GCmdD(g,"TTA",&OutputV1);
  rc = GCmdD(g,"TTC",&OutputV2);
  rc = GCmdD(g,"MG @AN[3]",&Barcode1);
  rc = GCmdD(g,"MG @AN[4]",&Barcode2);
  rc = GCmdD(g,"MG @AN[5]",&Barcode3);
  TensionArray[0] = float(Tension1);
  TensionArray[1] = float(Tension2);
  VelocityArray[0] = float(Velocity1);
  VelocityArray[1] = float(Velocity2);
  PositionArray[0] = float(Position1);
  PositionArray[1] = float(Position2);
  OutputVArray[0] = float(OutputV1);
  OutputVArray[1] = float(OutputV2);
  BarcodeArray[0] = float(Barcode1);
  BarcodeArray[1] = float(Barcode2);
  BarcodeArray[2] = float(Barcode3);
/*
  rc = GMessage(g, buf1, sizeof(buf1));
//  cout << buf1 << endl;

  stringstream iss (buf1);
  // output returned by Galil is stored in the following variables
  iss >> axes[0];
  iss >> speed[0];
  iss >> torque[0];
  iss >> Analog[0];
  iss >> Analog[1];
  iss >> Analog[2];

  cm_get_experiment_database(&hDB, NULL);
  db_set_value(hDB, 0, "/Equipment/Galil/Variables/Position",&axes, sizeof(axes), 3, TID_FLOAT);
  db_set_value(hDB,0,"/Equipment/Galil/Variables/Speed",&speed,sizeof(speed),3,TID_FLOAT);
  db_set_value(hDB,0,"/Equipment/Galil/Variables/Acceleration",&acceleration,sizeof(acceleration),3,TID_FLOAT);
  db_set_value(hDB,0,"/Equipment/Galil/Variables/Torque",&torque,sizeof(torque),3,TID_FLOAT);
*/

  bk_init32(pevent);

  GalilOutFile<<TensionArray[0]<<" "<<TensionArray[1]<<" "<<VelocityArray[0]<<" "<<VelocityArray[1]<<" "<<PositionArray[0]<<" "<<PositionArray[1]<<" "<<OutputVArray[0]<<" "<<OutputVArray[1]<<" "<<BarcodeArray[0]<<" "<<BarcodeArray[1]<<" "<<BarcodeArray[2]<<endl;
  cout<<TensionArray[0]<<" "<<TensionArray[1]<<" "<<VelocityArray[0]<<" "<<VelocityArray[1]<<" "<<PositionArray[0]<<" "<<PositionArray[1]<<" "<<OutputVArray[0]<<" "<<OutputVArray[1]<<" "<<BarcodeArray[0]<<" "<<BarcodeArray[1]<<" "<<BarcodeArray[2]<<endl;

  /* create banks */
  bk_create(pevent, "GALI", TID_FLOAT, (void **)&pdata);
  for (int j=0;j<2;j++){
    *pdata++ = TensionArray[j];
  }
  for (int j=0;j<2;j++){
    *pdata++ = VelocityArray[j];
  }
  for (int j=0;j<2;j++){
    *pdata++ = PositionArray[j];
  }
  for (int j=0;j<2;j++){
    *pdata++ = OutputVArray[j];
  }
  for (int j=0;j<2;j++){
    *pdata++ = BarcodeArray[j];
  }

  bk_close(pevent,pdata);
  /* 
     bk_create(pevent, "SPID", TID_FLOAT, (void **)&pspid);
     for (int j=0;j<3;j++){
   *pspid++ = speed[j];
   }
   bk_close(pevent,pspid);

   bk_create(pevent,"ACCL", TID_FLOAT, (void **)&pacc);
   for(int j=0;j<3;j++){
   *pacc++ = acceleration[j];
   }
   bk_close(pevent,pacc);
   */
  return bk_size(pevent);
}


