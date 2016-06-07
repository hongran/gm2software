/********************************************************************\

  Name:         gm2trolleyfe.cpp
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
  INT allow;

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


/*-- Globals -------------------------------------------------------*/

/* The frontend name (client name) as seen by other MIDAS clients   */
char *frontend_name = "gm2trolleyfe";
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
     1000,                  /* read every 1 sec */
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
 myfile.open("/Users/rhong/gm2software/experiments/gm2Trolley/galilmove.dmc");
 myfile <<"#MOVE\nKIA=0.1\nKPA=103\nKDA=2268\nSPA=700\nACA=20000\nTLA=3\nKIB=1\nKPB=103\nKDB=2268\nSPB=8000\nACB=20000\nTLB=8\nKIC=1\nKPC=103\nKDC=2268\nSPC=8000\nACC=20000\nTLC=8\n#A\npos=_TPA\npos1=_TPB\npos2=_TPC\nsp=_TVA\nsp1=_TVB\nsp2=_TVC\nacc=_ACA\nacc1=_ACB\nacc2=_ACC\ntor=_TTA\ntor1=_TTB\ntor2=_TTC\nMG pos, pos1, pos2, sp, sp1, sp2, acc, acc1, acc2, tor, tor1, tor2\nWT500\nJP#A\nEN\n";
 myfile.close();
 


	//b=GOpen("/dev/ttyUSB0 -t 1000 -s MG -d", &g);
	b=GOpen("192.168.1.13 -s ALL -t 1000 -d",&g);
	//GOpen("00:50:4c:38:19:AA -s ALL -t 1000 -d", &g);
	GInfo(g, buf, sizeof(buf)); //grab connection string
	cout << "buf is" << " "<<  buf << "\n";
        if (b==G_NO_ERROR){
        cout << "connection successful\n";
        }
	else {cout << "connection failed \n";}

	GProgramDownload(g,"",0); //to erase prevoius programs
         b=GProgramDownloadFile(g,"/Users/rhong/gm2software/experiments/gm2Trolley/galilmove.dmc",0);
	GCmd(g, "XQ");
       	
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



   return SUCCESS;
}

/*-- Resuem Run ----------------------------------------------------*/

INT resume_run(INT run_number, char *error)
{ 



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
  cm_get_experiment_database(&hDB, NULL);
  int size2 = sizeof(getaxes);
  INT size3 = sizeof(allow);

  //db_find_key(hDB,0,"/Equipment/Galil/Variables",&hkeyclient);
  db_get_value(hDB,hkeyclient,"/Equipment/Galil/Variables/Setting",&getaxes,&size2,TID_INT, TRUE); // read values from Setting in the ODB and store it in vector getaxes
  db_get_value(hDB,hkeyclient,"/Equipment/Galil/Variables/Condition",&allow,&size3,TID_INT, 0);   //read values from Condition and store it in variable allow   

  //the variable allow is controlled by the user. Movement only starts if this variable is set to 1. As soon as commands are sent to Galil and the motion has begun this variable is set to 0 again
   sprintf(buffer,"PAA=%d",getaxes[0]);
   //sprintf(buffer1, "PAB=%d",getaxes[1]);
   //sprintf(buffer2, "PAC=%d",getaxes[2]);
   if(allow==1){
     //send command to Galil 
     GCmd(g,buffer); 
     GCmd(g,"BGA");
     /*
     GCmd(g,buffer1);  
     GCmd(g,"BGB");
     GCmd(g,buffer2);
     GCmd(g,"BGC");
     */

      //reset
      allow=0;
      db_set_value(hDB,0,"/Equipment/Galil/Variables/Condition",&allow,sizeof(allow),1,TID_INT);
    }

  rc = GMessage(g, buf1, sizeof(buf1));
  cout << buf1 << endl;
  
  stringstream iss (buf1);
  // output returned by Galil is stored in the following variables
  iss >> axes[0];
  iss >> axes[1];
  iss >> axes[2];
  iss >> speed[0];
  iss >> speed[1];
  iss >> speed[2];
  iss >> acceleration[0];
  iss >> acceleration[1];
  iss >> acceleration[2];
  iss >> torque[0];
  iss >> torque[1];
  iss >> torque[2];
 
  //update ODB
 cm_get_experiment_database(&hDB, NULL);

 db_set_value(hDB, 0, "/Equipment/Galil/Variables/Position",&axes, sizeof(axes), 3, TID_FLOAT);
 db_set_value(hDB,0,"/Equipment/Galil/Variables/Speed",&speed,sizeof(speed),3,TID_FLOAT);
 db_set_value(hDB,0,"/Equipment/Galil/Variables/Acceleration",&acceleration,sizeof(acceleration),3,TID_FLOAT);
 db_set_value(hDB,0,"/Equipment/Galil/Variables/Torque",&torque,sizeof(torque),3,TID_FLOAT);


   bk_init32(pevent);

   cout<<" positions 1,2,3: "<<axes[0]<<","<<axes[1]<<","<<axes[2]<<endl;

   /* create banks */
   bk_create(pevent, "GALI", TID_FLOAT, (void **)&pdata);
   for (int j=0;j<3;j++){
     *pdata++ = axes[j];
   }
   for (int j=0;j<3;j++){
     *pdata++ = speed[j];
   }
   for (int j=0;j<3;j++){
     *pdata++ = acceleration[j];
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
   

