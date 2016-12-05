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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/filesystem.hpp>

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
  char *frontend_name = "TestMagnetProbeFe";
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
  inline int get_odb_key_json(char *keyname, boost::property_tree::ptree& conf);


  /*-- Equipment list ------------------------------------------------*/


  EQUIPMENT equipment[] = {

    {"NMRProbes",                /* equipment name */
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

//Flags
bool ReadyToMove = false;
bool ReadyToRead = false;

//Globals
boost::property_tree::ptree tSettings;
string EOFstr = "end_of_file";
bool IsDebug = false;
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

  char key[256];
  snprintf(key,256,"/Equipment/NMRProbes/Settings");
  // Grab the odb data as json.
  int rc = get_odb_key_json(key, tSettings);

  //Print to setup files
  string global_fn     = "global_on_off";
  string delay_time_fn = "delay-time";
  string fpga_fn       = "pulse-data";
  string fg_fn         = "sg382";
  string adc_fn        = "struck_adc";
  string util_fn       = "utilities";
  string com_fn        = "comments";
  cout <<"XXXXXXXXXXXXXX"<<endl;
  string config_tag = tSettings.get<std::string>("Utilities.Config Label");
  string conf_path  = string("/Users/rhong/gm2software/TestMagnetFe/NMRProbes/") + string("input/configs/") + config_tag + string(".cnf");
  string prefix     = string("/Users/rhong/gm2software/TestMagnetFe/NMRProbes/input/configs/files/");
  string global_path= prefix + global_fn + "_" + config_tag + ".dat";
  string fpga_path  = prefix + fpga_fn   + "_" + config_tag + ".dat";
  string fg_path    = prefix + fg_fn     + "_" + config_tag + ".dat";
  string adc_path   = prefix + adc_fn    + "_" + config_tag + ".dat";
  string util_path  = prefix + util_fn   + "_" + config_tag + ".dat";
  string com_path   = prefix + com_fn    + "_" + config_tag + ".txt";

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

  if (ReadyToRead)return 1;
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
  ReadyToMove = true;
  return bk_size(pevent);
}

inline int get_odb_key_json(char *keyname, boost::property_tree::ptree& conf)
{
  // Data part
  HNDLE hDB, hkey;
  int size = 0;
  int bytes_written = 0;
  int rc = 0;

  char *json_buf = new char[0x8000]; // This should be large enough
  boost::property_tree::ptree pt;
  std::stringstream ss;

  // Get the experiment database handle.
  cm_get_experiment_database(&hDB, NULL);

  db_find_key(hDB, 0, keyname, &hkey);

  if (hkey) {
    db_copy_json_values(hDB, hkey, &json_buf, &size, &bytes_written, 1, 1, 1);

  } else {

    cm_msg(MERROR, "failed to load \"%s\" from ODB", keyname);
    return FE_ERR_ODB;
  }

  // Store it in a property tree.
  ss << json_buf;
  boost::property_tree::read_json(ss, conf);

  return SUCCESS;
}

int PrintToFileGlobal(string fn)
{
  string global_tag   = "global_on_off";
  string header       = "#    ID    value";
  bool global_state = tSettings.get("Global On",false);
  string global_str = global_tag + "\t";
  if (global_state){
    global_str = global_tag + "1";
  }else{
    global_str = global_tag + "0";
  }
  string eof_str = EOFstr+"\t99";
  if (IsDebug==0){
    ofstream globalFile;
    globalFile.open(fn.c_str(),ios::out);
    globalFile<<header<<endl;
    globalFile<<global_str<<endl;
    globalFile<<eof_str<<endl;
    globalFile.close();
  }else if (IsDebug==1){
    cout << fn<<endl;
    cout << header<<endl;
    cout << global_str<<endl;
    cout << eof_str<<endl;
  }
  return 0;
}
