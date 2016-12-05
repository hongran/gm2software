/********************************************************************\

  Name:         experim.h
  Created by:   ODBedit program

  Contents:     This file contains C structures for the "Experiment"
                tree in the ODB and the "/Analyzer/Parameters" tree.

                Additionally, it contains the "Settings" subtree for
                all items listed under "/Equipment" as well as their
                event definition.

                It can be used by the frontend and analyzer to work
                with these information.

                All C structures are accompanied with a string represen-
                tation which can be used in the db_create_record function
                to setup an ODB structure which matches the C structure.

  Created on:   Tue Apr  5 18:16:11 2016

\********************************************************************/

#define EXP_PARAM_DEFINED

typedef struct {
  char      comment[80];
} EXP_PARAM;

#define EXP_PARAM_STR(_name) const char *_name[] = {\
"[.]",\
"Comment = STRING : [80] Test",\
"",\
NULL }

#ifndef EXCL_GLOBAL

#define GLOBAL_PARAM_DEFINED

typedef struct {
} GLOBAL_PARAM;

#define GLOBAL_PARAM_STR(_name) const char *_name[] = {\
"",\
NULL }

#endif

#ifndef EXCL_GALIL

#define GALIL_COMMON_DEFINED

typedef struct {
  WORD      event_id;
  WORD      trigger_mask;
  char      buffer[32];
  INT       type;
  INT       source;
  char      format[8];
  BOOL      enabled;
  INT       read_on;
  INT       period;
  double    event_limit;
  DWORD     num_subevents;
  INT       log_history;
  char      frontend_host[32];
  char      frontend_name[32];
  char      frontend_file_name[256];
  char      status[256];
  char      status_color[32];
  BOOL      hidden;
} GALIL_COMMON;

#define GALIL_COMMON_STR(_name) const char *_name[] = {\
"[.]",\
"Event ID = WORD : 3",\
"Trigger mask = WORD : 0",\
"Buffer = STRING : [32] SYSTEM",\
"Type = INT : 1",\
"Source = INT : 0",\
"Format = STRING : [8] MIDAS",\
"Enabled = BOOL : y",\
"Read on = INT : 377",\
"Period = INT : 1000",\
"Event limit = DOUBLE : 0",\
"Num subevents = DWORD : 0",\
"Log history = INT : 1",\
"Frontend host = STRING : [32] localhost",\
"Frontend name = STRING : [32] Sample Frontend",\
"Frontend file name = STRING : [256] frontend.c",\
"Status = STRING : [256] Sample Frontend@localhost",\
"Status color = STRING : [32] greenLight",\
"Hidden = BOOL : n",\
"",\
NULL }

#endif

