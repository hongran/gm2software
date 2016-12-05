mkdir "Equipment/NMRProbes/Settings"
cd "Equipment/NMRProbes/Settings"
create BOOL "Global On"
set "Global On" false
mkdir "FPGA Timing"
mkdir "LO"
mkdir "Digitizer"
mkdir "Utilities"

cd "FPGA Timing"
mkdir "S1"
cd "S1"

create STRING "Configuration"
create DOUBLE "Mechanical Offset"
create DOUBLE "Mech Switch Duration"
create DOUBLE "RF Transmit Offset"
create DOUBLE "RF Trans Switch Duration"
create DOUBLE "Tomco Offset"
create DOUBLE "Amplifier Duration"
create BOOL "Tomco Enable"
create DOUBLE "RF Receive Offset"
create DOUBLE "RF Rec Switch Duration"
set "Configuration" "OFF"
set "Mechanical Offset" 0.0
set "Mech Switch Duration" 0.0
set "RF Transmit Offset" 0.0
set "RF Trans Switch Duration" 0.0
set "Tomco Offset" 0.0
set "Amplifier Duration" 0.0
set "Tomco Enable" false
set "RF Receive Offset" 0.0
set "RF Rec Switch Duration" 0.0
cd ..
copy "S1" "S2"
copy "S1" "S3"
copy "S1" "S4"

cd ..
cd  "LO"
create DOUBLE "Frequency"
create DOUBLE "BNC Voltage"
create STRING "BNC V-choice"
create BOOL "BNC Switch"
create DOUBLE "N-Type Voltage"
create STRING "N-Type V-choice"
create BOOL "N-Type Switch"
set "Frequency" 0
set "BNC Voltage" 0
set "BNC V-choice" "Vpp"
set "BNC Switch" false
set "N-Type Voltage" 0
set "N-Type V-choice" "Vpp"
set "N-Type Switch" false

cd ..
cd "Digitizer"
create INT "Struck ID"
create INT "Channel Number"
create INT "Number of Pulses"
create DOUBLE "Sampling Frequency"
create BOOL "External"
set "Struck ID" 0
set "Channel Number" 1
set "Number of Pulses" 1
set "Sampling Frequency" 1
set "External" false

cd ..
cd "Utilities"
create STRING "Debug Mode"
create INT "Verbosity"
create INT "Test Mode"
create DOUBLE "Delay Time"
create DOUBLE "PTS160 Frequency"
create STRING "Config Label"
create STRING "Comments"
set "Debug Mode" "OFF"
set "Verbosity" 0
set "Test Mode" 0
set "Delay Time" 0
set "PTS160 Frequency" 0
set "Config Label" "Label"
set "Comments" "Comments"
