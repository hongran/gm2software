#ifndef __TrolleyHeader_h
#define __TrolleyHeader_h

typedef struct __attribute__((__packed__)) _Trolley_Header {	
  // NOTE: Group fields are listed from MSB to LSB, need packing attribute to avoid padding and hence misalignment
  unsigned int	 fixed_value1;	// Should be 0x80007FFF, if this is properly aligned.
  unsigned int	 fixed_value2;	// Should be 0x80007FFF, if this is properly aligned.
  unsigned int	 fixed_value3;	// Should be 0x80007FFF, if this is properly aligned.
  unsigned short fixed_value4; // Should be 0xAAAA
  unsigned int	 packet_length;	// Packet Length in bytes (including header)
  unsigned int	 packet_number;	// Packet Number
  char		 status;			// Status bits. LSB is clock source
  char		 probe_id;		// Probe ID (only lower 5 bits are used)
  unsigned short nmr_samples;
  unsigned short barcode_samples;
  unsigned short barcode_offset;
  unsigned short rest[48];
} Trolley_Header;


typedef struct _Trolley_Packet {
  Trolley_Header	header;		// Packet header
  short*        waveform;	// Array of waveform data
} Trolley_Packet;

#endif
