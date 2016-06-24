typedef struct __attribute__((__packed__)) _Packet_Header {	// NOTE: Group fields are listed from MSB to LSB, need packing attribute to avoid padding and hence misalignment
  unsigned int	 fixed_value1;	// Should be 0x80007FFF, if this is properly aligned.
  unsigned int	 fixed_value2;	// Should be 0x80007FFF, if this is properly aligned.
  unsigned int	 fixed_value3;	// Should be 0x80007FFF, if this is properly aligned.
  unsigned short fixed_value4; // Should be 0xAAAA
  unsigned int	 packet_length;	// Packet Length in bytes (including header)
  unsigned int	 packet_number;	// Packet Number
  char		 status;			// Status bits. LSB is clock source
  char		 probe_id;		// Probe ID (only lower 5 bits are used)
  unsigned short nmr_samples;
  unsigned short rest[50];
} Packet_Header;

typedef struct _Trolley_Packet {
  Packet_Header	header;		// Packet header
  short*        waveform;	// Array of waveform data
} Trolley_Packet;

