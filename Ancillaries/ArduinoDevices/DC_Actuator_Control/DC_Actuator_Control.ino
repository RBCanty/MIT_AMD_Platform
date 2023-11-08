const byte numChars = 8;
char receivedChars[numChars];
char tempChars[numChars]; // temporary array for use when parsing

char DeviceID[1] = {0};
int onSwitchState = 0;
int offSwitchState = 0;
boolean newData = false;

#define onSwitchDirPin 5 // digital
#define onSwitchPwmPin 6 // analog
#define offSwitchDirPin 8
#define offSwitchPwmPin 9

void setup() {
  pinMode(onSwitchDirPin,OUTPUT);
  pinMode(onSwitchPwmPin,OUTPUT);
  pinMode(offSwitchDirPin,OUTPUT);
  pinMode(offSwitchPwmPin,OUTPUT);

  Serial.begin(9600);
}

void loop() {
  recvWithStartEndMarkers();
  if (newData == true) {
      strcpy(tempChars, receivedChars);
      parseData();
      showParsedData();
      newData = false;
  }
}
