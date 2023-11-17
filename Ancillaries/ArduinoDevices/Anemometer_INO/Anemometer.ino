const byte numChars = 8;
char receivedChars[numChars];
char tempChars[numChars]; // temporary array for use when parsing

char DeviceID[1] = {0};
float anemoVal = 0;
int readState = 0;
boolean newData = false;

#define anemoPin A0

void setup() {
  pinMode(anemoPin,INPUT);
  Serial.begin(9600);
}

void loop() {
  recvWithStartEndMarkers();
  if (newData == true) {
      strcpy(tempChars, receivedChars);
      parseData();
      //showParsedData();
      newData = false;
  }
}
