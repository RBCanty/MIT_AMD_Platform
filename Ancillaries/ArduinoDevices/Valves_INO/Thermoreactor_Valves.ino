#include <AccelStepper.h>

const byte numChars = 8;
char receivedChars[numChars];
char tempChars[numChars]; // temporary array for use when parsing

char DeviceID[1] = {0};
int Valve1State = 0;
int Valve2State = 0;
int Valve3State = 0;
boolean newData = false;

//Stepper Driver Pins and Init
#define MOT1_STEP 7
#define MOT1_DIR 6
#define MOT2_STEP 5
#define MOT2_DIR 4
#define MOT3_STEP 3
#define MOT3_DIR 2

#define MOT1_SLEEP 8
#define MOT2_SLEEP 9
#define MOT3_SLEEP 10

int stepsPerRevolution = 3400;

void setup() {
  pinMode(MOT1_STEP,OUTPUT);
  pinMode(MOT1_DIR,OUTPUT);
  pinMode(MOT2_STEP,OUTPUT);
  pinMode(MOT2_DIR,OUTPUT);
  pinMode(MOT3_STEP,OUTPUT);
  pinMode(MOT3_DIR,OUTPUT);
  
  pinMode(MOT1_SLEEP,OUTPUT);
  digitalWrite(MOT1_SLEEP,LOW);
  pinMode(MOT2_SLEEP,OUTPUT);
  digitalWrite(MOT2_SLEEP,LOW);
  pinMode(MOT3_SLEEP,OUTPUT);
  digitalWrite(MOT3_SLEEP,LOW);

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
