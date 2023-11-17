#include <AccelStepper.h>

const byte numChars = 8;
char receivedChars[numChars];
char tempChars[numChars]; // temporary array for use when parsing

char DeviceID[1] = {0};
int ledVal = 0;
int doorMotState = 0;
int peltierState = 0;
int photoState = 0;
boolean newData = false;

//Define pins
#define ledPWM 3
#define doorStep 6
#define doorDir 7
#define peltier 8
#define photo A0
#define ledFan A1
#define peltFan A2

int stepsPerRevolution = 1000;

void setup() {
  pinMode(ledPWM,OUTPUT);
  pinMode(doorStep,OUTPUT);
  pinMode(doorDir,OUTPUT);
  pinMode(peltier,OUTPUT);
  pinMode(ledFan,OUTPUT);
  pinMode(peltFan,OUTPUT);
  pinMode(photo,INPUT);

  analogWrite(ledPWM,255);
  digitalWrite(doorStep,LOW);
  digitalWrite(doorDir,LOW);
  digitalWrite(peltier,LOW);
  digitalWrite(ledFan,LOW);
  digitalWrite(peltFan,LOW);

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
