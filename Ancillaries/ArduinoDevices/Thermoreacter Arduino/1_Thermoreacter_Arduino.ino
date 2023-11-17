#include <AccelStepper.h>

const byte numChars = 8;
char receivedChars[numChars];
char tempChars[numChars]; // temporary array for use when parsing

// variables to hold the parsed data
char DeviceID[1] = {0};
int DoorState = 0;
int FanNum = 0;
int FanSpeed = 0;
int HeatVal = 0;
int LAState = 0;

boolean newData = false;

//Stepper Driver Pins and Init
#define ST_STEP 2
#define ST_DIR 3
#define ST_EN 4
AccelStepper doorStepper = AccelStepper(1, ST_STEP, ST_DIR);
const int open_pos = 100; //Open position of door in steps
const int closed_pos = 0; //Closed position of door in steps

//Linear Actuator Pins and Settings
#define LA_PWM 5
#define LA_INA 7
#define LA_INB 8
#define LA_POT A7



#define actuation_speed 50 //Speed to move linear actuator (0-255)
#define pressed_pos 327  //Closed position. Higher number makes piston extend further. 350 for no springs. 1 step = 0.0042"
#define released_pos 40  //Open position. Should not need adjustment

// start with 290 for spings

byte actuatorState = 0;

//Mosfet Pins
#define Q1_PWM 6
#define Q2_PWM 9
#define Q3_PWM 10

void setup() {
    //Stepper Driver Setup
    pinMode(ST_EN,OUTPUT);
    digitalWrite(ST_EN,HIGH);
    doorStepper.setMaxSpeed(1000);
    doorStepper.setAcceleration(100);

    //Drive stepper to closed position to get reference position
    doorStepper.setSpeed(-200);
    doorStepper.runSpeed();
    delay(2000);
    doorStepper.setCurrentPosition(0);
    
    //Linear Actuator Pins
    pinMode(LA_PWM,OUTPUT);
    pinMode(LA_INA,OUTPUT);
    pinMode(LA_INB,OUTPUT);
    digitalWrite(LA_INA,LOW);
    digitalWrite(LA_INB,LOW);
    analogWrite(LA_PWM,0);
    
    pinMode(LA_POT,INPUT);
    //Mosfet Pins
    pinMode(Q1_PWM,OUTPUT);
    pinMode(Q2_PWM,OUTPUT);
    pinMode(Q3_PWM,OUTPUT);
    
    pinMode(LED_BUILTIN, OUTPUT);
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
    runActuator();
}
