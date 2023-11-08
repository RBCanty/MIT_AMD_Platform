void RunFan(int FanNum,int FanSpeed){
  if (FanNum == 1){
    analogWrite(Q1_PWM, FanSpeed);
  }
  if (FanNum == 2){
    analogWrite(Q2_PWM, FanSpeed);
  }
}

void runActuator(){
  switch(actuatorState){
    
    case 0://Do Nothing
      analogWrite(LA_PWM, 0);
    break;
    
    case 1://Release
      digitalWrite(LA_INA,LOW);
      digitalWrite(LA_INB,HIGH);
      if (analogRead(LA_POT) > released_pos){
        analogWrite(LA_PWM, actuation_speed);
      }
      else{
        analogWrite(LA_PWM, 0);
        actuatorState = 0;
        digitalWrite(LED_BUILTIN,LOW);
      }
    break;
    
    case 2://Press
      digitalWrite(LA_INA,HIGH);
      digitalWrite(LA_INB,LOW);
      if (analogRead(LA_POT) < pressed_pos){
        analogWrite(LA_PWM, actuation_speed);
      }
      else{
        analogWrite(LA_PWM, 0);
        actuatorState = 0;
        digitalWrite(LED_BUILTIN,LOW);
      }
    break;
  }
  
}

void OpenDoor(){
  digitalWrite(LED_BUILTIN,HIGH);
  digitalWrite(ST_EN,LOW);
  doorStepper.moveTo(open_pos);
  doorStepper.runToPosition();
  digitalWrite(ST_EN,HIGH);
  digitalWrite(LED_BUILTIN,LOW);
}

void CloseDoor(){
  digitalWrite(LED_BUILTIN,HIGH);
  digitalWrite(ST_EN,LOW);
  doorStepper.moveTo(closed_pos);
  doorStepper.runToPosition();
  digitalWrite(ST_EN,HIGH);
  digitalWrite(LED_BUILTIN,LOW);
}

void RunHeat(int HeatVal){
  analogWrite(Q3_PWM, HeatVal);
}
