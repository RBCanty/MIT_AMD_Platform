void OpenValve1(){
  // Turn on the motor:
  digitalWrite(MOT1_SLEEP, HIGH);
  
  // Set the spinning direction clockwise to open:
  digitalWrite(MOT1_DIR, LOW);

  // Spin the stepper motor 1 revolution slowly:
  for (int i = 0; i < stepsPerRevolution; i++) {
    // These four lines result in 1 step:
    digitalWrite(MOT1_STEP, HIGH);
    delayMicroseconds(1500);
    digitalWrite(MOT1_STEP, LOW);
    delayMicroseconds(1500);
  }

  //Turn off the motor:
  digitalWrite(MOT1_SLEEP, LOW);
}

void CloseValve1(){
  // Turn on the motor:
  digitalWrite(MOT1_SLEEP, HIGH);
  
  //Set spinning direction counterclockwise to close:
  digitalWrite(MOT1_DIR, HIGH);

  //Spin the stepper motor 1 revolution slowly:
  for (int i = 0; i < stepsPerRevolution; i++) {
    digitalWrite(MOT1_STEP, HIGH);
    delayMicroseconds(1500);
    digitalWrite(MOT1_STEP, LOW);
    delayMicroseconds(1500);
  }

  //Turn off the motor:
  digitalWrite(MOT1_SLEEP, LOW);
}

void OpenValve2(){
  // Turn on the motor:
  digitalWrite(MOT2_SLEEP, HIGH);
  
  // Set the spinning direction clockwise to open:
  digitalWrite(MOT2_DIR, LOW);

  // Spin the stepper motor 1 revolution slowly:
  for (int i = 0; i < stepsPerRevolution; i++) {
    // These four lines result in 1 step:
    digitalWrite(MOT2_STEP, HIGH);
    delayMicroseconds(1500);
    digitalWrite(MOT2_STEP, LOW);
    delayMicroseconds(1500);
  }

  //Turn off the motor:
  digitalWrite(MOT2_SLEEP, LOW);
}

void CloseValve2(){
  // Turn on the motor:
  digitalWrite(MOT2_SLEEP, HIGH);
  
  //Set spinning direction counterclockwise to close:
  digitalWrite(MOT2_DIR, HIGH);

  //Spin the stepper motor 1 revolution slowly:
  for (int i = 0; i < stepsPerRevolution; i++) {
    digitalWrite(MOT2_STEP, HIGH);
    delayMicroseconds(1500);
    digitalWrite(MOT2_STEP, LOW);
    delayMicroseconds(1500);
  }

  //Turn off the motor:
  digitalWrite(MOT2_SLEEP, LOW);
}

void OpenValve3(){
  // Turn on the motor:
  digitalWrite(MOT3_SLEEP, HIGH);
  
  // Set the spinning direction clockwise to open:
  digitalWrite(MOT3_DIR, LOW);

  // Spin the stepper motor 1 revolution slowly:
  for (int i = 0; i < stepsPerRevolution; i++) {
    // These four lines result in 1 step:
    digitalWrite(MOT3_STEP, HIGH);
    delayMicroseconds(1500);
    digitalWrite(MOT3_STEP, LOW);
    delayMicroseconds(1500);
  }

  //Turn off the motor:
  digitalWrite(MOT3_SLEEP, LOW);
}

void CloseValve3(){
  // Turn on the motor:
  digitalWrite(MOT3_SLEEP, HIGH);
  
  //Set spinning direction counterclockwise to close:
  digitalWrite(MOT3_DIR, HIGH);

  //Spin the stepper motor 1 revolution slowly:
  for (int i = 0; i < stepsPerRevolution; i++) {
    digitalWrite(MOT3_STEP, HIGH);
    delayMicroseconds(1500);
    digitalWrite(MOT3_STEP, LOW);
    delayMicroseconds(1500);
  }

  //Turn off the motor:
  digitalWrite(MOT3_SLEEP, LOW);
}
