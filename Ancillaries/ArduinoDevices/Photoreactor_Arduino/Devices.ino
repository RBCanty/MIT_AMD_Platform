void ledOn(int ledVal){
  analogWrite(ledPWM,255-ledVal);
  digitalWrite(ledFan,HIGH);
}

void ledOff(){
  analogWrite(ledPWM,255);
  digitalWrite(ledFan,LOW);
}

void openDoor(){
  // Set the spinning direction clockwise to open:
  digitalWrite(doorDir, HIGH);

  // Spin the stepper motor 1 revolution slowly:
  for (int i = 0; i < stepsPerRevolution; i++) {
    // These four lines result in 1 step:
    digitalWrite(doorStep, HIGH);
    delayMicroseconds(1500);
    digitalWrite(doorStep, LOW);
    delayMicroseconds(1500);
  }
}

void closeDoor(){
  //Set spinning direction counterclockwise to close:
  digitalWrite(doorDir, LOW);

  //Spin the stepper motor 1 revolution slowly:
  for (int i = 0; i < stepsPerRevolution; i++) {
    digitalWrite(doorStep, HIGH);
    delayMicroseconds(1500);
    digitalWrite(doorStep, LOW);
    delayMicroseconds(1500);
  }
}

void peltierOn(){
  digitalWrite(peltier,HIGH);
  digitalWrite(peltFan,HIGH);
}

void peltierOff(){
  digitalWrite(peltier,LOW);
  digitalWrite(peltFan,LOW);
}

void photoRead(){
  float photoVal = analogRead(photo);
  Serial.println(photoVal);
}
