void onSwitch(){
  // Turn on the actuator:
  analogWrite(onSwitchPwmPin,127);
  
  // Set the direction pin low to push the switch:
  digitalWrite(onSwitchDirPin,LOW);
  delay(1200);

  // Set the direction pin high to return to extended position:
  digitalWrite(onSwitchDirPin,HIGH);
  delay(1200);

  //Turn off the actuator:
  analogWrite(onSwitchPwmPin,0);
}

void offSwitch(){
  // Turn on the actuator:
  analogWrite(offSwitchPwmPin,127);
  
  // Set the direction pin low to push the switch:
  digitalWrite(offSwitchDirPin,LOW);
  delay(1200);

  // Set the direction pin high to return to extended position:
  digitalWrite(offSwitchDirPin,HIGH);
  delay(1200);

  //Turn off the actuator:
  analogWrite(offSwitchPwmPin,0);
}
