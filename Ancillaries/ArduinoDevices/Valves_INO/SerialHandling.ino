void recvWithStartEndMarkers() {
  static boolean recvInProgress = false;
  static byte ndx = 0;
  char startMarker = '<';
  char endMarker = '>';
  char rc;

  while (Serial.available() > 0 && newData == false) {
    rc = Serial.read();

    if (recvInProgress == true) {
      if (rc != endMarker) {
        receivedChars[ndx] = rc;
        ndx++;
        if (ndx >= numChars) {
          ndx = numChars - 1;
        }
      }
      else {
        receivedChars[ndx] = '\0'; //terminate the string
        recvInProgress = false;
        ndx = 0;
        newData = true;
      }
    }
    else if (rc == startMarker) {
      recvInProgress = true;
    }
  }
}

void parseData() { //split the data into its parts
  char * strtokIndx; //this is used by strtok() as an index
  strtokIndx = strtok(tempChars,","); //get the first part - the string
  strcpy(DeviceID, strtokIndx); //copy it to messageFromPC
    
  if (strcmp(DeviceID, "L")==0){
    strtokIndx = strtok(NULL, ","); //this continues where the previous call left off
    Valve1State = atoi(strtokIndx); //convert this part to an integer
    if (Valve1State == 2){
      CloseValve1();
    }
    if (Valve1State == 1){
      OpenValve1();
    }
  }

  if (strcmp(DeviceID, "M")==0){
    strtokIndx = strtok(NULL, ","); //this continues where the previous call left off
    Valve2State = atoi(strtokIndx); //convert this part to an integer
    if (Valve2State == 2){
      CloseValve2();
    }
    if (Valve2State == 1){
      OpenValve2();
    }
  }

  if (strcmp(DeviceID, "N")==0){
    strtokIndx = strtok(NULL, ","); //this continues where the previous call left off
    Valve3State = atoi(strtokIndx); //convert this part to an integer
    if (Valve3State == 2){
      CloseValve3();
    }
    if (Valve3State == 1){
      OpenValve3();
    }
  }
}

void showParsedData() {
  Serial.print("Device ");
  Serial.println(DeviceID);
}
