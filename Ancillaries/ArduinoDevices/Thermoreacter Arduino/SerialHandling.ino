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
                receivedChars[ndx] = '\0'; // terminate the string
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


void parseData() {      // split the data into its parts

    char * strtokIndx; // this is used by strtok() as an index

    strtokIndx = strtok(tempChars,",");      // get the first part - the string
    strcpy(DeviceID, strtokIndx); // copy it to messageFromPC
    
    if (strcmp(DeviceID, "F")==0){
      strtokIndx = strtok(NULL, ","); // this continues where the previous call left off
      FanNum = atoi(strtokIndx);     // convert this part to an integer
      
      strtokIndx = strtok(NULL, ","); // this continues where the previous call left off
      FanSpeed = atoi(strtokIndx);     // convert this part to an integer

      RunFan(FanNum, FanSpeed);
    }

    if (strcmp(DeviceID, "L")==0){
      strtokIndx = strtok(NULL, ","); // this continues where the previous call left off
      LAState = atoi(strtokIndx);     // convert this part to an integer
      
      if (LAState == 2){//Release
        actuatorState = 1;
      }
      if (LAState == 1){//Press
        actuatorState = 2;
      }
    }

    if (strcmp(DeviceID, "D")==0){
      strtokIndx = strtok(NULL, ","); // this continues where the previous call left off
      DoorState = atoi(strtokIndx);     // convert this part to an integer

      if (DoorState == 2){
        CloseDoor();
      }
      if (DoorState == 1){
        OpenDoor();
      }
    }

    if (strcmp(DeviceID, "H")==0){
      strtokIndx = strtok(NULL, ","); // this continues where the previous call left off
      HeatVal = atoi(strtokIndx);     // convert this part to an integer

      RunHeat(HeatVal);
    }
}


void showParsedData() {
    Serial.print("Device ");
    Serial.println(DeviceID);
    Serial.print("Int 1 ");
    Serial.println(FanNum);
    Serial.print("Int 2 ");
    Serial.println(FanSpeed);
}
