void read(){
  anemoVal = analogRead(anemoPin);
  anemoVal = anemoVal * 0.0049;
  anemoVal = (anemoVal-0.4)/0.00494 * 0.1;
  Serial.println(anemoVal);
}
