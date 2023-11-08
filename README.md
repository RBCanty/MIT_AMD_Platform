# MIT Accelerated Molecular Discovery Platform
This repo contains a copy of the python files for platform control as well as a snapshot of the property prediction model code

The master control network (MCN) comprises four computers, each running a System.py (MC, AH, LC, and SP)
A System can communicate over the Network with a Server.py
Within a System, there are threads for processing user input, communication over the network, and interpreting messages.
Each non-atomic message spawns a thread which performs the actions specified by the message.  
This may include spawning an instance of an API/Wrapper/Shell to control a subsystem (Ra, Lh, Pr, Sw, Fs, Ss, and Th).

* MC - Master Controller
    * Tasked with running scripts and thus coordinating all other devices
    * Ra - Robotic Arm
* AH - Automation Handler
    * Lh - Liquid Handler
* LC - Liquid Chromatograph
    * Sw - Shimadzu wrapper for HPLC-MS
        * as and fc for autosampler and fraction collector
* SP - Special Processes
    * Fs - FTIR shell
    * Ss - LPX storage unit
    * Th - Thermo Reactor
    * Pr - Plate reader
