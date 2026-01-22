# MonaCalibration
Calibrating MoNA from e23033 and e23068.

ALL CODES written by Kevin Eisenberg.

This is a collection of some of my calibration codes. Others I have are not included here. 

Directory: 

kevinMonaQCal and kevinMonaQCalLOOP both are CHARGE CALIBRATION (assigning muon peaks to 20.5 MeVee). The parameters (linear fit) are stored in chargeCalibParams.

pulserCalib is PULSER CALIBRATION, using pulser data to extract ns/ch for each PMT, converting TDC into ns. The parameters (linear fit) are stored in pulserParams.

graphWallTDiffs finds the RELATIVE OFFSETS for a single wall (I ran this 9 times to fill the offsets csv file), ns relative to the top bar. The parameters (absolute offsets) are stored in wallOffsets.

posCalibrateWall is not a calibration, it is a visualization of a particular mona wall in a particular run, you should see what the wall looks like, illuminated by data.
