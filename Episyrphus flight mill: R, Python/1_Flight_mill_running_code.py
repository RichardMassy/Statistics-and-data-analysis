#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 11:23:47 2019
Edit 05.10.20: filename number of recording
Edit 07.10.20: removed os chdir as unneeded
Edit 19.07.21 3.0 added second channel, reduced number of mills to 4,
    csv output file now writing in separate columns, log.flush replaced by
    line buffering option 1
Edit 27.08.21 removed switch time requirement, added forward and reverse
    directions, now triggers on 0 as well as 1,
    columns = Mill, time, direction, updates
Edit 10.09.21 changed timeit package to time, added time.sleep to loop
Edit 18.09.21 reduced sleep time to 0.0003, added distance counter to print,
    added b2w or w2b to file output

@author: Richard Massy
"""
import glob, time, RPi.GPIO as GPIO

mill = int(input("Enter mill number (1-4): "))

"""Setup GPIO"""
pins = {1:{"F":31,"R":33},2:{"F":35,"R":37},
        3:{"F":32,"R":36},4:{"F":38,"R":40}}[mill] # GPI0 pins - mills
GPIO.setmode(GPIO.BOARD)
GPIO.setwarnings(False)
# disable RuntimeWarning for pins with physical pull up resistors (3 and 5)
GPIO.setup(pins["F"], GPIO.IN, pull_up_down=GPIO.PUD_UP)
GPIO.setup(pins["R"], GPIO.IN, pull_up_down=GPIO.PUD_UP)
# Select pin and pull up to high level(3.3V)

"""Filename elements"""
fileID = input("Input file name, leave blank for timestamp: ")
rec_num = len(glob.glob("*mill{}*.csv".format(mill)))+1 # Recording number
starttime = time.perf_counter()
if len(fileID) == 0:
    fileID = int(starttime)
    
"""Logging"""
log = open('{0}-mill{1}({2}).csv'.format(fileID,mill,rec_num),'a',buffering=1)
log.write("Mill,Time,Direction,Updates\n")
toggleF = toggleR = triggeredF = triggeredR = dist = 1 # toggle 0 = open

print("Mill {} recording: Happy flying!".format(mill))

while time.perf_counter() - starttime < 14400: # 14400s = 4 hours
    if not GPIO.input(pins["F"]) == toggleF: # Forward
        if not GPIO.input(pins["R"]) == toggleF and triggeredF == 0:
            zeit = str(time.perf_counter() - starttime)
            dist += 1
            print ("mill{0}: {1} F {2}".format(mill,dist,zeit))
            log.write("{0},{1},F,{2}\n".format(mill,zeit,
                                               ["b2w","w2b"][toggleF == 1]))
            triggeredF = 1
        toggleF = GPIO.input(pins["F"])
        triggeredR = 0
    if not GPIO.input(pins["R"]) == toggleR: # Reverse
        if not GPIO.input(pins["F"]) == toggleR and triggeredR == 0:
            zeit = str(time.perf_counter() - starttime)
            dist += 1
            print ("mill{0}: {1} R {2}".format(mill,dist,zeit))
            log.write("{0},{1},R,{2}\n".format(mill,zeit,
                                               ["b2w","w2b"][toggleR == 1]))
            triggeredR = 1
        toggleR = GPIO.input(pins["R"])
        triggeredF = 0
    time.sleep(0.0002) # Pause every cycle to reduce CPU usage

"""Shutdown"""
log.write(",,,,,Time limit reached")
log.close()
GPIO.cleanup() # Release resource

print ('Programme finished')