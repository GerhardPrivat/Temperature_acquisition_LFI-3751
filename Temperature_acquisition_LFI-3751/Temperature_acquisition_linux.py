# Coded by AO

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# +-+-+- PROGRAM FOR LINUX +-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

import numpy as np
import matplotlib.pyplot as plt
import serial
import time
import datetime
import os
import threading
import multiprocessing
import ctypes
import sys
import operator 
import warnings

warnings.filterwarnings("ignore","Item size computed from the PEP 3118 buffer format string does not match the actual item size.")


# +-+-+ GENERAL REMARKS +-+-+
# ---------------------------
# The plotting has to be carried out by the main process. The rendering is not working otherwise. This
# is an intrinsic Python problem.
# On Windows machines, the main process has to be included in an 'if __name__ == '__main__'' statement.
# This is not necessary on other systems.
# If the process of waiting for user input from the command line is outsourced into the subprocess,
# a copy of stdin has to be created via 'newstdin = os.fdopen(os.dup(sys.stdin.fileno()))'. Under Windows,
# this has to be done within the subprocess. It is not possible to forward the file descriptor as a
# parameter to the second subprocess. This cost me some hours.


# +-+-+- CHILD PROCESS +-+-+-
# ---------------------------
# The child process takes care of communicating with the temperature controller and user input.
# Both actions are distributed on two different threads.
# The child process is defined as a single function called readingwriting
def readingwriting(timeint, countmax, newstdin, continuous, datalist, flag, settemp, Pvalue, Ivalue, Dvalue, settempflag, PIDflag, Pflag, Iflag, Dflag, OUTPUTflag, OUTPUTOFFflag, OUTPUTONflag):
    
    # -- SERIAL CONNECTION --
    ser = serial.Serial(
        port = "/dev/ttyUSB0",
        baudrate = 19200,
        parity = serial.PARITY_NONE,
        stopbits = serial.STOPBITS_ONE,
        bytesize = serial.EIGHTBITS,
        timeout=1
    )
    
    
    
    # -- FUNCTIONS FOR CHILD PROCESS --
    # Function for calculation the Frame Check Sum accordingly to the manual of the temperature controller.
    # The FCS has to be appended at every command.
    def returnfcs(string):
        FCS = 0
        for ite in range(len(string)):
            FCS = operator.xor(FCS,ord(string[ite]))
        return FCS
    
    
    # Function that is executed in the first thread. But only if the continuous variable is set to "n" at the beginning in the main process.
    # This function (thread) runs forever except the flag is set. Since the size of the storage array is fixed, the oldest data points
    # are going to be deleted.
    # UPDATE: After switching to the datamanager object, an array fixed in size is not needed anymore. But still, old data is delete in order
    # to reduce the need for resources.
    def recordingtemperature1(timeint, countmax, ser, datalist):
        
        count = 0
        starttime1 = time.time()
        
		timearray = np.zeros(countmax,dtype="float")
		temperaturearray = np.zeros(countmax,dtype="float")
		timestamparray = np.zeros(countmax,dtype="float")
	
        while True:
            if (flag.value == 0):
                starttime2 = time.time()
				# If a different set temperature is entered in the second thread, this if statement changes the actual set temperature of the controller.
				# But this is only done if the entered set temperature differs from the previous set temperature.
				if (settempflag.value == 1):
					ser.write("!101103+000.00026")
					acsettemp = float(ser.read(21)[10:17])
					if (acsettemp != settemp.value):
						command = "!101203+%07.3f" %settemp.value
						ser.write(command+hex(returnfcs(command))[2:])
						ser.read(21)
					settempflag.value = 0
				# This if statement switches the output on or off.
				if (OUTPUTflag.value == 1):
					if (OUTPUTOFFflag.value == 1):
						ser.write("!101251+000.00022")
						ser.read(21)
						OUTPUTOFFflag.value = 0
					elif (OUTPUTONflag.value == 1):
						ser.write("!101251+000.00123")
						ser.read(21)
						OUTPUTONflag.value = 0
					OUTPUTflag.value = 0
				# This if statement changes the PID parameters. Here, a nested statement is used to avoid checking three if statements 
				# if no parameter has been altered. 
				if (PIDflag.value == 1):
					if (Pflag.value == 1):
						Pcommand = "!101210+%07.3f" %Pvalue.value
						ser.write(Pcommand+hex(returnfcs(Pcommand))[2:])
						ser.read(21)
						Pflag.value = 0
					elif (Iflag.value == 1):
						Icommand = "!101211+%07.3f" %Ivalue.value
						ser.write(Icommand+hex(returnfcs(Icommand))[2:])
						ser.read(21)
						Iflag.value = 0
					elif (Dflag.value == 1):
						Dcommand = "!101212+%07.3f" %Dvalue.value
						ser.write(Dcommand+hex(returnfcs(Dcommand))[2:])
						ser.read(21)
						Dflag.value = 0
					PIDflag.value = 0
		
				ser.write("!101101+000.00024")
				temp = float(ser.read(21)[10:17])                
					
				# The if statement is executed if there are less data points available than supposed to be shown, i.e. when the time window 
				# that is supposed to be shown is larger than the time for which the program has been running.
				# The else statement is executed if the time for which the program has been running is greater than the plotted time frame. 
				# In this case, the oldest datapoints are delete via np.roll and entering a new datapoint for the oldest one. 
                if (count < countmax):
                    timearray[count] = time.time()-starttime1
                    temperaturearray[count] = temp
		    		timestamparray[count] = starttime2
						
                    count += 1
								
                else:
                    timearray = np.roll(timearray,len(timearray)-1)
                    temperaturearray = np.roll(temperaturearray,len(temperaturearray)-1)
					timestamparray = np.roll(timestamparray,len(timestamparray)-1)
                    timearray[countmax-1] = time.time()-starttime1
                    temperaturearray[countmax-1] = temp
					timestamparray[countmax-1] = starttime2
			
				datalist.append([timearray,temperaturearray,timestamparray])
				
				if (count >= 3):
					del datalist[0]
							
                try:
                    time.sleep(timeint-time.time()+starttime2)
                except:
                    time.sleep(2*timeint-time.time()+starttime2)
            else:
                break
    
    
    # Function that is executed in the first thread. But only if the continuous variable is set to "y" at the beginning in the main process.
    # This function (thread) stops when the storage array is completely filled. No data points are going to be deleted.
    # This function (thread) can also be stopped by entering "exit".
    # For further comments refer to recordingtemperature1.
    def recordingtemperature2(timeint, countmax, ser, datalist):
        
        count = 0
        starttime1 = time.time()
        
        timearray = np.zeros(countmax,dtype="float")
        temperaturearray = np.zeros(countmax,dtype="float")
        timestamparray = np.zeros(countmax,dtype="float")
	
        while True:
            if (flag.value == 0):
                starttime2 = time.time()
                if (settempflag.value == 1):
                	ser.write("!101103+000.00026")
                	acsettemp = float(ser.read(21)[10:17])
                	if (acsettemp != settemp.value):
                        command = "!101203+%07.3f" %settemp.value
                        ser.write(command+hex(returnfcs(command))[2:])
                        ser.read(21)
                    settempflag.value = 0
                if (OUTPUTflag.value == 1):
                    if (OUTPUTOFFflag.value == 1):
                        ser.write("!101251+000.00022")
                        ser.read(21)
                        OUTPUTOFFflag.value = 0
                    elif (OUTPUTONflag.value == 1):
                        ser.write("!101251+000.00123")
                        ser.read(21)
                        OUTPUTONflag.value = 0
                    OUTPUTflag.value = 0
                if (PIDflag.value == 1):
                    if (Pflag.value == 1):
                        Pcommand = "!101210+%07.3f" %Pvalue.value
                        ser.write(Pcommand+hex(returnfcs(Pcommand))[2:])
                        ser.read(21)
                        Pflag.value = 0
                    elif (Iflag.value == 1):
                        Icommand = "!101211+%07.3f" %Ivalue.value
                        ser.write(Icommand+hex(returnfcs(Icommand))[2:])
                        ser.read(21)
                        Iflag.value = 0
                    elif (Dflag.value == 1):
                        Dcommand = "!101212+%07.3f" %Dvalue.value
                        ser.write(Dcommand+hex(returnfcs(Dcommand))[2:])
                        ser.read(21)
                        Dflag.value = 0
                    PIDflag.value = 0

                ser.write("!101101+000.00024")
                temp = float(ser.read(21)[10:17])
            
                if (count < countmax):
                    timearray[count] = time.time()-starttime1
                    temperaturearray[count] = temp
                    timestamparray[count] = starttime2
                
                    count += 1
                        
                else:
                    flag.value = 1
                    break
                
                datalist.append([timearray,temperaturearray,timestamparray])
		
                if (count >= 3):
                    del datalist[0]
		    
                try:
                    time.sleep(timeint-time.time()+starttime2)
                except:
                    time.sleep(2*timeint-time.time()+starttime2)
            else:
                break
        
        
    # Function that is executed in the second thread. The function waits for user input.
    # Possible input:
    #  - new set temperature 
    #  - switch output on or off
    #  - change PID parameters
    #  - increase or decrease the temperature in 0.1 or 0.01 steps
    #  - exit for quitting
    # The concept with flags is used so that certain code parts are triggered in the other thread.
    def readtemperature(newstdin):
    
        while True:
            print "   - Current set temperature: "+str(settemp.value)+" C"
            print("   - New set temperature: "),
            settemp_new = newstdin.readline()
            settemp_new = settemp_new[:len(settemp_new)-1]
                        
            if (settemp_new == "exit"):
                flag.value = 1
                newstdin.close()
                break
            elif (settemp_new == "OUTPUTOFF"):
                OUTPUTflag.value = 1
                OUTPUTOFFflag.value = 1
            elif (settemp_new == "OUTPUTON"):
                OUTPUTflag.value = 1
                OUTPUTONflag.value = 1
            elif (settemp_new[:1] == "P"):
                Pvalue.value = int(settemp_new[1:])
                PIDflag.value = 1
                Pflag.value = 1
            elif (settemp_new[:1] == "I"):
                Ivalue.value = float(settemp_new[1:])
                PIDflag.value = 1
                Iflag.value = 1
            elif (settemp_new[:1] == "D"):
                Dvalue.value = int(settemp_new[1:])
                PIDflag.value = 1
                Dflag.value = 1
            elif (settemp_new == "+"):
                settemp.value += 0.01
                settempflag.value = 1
            elif (settemp_new == "-"):
                settemp.value -= 0.01
                settempflag.value = 1
            elif (settemp_new == "++"):
                settemp.value += 0.1
                settempflag.value = 1
            elif (settemp_new == "--"):
                settemp.value -= 0.1
                settempflag.value = 1
            else:
                try:
                    settemp.value = float(settemp_new)
                    settempflag.value = 1
                except:
                    pass



    # -- ACTUAL CHILD PROCESS --
    # Read out the initial set temperature:
    ser.write("!101103+000.00026")
    settemp.value = float(ser.read(21)[10:17])
                
    # Start both actions in tow different threads:
    if (continuous == "n"):
        thread1 = threading.Thread(target=recordingtemperature1, args=(timeint,countmax,ser,datalist))
        thread2 = threading.Thread(target=readtemperature, args=(newstdin,))
        
        thread1.start()
        thread2.start()
        
        thread1.join()
        thread2.join()

    else:
        thread1 = threading.Thread(target=recordingtemperature2, args=(timeint,countmax,ser,datalist))
        thread2 = threading.Thread(target=readtemperature, args=(newstdin,))
        thread2.daemon = True
            
        thread1.start()
        thread2.start()
    
        thread1.join()
        
    ser.close()





# +-+-+- MAIN PROCESS +-+-+-
# --------------------------
# The main process is solely for plotting the read data.

# -- PATH FOR IMPLEMENTING LATEX (ON OSX) --
os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
os.environ['PATH'] = os.environ['PATH'] + ':/usr/local/texlive/2014'



# -- FUNCTIONS FOR MAIN PROCESS --
def timeconvert(x):
    if (x[len(x)-1] == "s"):
        x = float(x[:len(x)-1])
        return x
    
    if (x[len(x)-1] == "m"):
        x = float(x[:len(x)-1])*60
        return x
    
    if (x[len(x)-1] == "h"):
        x = float(x[:len(x)-1])*3600
        return x
    
    
def cm2inch(value):
    return value/2.54



# -- ACTUAL MAIN PROCESS --
os.system("clear")
print "\n\n\n"
print "-------------------------------------------"
print "--------- TEMPERATURE ACQUISITION ---------"
print "-------------------------------------------"
print "\n\n\n"

timeint = raw_input(" + Time-intervals (#m or #s; at least 0.08s): ")
timeint = timeconvert(timeint)

timemax = raw_input(" + Maximal time plotted (#h, #m or #s): ")
timemax = timeconvert(timemax)

countmax = int(timemax/timeint)

alldata = raw_input(" + Show all data (y/n): ")
if (alldata == "n"):
	plottimemax = raw_input(" + Maximal time shown in plot (#h, #m or #s): ")
	plottimemax = timeconvert(plottimemax)

continuous = raw_input(" + Stop recording after maximal time elapsed (y/n) ?: ")

plotlimit = float(raw_input(" + Plotlimits: "))

# Here, the datamanger is implementedand a proxy list created.
datamanager = multiprocessing.Manager()
datalist = datamanager.list()

# The following two lines are deprecated. Here, a multiprocessing array of fixed size was used.
# The size was determined by 'plotted time'/'time intervals'
#data = multiprocessing.Array(ctypes.c_double, 2*countmax)
#data_np = np.ctypeslib.as_array(data.get_obj()).reshape(2,countmax)

flag = multiprocessing.Value(ctypes.c_int, 0)
settempflag = multiprocessing.Value(ctypes.c_int, 0)
PIDflag = multiprocessing.Value(ctypes.c_int, 0)
Pflag = multiprocessing.Value(ctypes.c_int, 0)
Iflag = multiprocessing.Value(ctypes.c_int, 0)
Dflag = multiprocessing.Value(ctypes.c_int, 0)
OUTPUTflag = multiprocessing.Value(ctypes.c_int, 0)
OUTPUTOFFflag = multiprocessing.Value(ctypes.c_int, 0)
OUTPUTONflag = multiprocessing.Value(ctypes.c_int, 0)

settemp = multiprocessing.Value(ctypes.c_double, 0.5)
Pvalue = multiprocessing.Value(ctypes.c_int, 1)
Ivalue = multiprocessing.Value(ctypes.c_double, 0.5)
Dvalue = multiprocessing.Value(ctypes.c_int, 1)

#if (data_np.nbytes/1e6 > 1000):
#    print "\n + Storage array greater than 1GB!"
#else:
#    print "\n + "+str(data_np.nbytes/1e6)+"MB storage array and "+str(countmax)+" values.\n"

print " + End acquisition with \"exit\":"


newstdin = os.fdopen(os.dup(sys.stdin.fileno()))

# Start the child process:
readingwritingprocess = multiprocessing.Process(target=readingwriting, args=(timeint,countmax,newstdin,continuous,datalist,flag,settemp,Pvalue,Ivalue,Dvalue,settempflag,PIDflag,Pflag,Iflag,Dflag,OUTPUTflag,OUTPUTOFFflag,OUTPUTONflag))
readingwritingprocess.start()


plt.ion()
plt.xlabel("Time [s]")
plt.ylabel("Temperature [C]")
plt.grid(True)
line, = plt.plot([],[])
line2, = plt.plot([],[],"r-")
plotstart = 0
while True:
    try:
        data = datalist[-1]
        if (flag.value == 0):
            try:
                count = np.where(data[0]==0)[0][0]
                if (count == 0): count = np.where(data[0]==0)[0][1]
            except:
                count = countmax
            
            if (alldata == "y"):
                line.set_xdata(data[0][:count])
                line.set_ydata(data[1][:count])
                if (count != 0):
                    plt.xlim(np.amin(data[0][:count])-plotlimit,np.amax(data[0])+plotlimit)
                    plt.ylim(np.amin(data[1][:count])-plotlimit,np.amax(data[1])+plotlimit)
            else:
                if (np.amax(data[0]) >= plottimemax):
                    plotstart = abs(count-plottimemax/timeint)
                else:
                    plotstart = 0
                line.set_xdata(data[0][plotstart:count])
                line.set_ydata(data[1][plotstart:count])
                if (count != 0):
                    plt.xlim(np.amin(data[0][plotstart:count])-plotlimit,np.amax(data[0][plotstart:count])+plotlimit)
                    plt.ylim(np.amin(data[1][plotstart:count])-plotlimit,np.amax(data[1][plotstart:count])+plotlimit)
            #else:
            #line.set_xdata(data_np[0][:count])
            #line.set_ydata(data_np[1][:count])
            #if (count != 0):
            #plt.xlim(np.amin(data_np[0][:count])-plotlimit,np.amax(data_np[0][:count])+plotlimit)
            #plt.ylim(np.amin(data_np[1][:count])-plotlimit,np.amax(data_np[1][:count])+plotlimit)
                    
            line2.set_xdata(np.array([plt.xlim()[0],plt.xlim()[1]]))
            line2.set_ydata(np.array([settemp.value,settemp.value]))
            try:
                temperatureannotation.remove()
                averagetemperatureannotation.remove()
                temperaturestandarddeviationannotation.remove()
            except:
                pass
            
            temperatureannotation = plt.annotate("T = %.3fC" %data[1][count-1], xy=(0.5,0.9), xycoords="axes fraction", ha="center", va="center")
            #if (plotstart != 0):
            averagetemperatureannotation = plt.annotate("T_mean = %.3fC" %np.nanmean(data[1][plotstart:count]), xy=(0.5,0.85), xycoords="axes fraction", ha="center", va="center")
            temperaturestandarddeviationannotation = plt.annotate("DT = %.3fC" %np.nanstd(data[1][plotstart:count]), xy=(0.5,0.8), xycoords="axes fraction", ha="center", va="center")
            #else:
            #averagetemperatureannotation = plt.annotate("T_mean = %.3fC" %np.nanmean(data_np[1][:count]), xy=(0.5,0.85), xycoords="axes fraction", ha="center", va="center")
            #temperaturestandarddeviationannotation = plt.annotate("DT = %.3fC" %np.nanstd(data_np[1][:count]), xy=(0.5,0.8), xycoords="axes fraction", ha="center", va="center")
            plt.draw()
            time.sleep(1)
            
        else:
            break
    except:
        time.sleep(1)


print "\n + All threads ended!\n"

data = datalist[-1]
now = datetime.datetime.now()

storagepath = "/home/labadmin/Documents/Alex/Measurements/Temperature_measurements"
storagefolder = "/%d-%02d-%02d" %(now.year,now.month,now.day)
if (os.path.exists(storagepath+storagefolder) == False):
	os.mkdir(storagepath+storagefolder)
storagepath = storagepath+storagefolder+"/"

if (count != countmax):
    data_copy = np.zeros(shape=(3,count), dtype="float")
    data_copy[0] = np.copy(data[0][:count])
    data_copy[1] = np.copy(data[1][:count])
    data_copy[2] = np.copy(data[2][:count])
    np.savetxt(storagepath+"Temperature_%d_%02d_%02d_%02d_%02d.txt" % (now.year,now.month,now.day,now.hour,now.minute), np.transpose(data_copy), fmt="%.3f", delimiter=",", newline="\n")        
else:
    np.savetxt(storagepath+"Temperature_%d_%02d_%02d_%02d_%02d.txt" % (now.year,now.month,now.day,now.hour,now.minute), np.transpose(data), fmt="%.3f", delimiter=",", newline="\n")        


plt.close()
plt.rc("text", usetex=True)
plt.rc("font", **{"family":"sans-serif","sans-serif":["Helvetica"],"size":11})
plt.rcParams["text.latex.preamble"]=["\\usepackage{siunitx}","\\usepackage[helvet]{sfmath}","\\sisetup{math-rm=\mathsf,text-rm=\sffamily}"]
plt.rcParams["legend.fontsize"]=11

fig = plt.figure(tight_layout=True)
ax = fig.add_subplot(111)
ax.plot(data[0][:count],data[1][:count])
ax.set_xlim(np.amin(data[0][:count])-0.5,np.amax(data[0])+0.5)
ax.set_ylim(np.amin(data[1][:count])-0.5,np.amax(data[1])+0.5)
ax.plot((plt.xlim()[0],plt.xlim()[1]),(settemp.value,settemp.value),"r-")
ax.set_xlabel(r"Time $[\SI{}{\second}]$")
ax.set_ylabel(r"Temperature $[\SI{}{\degreeCelsius}]$")
ax.grid(True)
plt.draw()

if (raw_input(" + Save plot (y/n) ?: ") == "y"):
    fig.set_size_inches(cm2inch(13),cm2inch(8.67))
    plt.savefig(storagepath+"Temperature_%d_%02d_%02d_%02d_%02d.pdf" % (now.year,now.month,now.day,now.hour,now.minute), format="pdf")

plt.close()
        

# -- REMARKS --
# The plotting must be done by the main process.
# The usual raw_input command doesn't work in child processes either. Therefore I implemented a workaround with duplicating the stdin-inputfile.
