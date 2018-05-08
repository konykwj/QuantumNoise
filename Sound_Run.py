from file_io import *
from os import listdir
import sys
from scipy.io.wavfile import write
import numpy as np
import pylab as py


pathname=sys.path[0]

#Modifying to contain either short or long notes only.

#C_Dorian=[440.00,329.63,349.23,0.0,523.25,392.00,493.88,293.66,]
C_Dorian=[329.63,440.00,349.23,0.0,493.88,523.25,293.66,392.00,]



def tone(Pitch, Duration, amp, rate,C123,C231,C312):
    ASDR_on=1
    t = np.linspace(0, Duration, Duration * rate)
    #print Duration
    if ASDR_on==1:
        ASDR=np.zeros(len(t))

        leftover=1.0-(C123+C231+C312)/3.0

        height=(C123+C231+C312)/3.0
        p=0
        q=0
        r=0
        while t[p]<Duration*C123/3.0:
            ASDR[p]=3.0/(Duration*C123)*float(p)/rate
            p=p+1
        while t[p]<Duration*(C123+C231)/3.0:
            ASDR[p]=1-3.0*(1.0-height)/(Duration*C231)*(float(q))/rate
            p=p+1
            q=q+1
        while t[p]<Duration*(C123/3.0+C231/3.0+leftover):
            ASDR[p]=height
            p=p+1
        while t[p]<Duration:
            ASDR[p]=height-3.0*height/(Duration*C312)*(float(r))/rate
            p=p+1
            r=r+1
        if t[p]==Duration:
            ASDR[p]=0.0
        
    else:
        ASDR=np.ones(len(t))
    data = np.sin(2*np.pi*Pitch*t)*amp*ASDR


    
    return data.astype(np.int16) # two byte integers

def amp_calc(M):
    x=4.0/M-1
    Overall_amp=3*10**(1.5*x+2.0)
    return Overall_amp

def pitch_calc(S,S_old,octave_index,octave_list): 
    D_S=abs(S-S_old)
    
    octave=octave_list[octave_index]
    if D_S>.1:
        if octave_index==5:
            octave_index=2
            octave=octave_list[octave_index]
        else:
            octave=octave_list[octave_index+1]
    if D_S<.001:
        octave=octave_list[octave_index-1]
        
    y=(S-4.0)/(4*(np.sqrt(2)-1))
    
    Pitch_index=int(round(y*7))

    if S_old<S:
        if Pitch_index==7:
            Pitch_index=2
        else:
            Pitch_index=Pitch_index+1
        
    if S_old>S:
        Pitch_index=Pitch_index-1
    
    Pitch=C_Dorian[Pitch_index]*octave
    
    
        
    return Pitch,octave_index

def duration_calc(Tgl,dur_control):
    max_dur=1.0
    
    #Short Note Scheme
    if dur_control==1:
        Dur=int(round(Tgl*400,-2))
        
        if Dur==100:
            Duration=max_dur*(3.0/8.0)
        if Dur==200:
            Duration=max_dur*(.25)
        if Dur==400:
            Duration=max_dur*(1.0/16.0)
        if Dur==0:
            Duration=max_dur*(3.0/16.0)
        if Dur==300:
            Duration=max_dur*(1.0/8.0)

    #Longe Note Scheme
    if dur_control==2:
        Dur=int(round(Tgl*300,-2))

        if Dur==300:
            Duration=max_dur*(1.0)
        if Dur==0:
            Duration=max_dur*(.75)
        if Dur==100:
            Duration=max_dur*(.5)
        if Dur==200:
            Duration=max_dur*(.25)
        


    #Triplet note scheme
    if dur_control==3:
        Dur=int(round(Tgl*400,-2))
        max_dur=1

        if Dur==300:
            Duration=max_dur*(1.0/3.0)
        if Dur==100:
            Duration=max_dur*(.25)
        if Dur==0:
            Duration=max_dur*(1.0/8.0)
        if Dur==200:
            Duration=max_dur*(1.0/6.0)
        if Dur==400:
            Duration=max_dur*(1.0/9.0)
            
    return Duration

    
def Create_Music(Data,dur_control):
    octave_index=1
    octave_list=[.25,.5,1.0,2.0,3.0]
    rate=44100
    
    Note_lists=[[],
                [],
                [],
                [],
                [],
                [],
                [],
                [],]

    raw=Data.raw
    
    if raw!=0: #This averages the raw so if I ever have more than one trajectory the program can handle it.
        rawlist=np.zeros(len(raw[0]),dtype=np.ndarray)
        for i in range(len(raw[0])): #For Every timestep
            for j in range(len(raw)): #For Every trajectory
                rawlist[i]=rawlist[i]+raw[j,i]
        rawlist=rawlist/float(len(raw))

            
    for i in range(len(Data.avgSlist)):
        S=Data.avgSlist[i]
        S_old=Data.avgSlist[i-1]
        M=Data.avgMlist[i]
        Tgl=Data.avgtglist[i]
        C123=Data.avgcr123list[i]
        C231=Data.avgcr231list[i]
        C312=Data.avgcr312list[i]
        
        Overall_amp=amp_calc(M)
        Base_Pitch,octave=pitch_calc(S,S_old,octave_index,octave_list)
        Duration=duration_calc(Tgl,dur_control)
        
        for j in range(8):

            if raw==0:
                amp=Overall_amp*1.0/8.0
            else:
                psi=rawlist[i]

                amp=Overall_amp*np.real((psi[j].conjugate()*psi[j]))

                
            Pitch=Base_Pitch*float(j+1)
            
            Note=tone(Pitch, Duration, amp, rate,C123,C231,C312)
            Note_lists[j].append(Note)

    seq1 = np.concatenate((Note_lists[0]),axis=1)
    seq2 = np.concatenate((Note_lists[1]),axis=1)
    seq3 = np.concatenate((Note_lists[2]),axis=1)
    seq4 = np.concatenate((Note_lists[3]),axis=1)
    seq5 = np.concatenate((Note_lists[4]),axis=1)
    seq6 = np.concatenate((Note_lists[5]),axis=1)
    seq7 = np.concatenate((Note_lists[6]),axis=1)
    seq8 = np.concatenate((Note_lists[7]),axis=1)

    song = seq1+seq2+seq3+seq4+seq5+seq6+seq7+seq8

    if dur_control==1:
        dur_name="Short"
    if dur_control==2:
        dur_name="Long"
    if dur_control==3:
        dur_name="Triplet"
        
    write('PHY_Song_test'+str(Data.name)+str(Data.infolist)+dur_name+'.wav', 44100, song)

    return


for m in range(2):
    print m
    dur_control=2+m
    
    numlist=[1,2]
    for u in range(len(numlist)):
        filename="single high driving, W "+str(numlist[u])+"[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 0, 10, 300, 1, 0, 10, 0].txt"
        Data=call_from_file(filename)
        Create_Music(Data,dur_control)
        print filename," done."

    for u in range(len(numlist)):
        filename="single high driving, GHZ "+str(numlist[u])+"[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 0, 10, 300, 1, 0, 10, 0].txt"
        Data=call_from_file(filename)
        Create_Music(Data,dur_control)
        print filename," done."

    for u in range(len(numlist)):
        filename="single low damping, W "+str(numlist[u])+"[1.0, 1.0, 1.0, 0.1, 0.1, 0.1, 1.0, 0, 10, 300, 1, 0, 10, 0].txt"
        Data=call_from_file(filename)
        Create_Music(Data,dur_control)
        print filename," done."
        
    numlist=[1,2]
    for u in range(len(numlist)):
        filename="single low damping, GHZ "+str(numlist[u])+"[1.0, 1.0, 1.0, 0.1, 0.1, 0.1, 1.0, 0, 10, 300, 1, 0, 10, 0].txt"
        Data=call_from_file(filename)
        Create_Music(Data,dur_control)
        print filename," done."

    numlist=[1,2,3]

    for u in range(len(numlist)):
        filename="single high driving, Ground "+str(numlist[u])+"[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 0, 10, 300, 1, 0, 10, 0].txt"
        Data=call_from_file(filename)
        Create_Music(Data,dur_control)
        print filename," done." 

##    filename="single high driving, Ground 11[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 0, 10, 300, 1, 0, 10, 0].txt"
##    Data=call_from_file(filename)
##    dur_control=m+1
##    Create_Music(Data,dur_control)
##    print filename," done."
##
##    filename="Ground, high driving[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 0, 10, 300, 500, 0, 10, 0].txt"
##    Data=call_from_file(filename)
##    Create_Music(Data,dur_control)
##    print filename," done."
##    
##
##    filename="Ground, low damping[1.0, 1.0, 1.0, 0.1, 0.1, 0.1, 1.0, 0, 10, 300, 500, 0, 10, 0].txt"
##    Data=call_from_file(filename)
##    Create_Music(Data,dur_control)
##    print filename," done."
##
##    filename="W, low damping[1.0, 1.0, 1.0, 0.1, 0.1, 0.1, 1.0, 0, 10, 300, 500, 0, 10, 0].txt"
##    Data=call_from_file(filename)
##    Create_Music(Data,dur_control)
##    print filename," done."

    
    filename="W, high driving[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, 0, 10, 300, 500, 0, 10, 0].txt"
    Data=call_from_file(filename)
    Create_Music(Data,dur_control)
    print filename," done."
