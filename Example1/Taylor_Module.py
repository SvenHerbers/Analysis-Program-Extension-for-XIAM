import numpy as np
import subprocess

def binexp(Freqs,binsize,exp_freq,exp_linear_int):
    'This bins the experimental spectrum for you'
    Ints=np.copy(Freqs)*0
    Freqmin=Freqs[0]
    for j in range(0,len(exp_freq)):
        Freq=exp_freq[j]
        idx=int((Freq-Freqmin)/binsize)
        if idx >0 and idx <len(Freqs):
                Ints[idx]+=exp_linear_int[j]
    print("Done binning exp")
    return Ints
def binXIAMoutput(Freqs,binsize,theXI,mask):
    'bins the cat file as loaded from np.loadtxt after formatting'
    Freqmin=Freqs[0]
    Ints=np.copy(Freqs)*0
    theXInew=np.copy(theXI)
    theXInew[:,4]=(10**theXI[:,4])*mask
    for j in range(len(theXI[:,0])):
        Freq=theXInew[:,0][j]
        idx=int((Freq-Freqmin)/binsize)
        if idx >0 and idx <len(Freqs):
                Ints[idx]+=np.float64(theXInew[:,4][j])
    return Ints
def ReadXIAM2NQ(filename):
    'simple function to read out .xo files representing predictions'
    f=open(filename+".xo",'r')
    count=0
    marker=0
    case=0
    for line in f:
        if "d_pop       spcat" in line:
            marker=count
        count += 1
    if marker == 0:
       "ERROR in building mule_xiam"
    f.seek(0)
    g=open("mule_xiam",'w')
    lineno=0
    for line in f:
        if lineno>marker:
            if "." in line :
                    if "str**2" not in line:
                        g.write(line)
                        if ("F1" in line) and ("F " in line):
                            case=2
                        elif ("F" in line):
                            case=1
        lineno+=1
    g.close()
    if case == 0: # no quadrupole coupling
        SFreqintis = np.loadtxt("mule_xiam", usecols=(10, 0, 3, 7, 15, -2, -1))
    if case == 1: # one qudrupole coupling
        SFreqintis = np.loadtxt("mule_xiam", usecols=(13, 0, 3, 10, 18, 7, 8 , -2, -1))
    if case == 2: # two quadrupolar nuclei
        SFreqintis = np.loadtxt("mule_xiam", usecols=(16, 0, 3, 13, 21, 7, 8 ,10, 11, -2, -1))
    SFreqintis[:,0]=SFreqintis[:,0]*1000
    return case, SFreqintis
def ReadXIAM2NQ_list(case,filename):
    'simple function to read out .xo files representing fits'
    f=open(filename,'r')
    count=0
    marker=0
    for line in f:
        if "diff/MHz      obs/GHz" in line:
            marker=count
        count += 1
    if marker == 0:
       "ERROR in building mule_xiam_list"
    f.seek(0)
    g=open("mule_xiam_list",'w')
    lineno=0
    for line in f:
        if lineno>marker:
            if "." in line :
                if "str**2" in line:
                    if "spcat" not in line:
                        g.write(line)
        lineno+=1
    g.close()
    if case == 0:
        Calc_Obs=np.loadtxt("mule_xiam_list",usecols=(9,11,10))
    if case == 1:
        Calc_Obs=np.loadtxt("mule_xiam_list",usecols=(9+3,11+3,10+3))
    if case == 2:
        Calc_Obs=np.loadtxt("mule_xiam_list",usecols=(9+6,11+6,10+6))
    return Calc_Obs
def Run_and_Read_XIAM2NQ(input_name):
    'The routine first runs XIAM once, then reads it'
    input_string="XIAMi2NQ.exe <"+input_name+".xi >"+input_name+".xo"
    subprocess.run(input_string, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    case, output=ReadXIAM2NQ(input_name)
    return case, output
def Calc_taylor_coefsXIAM2NQ(input_name,nrot,linkv):
    'This function derives the coefficients of the 2nd order taylor expansion using repetitive XIAM calls'
    case, theXI=Run_and_Read_XIAM2NQ(input_name)
    list_of_spectra=[]
    list_of_idx_x=[]
    list_of_idx_y=[]
    Inputpar=open(input_name+".xi",'r')
    par_splitlines=[]
    par_lines=[]
    step=1.0/1000 #
    Vsteppercent=0.5 #
    V1n0_1=0.
    V1n0_2=0.
    V1n0_3=0.
    markerV=1
    for line in Inputpar:
        par_splitlines+=[line.split()]
        par_lines+=[line]
        if "BJ" in line:
            BJ=np.float64(line.split()[1])#
        if "BK" in line:
            BK=np.float64(line.split()[1])#
        if "B-" in line:
            Bm=np.float64(line.split()[1])#
        if "V1n" in line:
            if markerV==1:
                V1n0_1 = np.float64(line.split()[1])
                if linkv==1:
                    V1n0_2 = np.float64(line.split()[2])
                markerV += 1
            elif markerV==2:
                V1n0_2 = np.float64(line.split()[1])
                markerV += 1
            elif markerV ==3:
                V1n0_3 = np.float64(line.split()[1])
                markerV += 1
        if markerV>nrot+1 or markerV >4:
            pass
    A0=BJ+BK
    B0=BJ+Bm
    C0=BJ-Bm
    stepv1=(0.1+Vsteppercent/100*V1n0_1)**1.3
    stepv2=(0.1+Vsteppercent/100*V1n0_2)**1.3
    stepv3=(0.1+Vsteppercent/100*V1n0_3)**1.3
    stepv=[stepv1,stepv2,stepv3]
    for x in range(0,4+nrot):
        for y in range (0,4+nrot):
            dir="temp\\"
            label=str(x)+str(y)
            par = open(dir+label+".xi", 'w')
            Aconsti=np.copy(A0)
            if x==1:
                Aconsti += step
            if y==1:
                Aconsti += step
            Bconsti=np.copy(B0)
            if x==2:
                Bconsti += step
            if y==2:
                Bconsti += step
            Cconsti=np.copy(C0)
            if x==3:
                Cconsti += step
            if y==3:
                Cconsti += step
            V1n_1 = np.copy(V1n0_1)
            V1n_2 = np.copy(V1n0_2)
            if x==4:
                V1n_1 += stepv1
                if linkv == 1:
                    V1n_2 += stepv2
            if y==4:
                V1n_1 += stepv1
                if linkv == 1:
                    V1n_2 += stepv2
            if x==5:
                V1n_2 += stepv2
            if y==5:
                V1n_2 += stepv2
            V1n_3 = np.copy(V1n0_3)
            if x==6:
                V1n_3 += stepv3
            if y==6:
                V1n_3 += stepv3
            checkit=0
            BJnew = 0.5 * (Bconsti + Cconsti)
            BKnew = Aconsti - BJnew
            Bmnew = 0.5 * (Bconsti - Cconsti)
            Inputpar.seek(0)
            markerV = 1
            for line in Inputpar:
                if "BJ" in line:
                    par.write(" BJ {:20.15e}\n".format(BJnew))
                elif "BK" in line:
                    par.write(" BK {:20.15e}\n".format(BKnew))
                elif "B-" in line:
                    par.write(" B- {:20.15e}\n".format(Bmnew))
                elif "V1n" in line:
                    if markerV == 1:
                        if linkv == 1:
                            par.write(" V1n {:20.15e} {:20.15e}\n".format(V1n_1,V1n_1))
                        else:
                            par.write(" V1n {:20.15e}\n".format(V1n_1))
                        markerV += 1
                    elif markerV == 2:
                        par.write(" V1n {:20.15e}\n".format(V1n_2))
                        markerV += 1
                    elif markerV == 3:
                        par.write(" V1n {:20.15e}\n".format(V1n_3))
                        markerV += 1
                else:
                    par.write(line)
                if markerV > nrot + 1 or markerV > 4:
                    pass
                if " V " in line:
                    checkit=1
                    par.write("\n")
                    break
            if checkit==0:
                print(" Torsional state definitions not found - this might cause trouble! make a line with spaces around the V like ' V 0 ' or ' V 0 0 '")
            for k in range(0, len(theXI[:,0])):
                J1=theXI[k][1]
                J2=theXI[k][2]
                Freq=theXI[k][0]
                S = theXI[k][3]
                tau1=theXI[k][-2]
                tau2=theXI[k][-1]
                if case == 1:
                    F1=theXI[k][5]
                    F2=theXI[k][6]
                if case == 2:
                    F11=theXI[k][5]
                    F12=theXI[k][6]
                    F1=theXI[k][7]
                    F2=theXI[k][8]
                if tau1%2==0:
                        Ka1=tau1/2
                        Kc1=J1+1-Ka1
                else:
                        Ka1=(tau1-1)/2
                        Kc1=J1-Ka1
                if tau2%2==0:
                        Ka2=tau2/2
                        Kc2=J2+1-Ka2
                else:
                        Ka2=(tau2-1)/2
                        Kc2=J2-Ka2
                if case==0:
                    par.write("{:3.0f}{:3.0f}{:3.0f} {:3.0f}{:3.0f}{:3.0f} S {:3.0f} = {:15.7f}\n".format(J1,Ka1,Kc1,J2,Ka2,Kc2,S,Freq/1000))
                if case==1:
                    par.write("{:3.0f}{:3.0f}{:3.0f} {:3.0f}{:3.0f}{:3.0f} F {:3.0f} {:3.0f} S {:3.0f} = {:15.7f}\n".format(J1,Ka1,Kc1,J2,Ka2,Kc2,F1,F2,S,Freq/1000))
                if case==2:
                    par.write("{:3.0f}{:3.0f}{:3.0f} {:3.0f}{:3.0f}{:3.0f} F1 {:3.0f} {:3.0f} F {:3.0f} {:3.0f} S {:3.0f} = {:15.7f}\n".format(J1,Ka1,Kc1,J2,Ka2,Kc2,F11,F12,F1,F2,S,Freq/1000))
            par.close()
            xiamstring="XIAMi2NQ.exe <"+dir+label+".xi >"+dir+label+".xo"
            subprocess.run(xiamstring, shell=True, check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.STDOUT)
            list_of_spectra+=[ReadXIAM2NQ_list(case,dir+label+".xo")[:,0]]
            list_of_idx_x+=[x]
            list_of_idx_y+=[y]
            if x==0:  # only 00 is required as reference.
                break
    #### calculate taylor
    f0=list_of_spectra[0]
    fd=first_derivatives=[]
    fd_idx=[]
    for j in range(1,len(list_of_idx_x)):
        x=list_of_idx_x[j]
        y=list_of_idx_y[j]
        if y==0:
            if x < 4:
                fd+=[(list_of_spectra[j]-list_of_spectra[0])/step]
            else:
                fd+=[(list_of_spectra[j]-list_of_spectra[0])/(stepv[x-4])]
            fd_idx+=[j]
    sd=second_derivative=[]
    for j in range(1,len(list_of_idx_x)):
        x=list_of_idx_x[j]
        y=list_of_idx_y[j]
        if y==0:
         pass
        else:
            S_xy=list_of_spectra[j]
            S_x0=list_of_spectra[fd_idx[x-1]]
            S_y0=fd[y-1]
            if y < 4:
                grad1=(S_xy - S_x0) / step
            else:
                grad1 = (S_xy - S_x0) / (stepv[y-4])
            grad2=S_y0
            if x < 4:
                sd += [(grad1-grad2) / step]
            else:
                sd += [(grad1 - grad2)  / (stepv[x-4])]
    Rot0=np.array([A0,B0,C0,V1n0_1,V1n0_2,V1n0_3])
    Inputpar.close()
    return Rot0,f0,fd,sd,theXI
def TaylorXIAM(Rot0,Rot,f0,fd,sd,nrot):
    'This function calculates the values according to the 2nd order taylor expansion'
    [A0,B0,C0,V1n0_1,V1n0_2,V1n0_3]=Rot0
    [A,B,C,V1n_1,V1n_2,V1n_3]=Rot
    a=A-A0
    b=B-B0
    c=C-C0
    v1=V1n_1-V1n0_1
    v2=V1n_2-V1n0_2
    v3=V1n_3-V1n0_3
    f0=np.array(f0)
    if nrot==0:
        xDf=fd[0]*a+fd[1]*b+fd[2]*c
    if nrot==1:
        xDf=fd[0]*a+fd[1]*b+fd[2]*c+fd[3] * v1
    if nrot==2:
        xDf=fd[0]*a+fd[1]*b+fd[2]*c+fd[3] * v1 + fd[4] * v2
    if nrot==3:
        xDf=fd[0]*a+fd[1]*b+fd[2]*c+fd[3] * v1 + fd[4] * v2 + fd[5] * v3
    Sec=(3+nrot)*[0]
    for k in range(len(Sec)):
        if nrot==0:
            Sec[k] = sd[0 + k * (3 + nrot)] * a + sd[1 + k * (3 + nrot)] * b + sd[2 + k * (3 + nrot)] * c
        if nrot == 1:
            Sec[k] = sd[0 + k * (3 + nrot)] * a + sd[1 + k * (3 + nrot)] * b + sd[2 + k * (3 + nrot)] * c + sd[3
                        + k * (3 + nrot)] * v1
        if nrot == 2:
            Sec[k] = sd[0 + k * (3 + nrot)] * a + sd[1 + k * (3 + nrot)] * b + sd[2 + k * (3 + nrot)] * c + sd[3
                        + k * (3 + nrot)] * v1 + sd[4 + k * (3 + nrot)] * v2
        if nrot == 3:
            Sec[k] = sd[0 + k * (3 + nrot)] * a + sd[1 + k * (3 + nrot)] * b + sd[2 + k * (3 + nrot)] * c + sd[3
                        + k * (3 + nrot)] * v1 + sd[4 + k * (3 + nrot)] * v2 + sd[3 + k * (3 + nrot)] * v1 + sd[5 + k * (3 + nrot)] * v3
    if nrot==0:
        x2D2f=Sec[0]*a+Sec[1]*b+Sec[2]*c
    if nrot==1:
        x2D2f=Sec[0]*a+Sec[1]*b+Sec[2]*c+Sec[3]*v1
    if nrot==2:
        x2D2f=Sec[0]*a+Sec[1]*b+Sec[2]*c+Sec[3]*v1+Sec[4]*v2
    if nrot==3:
        x2D2f=Sec[0]*a+Sec[1]*b+Sec[2]*c+Sec[3]*v1+Sec[4]*v2+Sec[5]*v3
    T = np.array(f0) + xDf + 1/2 * x2D2f # zeros order
    return T


