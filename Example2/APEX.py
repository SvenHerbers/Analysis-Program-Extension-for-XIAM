import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import Taylor_Module as Taylor
import os
import glob
print("Version published 01.11.2025")
print("Script to allow for spectral shifting based on 2nd order Taylor expansion.")
print("Errors of less than 4MHz for 100 MHz shifts can be expected.")
print("APEX_output.txt will contain the last slider values after closing the plot.")
try:
    os.mkdir("temp")
except:
    pass
#========================================
#===         APEX INPUT SECTION       ===
#========================================
experiment      = np.loadtxt("Spectrum_3Bar_He_Ar_RT_800KAVG.txt",skiprows=15,delimiter=",")  # load experimental spectrum (2 columns: freq/intensity)
figureaspects   = (16,8)     # Change the intial plotsize to match your screen size, if it does not fit (you can also resize after plotting)

Variable_choice = 1          # slider parameters: 0 = A,B,C ; 1 = BJ,BK,B− (Ir assumption)
Arange          = 30.0       # slider range for A or BJ (MHz)
Brange          = 30.0       # slider range for B or BK (MHz)
Crange          = 5.0        # slider range for C or B− (MHz)
v1range         = 30.0       # torsional barrier range V1 (cm^-1)
v2range         = 30.0       # torsional barrier range V2 (cm^-1)
v3range         = 30.0       # torsional barrier range V3 (cm^-1)
freqlow         = 2000.0     # lower frequency limit to plot (MHz)
freqhigh        = 8000.0     # upper frequency limit to plot (MHz)
binsize         = 2.0        # bin width for experimental and prediction spectrum in plot (MHz)
nrot            = 0          # number of internal rotors that should get a slider (can be less than number of rotors in prediction file)
speciesplot     = [1,0,0,0,0,0]   # species selection for plot: index 0 = combined, 1–5 = individual S1,S2,S3,S4,S5
link_v          = 0          # link torsional potentials v3_1 and v3_2 (requires nrot == 1)

#========================================
#== End of Input, Starting Computation ==
#========================================


sp=speciesplot
Freqs=np.arange(freqlow,freqhigh,binsize,dtype=np.float64) # For binning
Inputfile="Input"#without extension
Rot0,f0,fd,sd,theXI=Taylor.Calc_taylor_coefsXIAM2NQ(Inputfile,nrot, link_v)
Conversion = C = np.array([1000,1000,1000,1/29.9792458,1/29.9792458,1/29.9792458]) # Convert XIAMs GHz to MHz and wavenumbers
Rot0=np.array(Rot0)*C
f0=np.array(f0)
fd=np.array(fd)
sd=np.array(sd)
BJ0=0.5*(Rot0[1]+Rot0[2])
BK0=Rot0[0]-0.5*(Rot0[1]+Rot0[2])
Bm0=0.5*(Rot0[1]-Rot0[2])
BJrot0=[BJ0,BK0,Bm0,Rot0[3],Rot0[4],Rot0[5]]
Rot=np.copy(Rot0)
T0= Taylor.TaylorXIAM(Rot0/C,Rot/C,f0,fd,sd,nrot)
T= np.copy(T0)

mask0=1                                           *sp[0]
mask1=np.array([theXI[:,3]==1.0],dtype=np.float64)*sp[1] # the various symmetry species - changing the last multiplyer from 1 to 0 removes the respective symmetry species
mask2=np.array([theXI[:,3]==2.0],dtype=np.float64)*sp[2]
mask3=np.array([theXI[:,3]==3.0],dtype=np.float64)*sp[3]
mask4=np.array([theXI[:,3]==4.0],dtype=np.float64)*sp[4]
mask5=np.array([theXI[:,3]==5.0],dtype=np.float64)*sp[5]
Scolors=["tab:blue","tab:red","tab:orange","tab:purple","tab:green"]
Parameters = [Rot0, f0, fd, sd, theXI]
#print(T)
def update(val):
    [Rot0, f0, fd, sd, theXI]=Parameters
    if Variable_choice ==  0:
        Rot[:]=np.array([freq_slider1.val,freq_slider2.val,freq_slider3.val,freq_slider4.val,freq_slider5.val,freq_slider6.val])
    else:
        Rot[:]=np.array([freq_slider2.val+freq_slider1.val,freq_slider1.val+freq_slider3.val,freq_slider1.val-freq_slider3.val,freq_slider4.val,freq_slider5.val,freq_slider6.val])
    T = Taylor.TaylorXIAM(Rot0/C, (Rot)/C, f0, fd, sd, nrot)
    newcat=np.copy(theXI)
    newcat[:,0]=T*1000
    uod=button2.val#upordown
    line0.set_ydata(uod*amp_slider.val * Taylor.binXIAMoutput(Freqs, binsize, newcat, mask0))
    line1.set_ydata(uod*amp_slider.val * Taylor.binXIAMoutput(Freqs, binsize, newcat, mask1))
    line2.set_ydata(uod*amp_slider.val * Taylor.binXIAMoutput(Freqs, binsize, newcat, mask2))
    line3.set_ydata(uod*amp_slider.val * Taylor.binXIAMoutput(Freqs, binsize, newcat, mask3))
    line4.set_ydata(uod*amp_slider.val * Taylor.binXIAMoutput(Freqs, binsize, newcat, mask4))
    line5.set_ydata(uod*amp_slider.val * Taylor.binXIAMoutput(Freqs, binsize, newcat, mask5))
    fig.canvas.draw_idle()
    ax.set_ylim(ax.get_ylim()[1]*0.5*(uod-1),ax.get_ylim()[1])
def refreshslider(slider,srange):
    slider.valinit=slider.val
    slider.valmin=slider.valinit-srange,
    slider.valmax=slider.valinit+srange,
    slider.ax.set_xlim(slider.valmin,slider.valmax)
def reset(event):
    freq_slider1.reset()
    freq_slider2.reset()
    freq_slider3.reset()
    freq_slider4.reset()
    freq_slider5.reset()
    freq_slider6.reset()
    amp_slider.reset()
def mirror(event):
    x=button2.val
    button2.set_val(-x)
    if x == 1:
        button3.color = "tab:green"
    else:
        button3.color = "0.85"

def zoom_factory(ax,base_scale = 2.): # Stackoverflow user tacaswell - 2012 modified #NOT USE AT THE MOMENT
    def zoom_fun(event):
        # get the current x and y limits
        cur_xlim = ax.get_xlim()
        cur_ylim = ax.get_ylim()
        cur_xrange = (cur_xlim[1] - cur_xlim[0])*.5
        cur_xcenter = (cur_xlim[1] + cur_xlim[0])*.5
        #cur_yrange = (cur_ylim[1] - cur_ylim[0])*.5
        xlimits = ax.get_xlim()
        #xdata = (xlimits[1]-xlimits[1])/2
        xdata = event.xdata # get event x location
        ydata = event.ydata # get event y location
        if event.button == 1:
            # deal with zoom in
            scale_factor = 1/base_scale
        elif event.button == 3:
            # deal with zoom out
            scale_factor = base_scale
        else:
            # deal with something that should never happen
            scale_factor = 1
            print(event.button)
        # set new limits
        if xdata != None: # Only zoom when clicking on graph
         if ydata != None:
             if ax.get_xlabel() == 'Frequency / MHz':
                ax.set_xlim([cur_xcenter - cur_xrange*scale_factor,
                     cur_xcenter + cur_xrange*scale_factor])
        plt.draw() # force re-draw

    fig = ax.get_figure() # get the figure of interest
    # attach the call back
    fig.canvas.mpl_connect('button_press_event',zoom_fun)
    return zoom_fun

def snap_zoom(event):
    uod = button2.val
    ax.set_ylim(ax.get_ylim()[1] * 0.5 * (uod - 1), ax.get_ylim()[1])
    plt.draw() # force re-draw
# Create the figure and the line that we will manipulate
fig, ax = plt.subplots(figsize=figureaspects)
f=Taylor.binXIAMoutput(Freqs,binsize,theXI,np.ones(len(theXI[:,0])))

expintbin=Taylor.binexp(Freqs,binsize,experiment[:,0],experiment[:,1])
if np.max(f) == 0:
    print("ERROR - XIAM - Prediction gave zero intensities")
    ax.plot(Freqs,expintbin/np.max(expintbin), lw=0.5,color='k')
else:
    ax.plot(Freqs,expintbin/np.max(expintbin)*np.max(f), lw=0.5,color='k')
#
line0, = ax.plot(Freqs, Taylor.binXIAMoutput(Freqs,binsize,theXI,mask0), color=Scolors[0],lw=1,alpha=0.5)
line1, = ax.plot(Freqs, Taylor.binXIAMoutput(Freqs,binsize,theXI,mask1), color=Scolors[0],lw=1,alpha=0.5)
line2, = ax.plot(Freqs, Taylor.binXIAMoutput(Freqs,binsize,theXI,mask2), color=Scolors[1],lw=1,alpha=0.5)
line3, = ax.plot(Freqs, Taylor.binXIAMoutput(Freqs,binsize,theXI,mask3), color=Scolors[2],lw=1,alpha=0.5)
line4, = ax.plot(Freqs, Taylor.binXIAMoutput(Freqs,binsize,theXI,mask4), color=Scolors[3],lw=1,alpha=0.5)
line5, = ax.plot(Freqs, Taylor.binXIAMoutput(Freqs,binsize,theXI,mask5), color=Scolors[4],lw=1,alpha=0.5)
ax.set_xlabel('Frequency / MHz')

# adjust the main plot to make room for the sliders

fig.subplots_adjust(left=0.25, bottom=0.25+nrot*0.05)
xslidecoords=xsc=np.array([0.25,0.25,0.25,0.25,0.25,0.25])
yslidecoords=ysc=np.array([0.15,0.1,0.05,-11,-12,-13])
xslidewidth=xsw=np.array([0.65,0.65,0.65,0.65,0.65,0.65])
yslidewidth=ysw=np.array([0.03,0.03,0.03,0.03,0.03,0.03])
if nrot==0:
    pass
elif nrot==1:
    ysc+=0.05
    ysc[3]=0.05
elif nrot==2:
    ysc += 0.1
    ysc[3] = 0.1
    ysc[4] = 0.05
elif nrot==3:
    ysc += 0.15
    ysc[3] = 0.15
    ysc[4] = 0.1
    ysc[5] = 0.05
elif nrot>3:
    print("Error - your nrot parameter is too high - only up to three rotors are accepted")
# Make a horizontal slider to control the frequency.

if Variable_choice==0:
    s = slidervalues = Rot0
    sl=sliderlabels=[r'$A/\mathrm{MHz}$',r'$B/\mathrm{MHz}$',r'$C/\mathrm{MHz}$',r'$V_{3,1}/\mathrm{cm^{-1}}$',r'$V_{3,2}/\mathrm{cm^{-1}}$',r'$V_{3,3}/\mathrm{cm^{-1}}$']
else:
    s=slidervalues=BJrot0
    sl = sliderlabels = [r'$B_J/\mathrm{MHz}$',r'$B_K/\mathrm{MHz}$',r'$B_{-}/\mathrm{MHz}$',r'$V_{3,1}/\mathrm{cm^{-1}}$',r'$V_{3,2}/\mathrm{cm^{-1}}$',r'$V_{3,3}/\mathrm{cm^{-1}}$']
m=0
axfreq = fig.add_axes([xsc[m], ysc[m], xsw[m], ysw[m]])
freq_slider1 = Slider(
    ax=axfreq,
    label=sl[m],
    valmin=s[m]-Arange,
    valmax=s[m]+Arange,
    valinit=s[m],
    valfmt="%.4f",
)
m+=1
axfreq = fig.add_axes([xsc[m], ysc[m], xsw[m], ysw[m]])
freq_slider2 = Slider(
    ax=axfreq,
    label=sl[m],
    valmin=s[m]-Brange,
    valmax=s[m]+Brange,
    valinit=s[m],
    valfmt="%.4f",
)
m+=1
axfreq = fig.add_axes([xsc[m], ysc[m], xsw[m], ysw[m]])
freq_slider3 = Slider(
    ax=axfreq,
    label=sl[m],
    valmin=s[m]-Crange,
    valmax=s[m]+Crange,
    valinit=s[m],
    valfmt="%.4f",
)
m+=1
axfreq = fig.add_axes([xsc[m], ysc[m], xsw[m], ysw[m]])
freq_slider4 = Slider(
    ax=axfreq,
    label=sl[m],
    valmin=s[m]-v1range,
    valmax=s[m]+v1range,
    valinit=s[m],
    valfmt="%.4f",
)
m+=1
if link_v==1:
    axfreq = fig.add_axes([xsc[m], ysc[m]-10., xsw[m], ysw[m]])
else:
    axfreq = fig.add_axes([xsc[m], ysc[m], xsw[m], ysw[m]])


freq_slider5 = Slider(
    ax=axfreq,
    label=sl[m],
    valmin=s[m]-v2range,
    valmax=s[m]+v2range,
    valinit=s[m],
    valfmt="%.4f",
)
m+=1
axfreq = fig.add_axes([xsc[m], ysc[m], xsw[m], ysw[m]])
freq_slider6 = Slider(
    ax=axfreq,
    label=sl[m],
    valmin=s[m]-v3range,
    valmax=s[m]+v3range,
    valinit=s[m],
    valfmt="%.4f",
)
# Make a vertically oriented slider to control the amplitude
axamp = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
amp_slider = Slider(
    ax=axamp,
    label="Amplitude",
    valmin=0,
    valmax=5,
    valinit=1.0,
    orientation="vertical"
)

# register the update function with each slider
freq_slider1.on_changed(update)
freq_slider2.on_changed(update)
freq_slider3.on_changed(update)
freq_slider4.on_changed(update)
freq_slider5.on_changed(update)
freq_slider6.on_changed(update)
amp_slider.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = fig.add_axes([0.85, 0.84, 0.05, 0.04])
mirrorax = fig.add_axes([0.85, 0.80, 0.05, 0.04])
mirroraxb = fig.add_axes([0.85, -10.80, 0.05, 0.01])
#retaylorax = fig.add_axes([0.80, 0.84, 0.05, 0.04])

button2 = Slider( # Fake button to store value of 1 or -1.0 to pass to update when mirroring
    ax=mirroraxb,
    label="",
    valmin=0,
    valmax=1.0,
    valinit=1.0,
    orientation="horizontal"
)
button2.valtext.set_visible(False)
button3 = Button(mirrorax, 'Mirror', hovercolor='0.975')
button = Button(resetax, 'Reset', hovercolor='0.975')
button.on_clicked(reset)
button3.on_clicked(mirror)
button2.on_changed(update)
#button2.on_clicked(mirror)
ax.set_ylim(0,ax.get_ylim()[1])
scale = 1.5
#f = snap_zoom('button_release_event')
fig = ax.get_figure()  # get the figure of interest
fig.canvas.mpl_connect('button_release_event', snap_zoom)
#f = zoom_factory(ax, base_scale=scale)
plt.show()
g=open("APEX_output.txt",'w')
if Variable_choice           ==   1 :
    g.write("A    / MHz  {:25.10f}\n".format(freq_slider2.val+freq_slider1.val))
    g.write("B    / MHz  {:25.10f}\n".format(freq_slider1.val+freq_slider3.val))
    g.write("C    / MHz  {:25.10f}\n".format(freq_slider1.val-freq_slider3.val))
    g.write("BJ   / GHz  {:25.13f}\n".format(freq_slider1.val/1000))
    g.write("BK   / GHz  {:25.13f}\n".format(freq_slider2.val/1000))
    g.write("B-   / GHz  {:25.13f}\n".format(freq_slider3.val/1000))
else:
    g.write("A    / MHz  {:25.10f}\n".format(freq_slider1.val))
    g.write("B    / MHz  {:25.10f}\n".format(freq_slider2.val))
    g.write("C    / MHz  {:25.10f}\n".format(freq_slider3.val))
    g.write("BJ   / GHz  {:25.13f}\n".format(0.5*(freq_slider2.val+freq_slider3.val)/1000))
    g.write("BK   / GHz  {:25.13f}\n".format((freq_slider1.val-0.5*(freq_slider2.val+freq_slider3.val))/1000))
    g.write("B-   / GHz  {:25.13f}\n".format(0.5*(freq_slider2.val-freq_slider3.val)/1000))

if nrot == 1:
    g.write("V3   / GHz  {:17.5f}             V3   / cm-1   {:25.10f}\n".format(freq_slider4.val*29.9792458,freq_slider4.val))
if nrot == 2:
    g.write("V3_1 / GHz  {:17.5f}             V3_1 / cm-1   {:25.10f}\n".format(freq_slider4.val*29.9792458,freq_slider4.val))
    g.write("V3_2 / GHz  {:17.5f}             V3_2 / cm-1   {:25.10f}\n".format(freq_slider5.val*29.9792458,freq_slider5.val))
if nrot == 3:
    g.write("V3_1 / GHz  {:17.5f}             V3_1 / cm-1   {:25.10f}\n".format(freq_slider4.val*29.9792458,freq_slider4.val))
    g.write("V3_2 / GHz  {:17.5f}             V3_2 / cm-1   {:25.10f}\n".format(freq_slider5.val*29.9792458,freq_slider5.val))
    g.write("V3_3 / GHz  {:17.5f}             V3_3 / cm-1   {:25.10f}\n".format(freq_slider6.val*29.9792458,freq_slider6.val))
g.close()
cleanup=1
templist=["mule_xiam","mule_xiam_list","tori.b1","recalc.xi","recalc.xo"]
tempxi=glob.glob("temp\\*.xi")
tempxo=glob.glob("temp\\*.xo")
temptori=["temp\\tori.b1"]


if cleanup==1:
    for file in templist:
        try:
            os.remove(file)
        except:
            pass
    for file in tempxi:
        try:
            os.remove(file)
        except:
            pass
    for file in tempxo:
        try:
            os.remove(file)
        except:
            pass
    for file in temptori:
        try:
            os.remove(file)
        except:
            pass
    try:
        os.rmdir("temp")
    except:
        pass


    try:
        os.remove("mule_xiam_list")
    except:
        pass

