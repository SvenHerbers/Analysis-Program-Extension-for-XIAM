# Analysis-Program-Extension-for-XIAM (APEX)
**APEX** is a Python script that uses `matplotlib` and `matplotlib.widgets` to interactively visualize how changes in rotational (and optionally torsional) parameters affect a predicted rotational spectrum. With APEX, you can adjust parameters in real time using sliders and immediately compare the calculated spectrum with experimental data, making spectral assignment faster and more intuitive.

Intended use: The use of APEX is to find very first inital assignments for a rotational spectrum. If you have an experimental spectrum and a computational prediction to compare it with, you can plot them together using APEX. It is then possible to move patterns of lines by continously varying rotational constants or internal rotation barriers, which will help Identify matching groups of transitions.
If a good matching as achived, APEX can be closed and the new parameters can be used as a refined intial guess to carry out the fitting procedure.

This README provides the required information to run APEX, along with a short tutorial based on two example datasets.

## Required Files

All of the following files must be located in the **same folder**:

| File | Description |
|------|-------------|
| Experimental spectrum file (e.g., `Spectrum.txt`) | Two-column text file: frequency (MHz), intensity (arb. units); plotted in the background. |
| `Input.xi` | XIAM input file used by `xiam2nq`. |
| `APEX.py` | Main program — run this file to start APEX. |
| `Taylor_Module.py` | Contains helper functions used by `APEX.py`. |
| `XIAMi2NQ.exe` | Backend executable performing the XIAM/xiam2nq calculation. |

## Documentation and How-to-use



### Example 1 - Synthetic two Rotor case

Several steps have to be followed to get APEX running. 

#### Step 1 - Set up Input.xi
In a first step, a reference prediction input file called `Input.xi` has to be provided.
The structure of `Input.xi` follows the general form used in XIAM, however, there are a few additional restrictions. On the one hand comments should be avoided, since APEX checks for strings. If comments contain strings by accident, this could cause malfunction. On the other hand Barrier `V1n` must be provided in seperate lines if the barriers are not linked. Or they must be provided in the same line, if the two barrier are to be linked to one each other (equivalent top case).

The following input file shows a typical two rotor case with no quadrupole coupling

```
~Synthetic Testmolecule. C1, near prolate
 
ntop  2   
freq  4.5 7.5
int 2
limit -6  
temp  6         
qsum  3670.1  
reduc 1 
print 4 
                        
 BJ        0.6640387722137
 BK        1.9618593342137
 B-        0.0234845941710
 DJ            0.237760E-6  
 DJK          -2.156818E-6  
 DK            6.147716E-6  
 dj           -0.010541E-6  
 dk           -0.001420E-6  
 V1n          11162.94423  
 V1n           5241.09574
 F0          161.024267987  160.249804221 
 epsil         0.936243809    2.748092877 
 delta         2.651195789    0.750018747 
  mu_x           2.1417731                        
  mu_y           2.2346126             
  mu_z           0.3338271              
                      
                      
   S     0   0   !S 1 
   S     1   0   !S 2 
   S     0   1   !S 3 
   S     1   1   !S 4 
   S    -1   1   !S 5 
                      
 V 0 0            ! ground state    

J 15

```
The parameters have the following meaning 
| Parameter | Meaning / Description |
|----------|------------------------|
| `ntop` | Number of internal rotors in the molecule. |
| `freq` | Frequency prediction window (lower / upper limit in GHz). |
| `int` | Intensity calculation setting (type of intensity treatment). |
| `limit` | Logarithmic intensity limit referring to `log10_total` in prediction outputs. |
| `temp` | Rotational temperature (K) used for intensity weighting. |
| `qsum` | Partition function value used for SPCAT type intensities. |
| `reduc` | Type of reduction for centrifugal distortion (0: Watson A, 1: Watson S) |
| `print` | Output verbosity (higher = more printed internal info). |
| `BJ` |`0.5 (B_x + B_y)`|
| `BK` |`B_z - 0.5 (B_x + B_y)`|
| `B-` |`0.5 (B_x - B_y)`|
| `DJ`, `DJK`, `DK`, `dj`, `dk` | Centrifugal distortion constants (quartic). |
| `V1n` | Barrier height of torsional potential for rotor 1 and rotor 2 (in GHz). |
| `F0` | Torsional rotational constant (one value per rotor). |
| `epsil` | The angle (radian) between the principal axis x and the projection of the internal rotation axis onto xy-plane. |
| `delta` | The angle (radian) between the internal rotation axis and the principal axis z. |
| `mu_x`, `mu_y`, `mu_z` | Dipole moment components (Debye). |
| `S 0 0`, `S 1 0`, … | Defines torsional symmetry species for internal rotors. |
| `V 0 0` | Specifies vibrational state (v = 0 ground state here). |
| `J 15` | Maximum rotational quantum number used in the calculation. |

##### Notes and Things to Consider

- The identification of the molecular axes `x`, `y`, `z` with the principal axes `a`, `b`, `c` determines the chosen representation.  
  Example: `z = a`, `x = b`, `y = c` corresponds to the **Ir representation**.

- Centrifugal distortion constants (and other additional fixed parameters) are **not controlled by sliders in APEX**.  
  From experience predicted centrifugal distortion coefficients using hybrid DFT methods are fairly accurate and should also be included in initial predictions if available.

- The `qsum` parameter corresponds to the rotational partition function.  
  It mostly affects the **absolute intensity scale**, not the **relative intensities**.  
  You may:
  - set `qsum` to an estimate, or
  - simply assign `qsum  1000`, or
  - omit the line entirely (XIAM defaults to **1000**).

  Example function to compute `qsum` (σ = 1 for C₁ / Cₛ symmetry):

```python
import numpy as np
import scipy.constants as const
def calc_qr(A,B,C,sigma,T):
    "Rotational constants in MHz, symmetry number and Temperatur are used to calculate qr"
    A=A*10**6*const.Planck/const.k
    B=B*10**6*const.Planck/const.k
    C=C*10**6*const.Planck/const.k
    qr=(np.pi**(1/2)/sigma) * np.sqrt((T**3/(A*B*C)))
    return qr
```

- APEX requires the spectrum to be calculated using `int 2`.  
  **No other intensity mode is supported.**

- For torsional barriers (`V1n`), independent rotors must be written on **separate lines**, allowing APEX to detect non-equivalent, unlinked rotors correctly.

- Blank lines in the `.xi` input **matter**.  
  XIAM interprets empty lines as block separators — removing them will break the file structure and XIAM will fail to run.

Run the following command in a terminal (cmd):

```bash
XIAMi2NQ.exe < Input.xi > Input.xo
```

If everything is set up correctly, the output (`Input.xo`) should begin similar to this:
```
XIAMi2NQ.exe < Input.xi > Input.xo
```
The output should look sth like this:

```
 Rotational, Centrifugal Distortion, Internal Rotation Calculation (V2.5e)
                       Holger Hartwig 08-Nov-96 (hartwig@phc.uni-kiel.de)

 Please cite: H.Hartwig and H.Dreizler, Z.Naturforsch, 51a (1996) 923.

 
Modified Version: XIAM-2NQ v0.24b -By Sven Herbers 17-June-2025

...
...
...

  spcat is intensity as expected from spcat
-- B 1                                Freq log10_str**2 log10_total     stat.w.       d_pop       spcat
  2  1  2   1  0  1   S 1  V 1     4.549426      0.8370     -3.7020      3.0000     -4.5391     -4.9019  B 1  K -1  0  t  2  1
  2  1  1   1  0  1   S 1  V 1     4.690542      0.8750     -3.6510      3.0000     -4.5260     -4.8376  B 1  K  1  0  t  3  1
  2  2  1   2  1  2   S 1  V 1     5.961392      0.5818     -3.6363      5.0000     -4.2181     -4.9406  B 1  K -2 -1  t  4  2
  2  2  0   2  1  2   S 1  V 1     5.962236      0.6040     -3.6140      5.0000     -4.2180     -4.9183  B 1  K  2 -1  t  5  2
  2  2  1   2  1  1   S 1  V 1     5.820275      0.6198     -3.6090      5.0000     -4.2287     -4.9236  B 1  K -2  1  t  4  3
  2  2  0   2  1  1   S 1  V 1     5.821120      0.5972     -3.6315      5.0000     -4.2287     -4.9461  B 1  K  2  1  t  5  3
...
...
...
```


If **no intensity lines** are printed (lines containing `log10_str**2` / `spcat` values), something is incorrect in the input.

If **too many** lines appear (especially many very weak transitions):

- Increase the value of the `limit` parameter (this removes weak transitions).
- Keep the `freq` range as narrow as possible.

The more lines there are in the prediction file the longer APEX will take to set up.

#### Step 2 - Set up APEX.py

`Apex.py` has an input section starting at line 15 and ending at line 36. The following parameters are available:
```
#========================================
#===         APEX INPUT SECTION       ===
#========================================
experiment      = np.loadtxt("synthetic_spectrum.txt")  # load experimental spectrum (2 columns: freq/intensity)
figureaspects   = (16,8)     # Change the intial plotsize to match your screen size, if it does not fit (you can also resize after plotting)

Variable_choice = 1          # slider parameters: 0 = A,B,C ; 1 = BJ,BK,B− (Ir assumption)
Arange          = 30.0       # slider range for A or BJ (MHz)
Brange          = 30.0       # slider range for B or BK (MHz)
Crange          = 10.0        # slider range for C or B− (MHz)
v1range         = 30.0       # torsional barrier range V1 (cm^-1)
v2range         = 30.0       # torsional barrier range V2 (cm^-1)
v3range         = 30.0       # torsional barrier range V3 (cm^-1)
freqlo          = 4500.0     # lower frequency limit to plot (MHz)
freqhigh        = 7500.0     # upper frequency limit to plot (MHz)
binsize         = 0.15       # bin width for experimental and prediction spectrum in plot (MHz)
nrot            = 2          # number of internal rotors that should get a slider (can be less than number of rotors in prediction file)
speciesplot     = [0,1,1,1,1,1]   # species selection for plot: index 0 = combined, 1–5 = individual S1,S2,S3,S4,S5
link_v          = 0          # link torsional potentials v3_1 and v3_2 (requires nrot == 1)

#========================================
#== End of Input, Starting Computation ==
#========================================
```
##### Notes and Things to Consider
`experiment` can be either a complete spectrum or a line list with line intensities. Since rebinning is carried out, both will work the same way.
It is up to the user, how to load experiment as a 2D array into the program. I usually use the function np.loadtxt, which allows also to specify `delimiter` if there are special delimiters between values (default is space) and `skiprows` if there is a header.
`binsize` should be chosen rather large. Since the purpose of this program is to carry out initial assignments for better starting points.
`freqlo` and `freqhigh` can be set to the same limits as in the `Input.xi`, however, here they are defined in MHz and not GHz then. They can also be chosen to be more constrained, then part of the predictions are ignored.
`speciesplot` is just a mask, for which things to show in the plot. `speciesmask[1:5]` represents symmetryspeices S1,S2,S3,S4,S5 which can each be plot individually this way. `speciesmask[0]` represents the combined spectrum, which also works if the number of species exceeds 5.
`link_v` is 0 for non-equivalent rotors, and 1 for two equivalent rotors. If set to 1, both barriers should be defined in one line in `Input.xi` (e.g. `V1n 1162.94423 1162.94423`)

#### Step 3 - Run APEX.py

# Citation
If you use **APEX** in a publication or presentation, please cite this repository and/or the following ISMS conference abstract:

Sven Herbers and Ha Vinh Lam Nguyen, *A New Spectral Assignment Tool and Rotational Characterization of 3,4-Dihydro-2H-Pyran-2-Carboxaldehyde (C₆H₈O₂).* 77th International Symposium on Molecular Spectroscopy (ISMS), 19 June 2024, University of Illinois at Urbana–Champaign, Urbana, IL, USA.

See ISMS2024_Abstract.png for the conference abstract in this repository


