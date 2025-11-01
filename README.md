# Analysis-Program-Extension-for-XIAM (APEX)
**APEX** is a Python script that uses `matplotlib` and `matplotlib.widgets` to interactively visualize how changes in rotational (and optionally torsional) parameters affect a predicted rotational spectrum. With APEX, you can adjust parameters in real time using sliders and immediately compare the calculated spectrum with experimental data, making spectral assignment faster and more intuitive.
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

# Citation
If you use **APEX** in a publication or presentation, please cite this repository and/or the following ISMS conference abstract:

Sven Herbers and Ha Vinh Lam Nguyen, *A New Spectral Assignment Tool and Rotational Characterization of 3,4-Dihydro-2H-Pyran-2-Carboxaldehyde (C₆H₈O₂).* 77th International Symposium on Molecular Spectroscopy (ISMS), 19 June 2024, University of Illinois at Urbana–Champaign, Urbana, IL, USA.


