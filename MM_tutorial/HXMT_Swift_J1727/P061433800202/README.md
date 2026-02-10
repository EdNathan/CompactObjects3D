
Here you'll find six spectra (with corresponding matrices) and two .xcm files for Xspec to fit simultaneously:

Power spectrum band 1 (PS1)
Power spectrum band 2 (PS2)
Real part of the cross spectrum of band 1 vs. band 2 (RE)
Imaginary part of the cross spectrum of band 1 vs. band 2 (IM)

The file [P061433800202_pei_9lor.xcm](https://github.com/EdNathan/CompactObjects3D/tree/main/MM_tutorial/HXMT_Swift_J1727/P061433800202/P061433800202_pei_9lor.xcm) uses a linear combination of Lorentzians to fit PS1 and PS2 and a combination of Lorentzians multiplied by the cosine (sine) of a parameter (the phase lag) to fit the RE (IM) part of the cross spectrum.

The .xcm file also reads the phase-lag (PL) spectrum and the coherence function (CF)
and computes the derived model of these two on the basis of the model fitted to
PS1, PS2, RE and IM.

See for [Méndez et al. 2024](https://ui.adsabs.harvard.edu/abs/2024MNRAS.527.9405M/abstract) for details


