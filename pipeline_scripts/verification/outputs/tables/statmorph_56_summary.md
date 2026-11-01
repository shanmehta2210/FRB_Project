# statmorph vs GALFIT, 56-host working sample

Inclinations from cos^2 i = (q^2 - q0^2)/(1 - q0^2), q0 = 0.2, q = 1 - ellipticity.

- Sersic-fit ellipticity, all hosts: N=56, mean=+0.61, median=+0.14, RMS=4.02, std=4.01, 1.4826*MAD=0.83, |d|>10 deg: 3
- Sersic-fit ellipticity, flag=0 & flag_sersic=0 (45): N=45, mean=+0.72, median=+0.15, RMS=3.57, std=3.54, 1.4826*MAD=0.76, |d|>10 deg: 1
- Moment (asymmetry-centre) ellipticity, all: N=56, mean=-11.53, median=-12.13, RMS=17.21, std=12.89, 1.4826*MAD=12.81, |d|>10 deg: 34
- Moment ellipticity, flag=0 (48): N=48, mean=-10.80, median=-10.55, RMS=14.97, std=10.48, 1.4826*MAD=11.52, |d|>10 deg: 27

## Classification counts

- Lotz+2008 G-M20: {'late': 52, 'early': 3, 'merger': 1}
- F(G,M20) bulge statistic < 0 (disk side): 53 / 56
- statmorph Sersic n < 2.5: 49; GALFIT n < 2.5: 45
- Concentration C < 3.5 (Conselice 2003 late-type side): 54
