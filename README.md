# pypngw2

Post-Newtonian GW waveform generator

Installation && compile

```
python3 setup.py build_py
python3 setup.py install
```

sometimes the lib cannot use in Apple core macbook.

Run in python

Run in python

```python
from pypngw2 import calculate_strain_SPA
m1 = 12
m2 = 10
chi1 = 0
chi2 = 0
iota = 45 * np.pi/180
phic = 0
e0 = 0.3
distance = 100
fmin = 5
deltaF = 1
detname = 'ET1'
psi = 0
ra = 0
dec = 0
freqs, htilde = calculate_strain_SPA(m1, m2, chi1, chi2, 
                iota, phic, e0, distance, fmin, deltaF,
                detname, psi, ra, dec)
```

see more details and description in \_\_init\_\_.py
