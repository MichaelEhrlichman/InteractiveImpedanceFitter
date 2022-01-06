# InteractiveImpedanceFitter
Fit data to detector-to-detector phase advance using scipy optimizers and an interactive matplotlib front end.  Python calls Bmad library through ctypes.

The blue bars are the strength of the impedance at the expected strongest impedance sources in the machine.  The impedance is in units k1/mA, or in other words, focusing per beam current.  Higher beam current yields a stronger focusing effect from that element.

The best result was found by adjusting the impedances manually using the front end, then using a local optimizer to refine the result.

![alt text](GUI.png)
