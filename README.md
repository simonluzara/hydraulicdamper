<img width="580" height="460" alt="TwinTube" src="https://github.com/user-attachments/assets/4da19c7a-c91f-49e0-bd7b-a86cd558fbb4" />

A Twin-tube Hydraulic Damper with organic oil as the compressible fluid.

The fluid flow is modeled in 1D compressible flow, including the isothermal bulk modulus of the fluid.

The compressed fluid flow from one chamber to the other through orifice and valve port flows.
The spring valves are preloaded at an initial displacement to give a certain resistance.


$$
\frac{\delta \rho}{\delta P}=\frac{\rho}{\beta_f}\rvert_{T(K)}
$$
$$
    \frac{dP_1}{dt} =\frac{\beta_1(P_1)}{V_{1_0}+A_{1}x}(\frac{\bar{\rho}}{\rho_1}(Q_{ori}+Q_{valve}+Q_{leak})-A_{1}\dot{x})
$$
$$
    \frac{dP_2}{dt} =\frac{\beta_2(P_2)}{V_{2_0}-A_{2}x}(-\frac{\bar{\rho}}{\rho_2}(Q_{ori}+Q_{valve}+Q_{base valve}+Q_{leak})+A_{2}\dot{x})
$$
