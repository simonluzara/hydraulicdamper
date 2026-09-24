<img width="580" height="460" alt="TwinTube" src="https://github.com/user-attachments/assets/4da19c7a-c91f-49e0-bd7b-a86cd558fbb4" />

A Twin-tube Hydraulic Damper with organic oil as the compressible fluid.

The fluid flow is modeled in 1D compressible flow, including the isothermal bulk modulus of the fluid.

The compressed fluid flow from one chamber to the other through orifice and valve port flows.
The spring valves are preloaded at an initial displacement to delay the flow and increase fluid compression time.


$$
\frac{\delta \rho}{\delta P}=\frac{\rho}{\beta_f}\rvert_{T(K)}
$$
$$
    \frac{dP_1}{dt} =\frac{\beta_1(P_1)}{V_{1_0}+A_{1}x}(\frac{\bar{\rho}}{\rho_1}(Q_{ori}+Q_{valve}+Q_{leak})-A_{1}\dot{x})
$$
$$
    \frac{dP_2}{dt} =\frac{\beta_2(P_2)}{V_{2_0}-A_{2}x}(-\frac{\bar{\rho}}{\rho_2}(Q_{ori}+Q_{valve}+Q_{base valve}+Q_{leak})+A_{2}\dot{x})
$$
$$  
    \frac{dP_{3}}{dt}=\frac{\bar{\rho}}{\rho_3}\frac{Q_{ori}+Q_{basevalve}+Q_{leak}}{((V_{3_0}+V_{g_0}-(\frac{P_{g_0}}{P_{3}})^{1/1.4}V_{g_0})/\beta_f(P_{3}))+\frac{P_{g_0}^{1/1.4}V_{g_0}}{1.4P_{3}^{1+1/1.4}}}
$$

The base valve is used to avoid cavitation to replace the extruded volume of the piston-rod from the working tube. It's opening displacement is modelled as:

<img width="358" height="159" alt="image" src="https://github.com/user-attachments/assets/4ffcdb55-82fa-4874-8b73-bf8cd8f2c9ac" /> <img width="308" height="176" alt="image" src="https://github.com/user-attachments/assets/7f760cfa-d49e-4e01-95f8-226941ca2241" />



This numerical model is developped to solve the four stiff ODE's and compute the instantaneous pressure in a period of oscillation.

For more details, I recommend you to read the online version of the PhD thesis : https://univ-rennes2.hal.science/ENTPE/tel-05413895v1
