function [dpdt,Force,Q_v,v_open,Q_leakage]=damper_comp(t,P,amp,freq,pars)
% Positive direction for rebound, negative direction for compression (as in MTS).
%--Pressure in the 3 Chambers--%
P_ex=P(1);                     % Pressure in the rebound chamber
P_co=P(2);                     % Pressure in the compression chamber
P_res=P(3);                    % Pressure in the reservoir
%------------------------------%
%--Call out the constant parameters from the function--%%

[d_p,d_r,A_p,A_r,A_tp,h_p,h_gu,Delta_p,d_gu,l_sealpist,delta_p,l_sealrod,delta_r,...
d_t,...
d_ori_co,d_ori_re,d_ori_ba,l_ori_co,l_ori_re,l_ori_ba,~,Cfmax_ori_co,Cfmax_ori_re,Cfmax_ori_ba,...
r_v_co,l_v_co,r_ex_co,A_ex_co,k_co,x0_co,stroke_co,l_s_co,delta_co,...
r_v_re,l_v_re,r_ex_re,A_ex_re,k_re,x0_re,stroke_re,l_s_re,delta_re,...
r_v_ba,l_v_ba,r_ex_ba,A_ex_ba,k_ba,x0_ba,stroke_ba,l_s_ba,delta_ba,...
d_o,l_open,Dh_o,Cfmax_o,d_film_1,d_film_2,A_bot,A_top,C_lang,k_in,h_0,...
L_co,Vco_0,L_re,Vre_0,...
Vg_0,Pg_0,Vres_0,...
V_a_0,P_cri,P_vap,x_0,alpha_T,rho_f_0,rho_a_0,T,nu_Pa_T]=pars{1:78};

N_ori=length(d_ori_co)+length(d_ori_ba);
%-----------------------------------------%

%--Imposed Sinusodial Excitation--%
x=round(amp*sin(2*pi*freq*t+pi/2),8);            % Displacement from the middle of the piston 

u=round(amp*2*pi*freq*cos(2*pi*freq*t+pi/2),8);    % Velocity of the piston rod
%---------------------------------%
%--Atmospheric Surrounding Temp and Pressure--%
T_atm=15;                                        % Fluid Initial Properties by the supplier
p_atm=101325;
%---------------------------------------------%Acsmccsm

dpdt=zeros(3,1);
Q_v=zeros(1,6);
v_open=zeros(1,6);

[rho1,beta1,mu1]=props(P_ex,T,nu_Pa_T);
[rho2,beta2,mu2]=props(P_co,T,nu_Pa_T);
[rho3,beta3,mu3]=props(P_res,T,nu_Pa_T);
   
[rho_bar,~,mu_bar]=props(P_ex/2+P_co/2,T,nu_Pa_T);

[rho_bar_bas,~,mu_bar_bas]=props(P_res/2+P_co/2,T,nu_Pa_T);

if (freq<0.3)              % Return flow through the orifices and leak
    zeta=1;
elseif (freq>0.3 && freq<0.7)
    zeta=0.6;
else
    zeta=0.8;
end

[Q_v(1),v_open(1)]=valve(P_co,P_ex,rho_bar,mu_bar,r_v_co(1),l_v_co,r_ex_co(1),A_ex_co(1),k_co(1),x0_co(1),stroke_co(1),l_s_co,delta_co,3);
[Q_v(2),v_open(2)]=valve(P_co,P_ex,rho_bar,mu_bar,r_v_co(2),l_v_co,r_ex_co(2),A_ex_co(2),k_co(2),x0_co(2),stroke_co(2),l_s_co,delta_co,3,d_ori_co(1),l_ori_co(1),Cfmax_ori_co(1));
[Q_v(3),v_open(3)]=valve(P_co,P_ex,rho_bar,mu_bar,r_v_co(3),l_v_co,r_ex_co(3),A_ex_co(3),k_co(3),x0_co(3),stroke_co(3),l_s_co,delta_co,3,d_ori_co(2),l_ori_co(2),Cfmax_ori_co(2));

[Q_v(4),v_open(4)]=valve_emb(P_co,P_res,rho_bar_bas,mu_bar_bas,r_v_ba(1),l_v_ba,r_ex_ba(1),A_ex_ba(1),k_ba(1),x0_ba(1),stroke_ba(1),l_s_ba,delta_ba,3,d_ori_ba,l_ori_ba,Cfmax_ori_ba);
[Q_v(5),v_open(5)]=valve_emb(P_co,P_res,rho_bar_bas,mu_bar_bas,r_v_ba(1),l_v_ba,r_ex_ba(2),A_ex_ba(2),k_ba(2),x0_ba(2),stroke_ba(2),l_s_ba,delta_ba,3);
[Q_v(6),v_open(6)]=valve_emb(P_co,P_res,rho_bar_bas,mu_bar_bas,r_v_ba(1),l_v_ba,r_ex_ba(3),A_ex_ba(3),k_ba(3),x0_ba(3),stroke_ba(3),l_s_ba,delta_ba,3);

Q_leakage=zeta*leak(P_co,P_ex,mu_bar,d_p,d_t,h_p,u)...
    +reflux(P_co,P_ex,rho_bar,mu_bar,r_v_re(2),l_v_re,r_ex_re(2),A_ex_re(2),k_re(2),x0_re(2),stroke_re(2),l_s_re,delta_re,3,d_ori_re,l_ori_re,Cfmax_ori_re);      %Leakage of the piston seals

Q_21=Q_v(1)+Q_v(2)+Q_v(3)+...
    Q_leakage;      %Leakage of the piston seals
Q_23=Q_v(4)+Q_v(5)+Q_v(6);

dpdt(1)=beta1/(Vre_0+A_tp*x)*(rho_bar/rho1*Q_21...              
-A_tp*u);

dpdt(2)=beta2/(Vco_0-A_p*x)*(-rho_bar/rho2*Q_21-rho_bar_bas/rho2*Q_23...              
+A_p*u);

dpdt(3)=(rho_bar_bas/rho3*Q_23)/...
(((Vres_0+Vg_0-((Pg_0/P_res)^(1/1.4))*Vg_0)/beta3)+(Pg_0^(1/1.4)*Vg_0/1.4/P_res^(1+1/1.4)));

Force_dyn=P_ex*A_tp - P_co*A_p;

if (freq<0.6 || freq>3)
    mu1=0.1;
    mu2=0.1;
else
    mu1=0.3;
    mu2=0.3;
end

F_fr=mu1*pi*d_p*l_sealpist*(P_ex-P_co)*0.5 - mu2*pi*d_r*l_sealrod*abs(p_atm-P_ex)*0.5;

if u==0
    F_dry=0;
else
    F_dry=5 + (60-5).*exp(-0.25*(u/0.01).^2);      %Dry friction from low speed experimental observation.
end
Force=Force_dyn+F_fr-F_dry;


    function [rho,Beta_f,mu_f]=props(Pre,Temp,nu_Pa_T)
        x_0=0.001;                                 % Air content percentage in the oil at STP
        alpha_T=7.2e-4;                                % Volumetric thermal expansion coeff.  (/K)
        rho_f_0=872.0;                                   % oil density at stp
        p_atm=1.01325e5;                             % 1 atm pressure of airbag at stp
        T_atm=15.0;                                  % Temperature at stp.°C
        
    %     bulk_liq=1.517e9;
        bulk_liq=(10^((0.3766*(log10(nu_Pa_T)^0.3307))-0.2766)...
                +(5.851-0.01382*Temp)*(Pre/1e9))*1e9;                              %Bulk at operating temperature and pressure in Pa [L. Duda]
        mu_0=(nu_Pa_T*1e-6)*rho_f_0*exp(-alpha_T*(Temp-T_atm));              %Dynamic viscosity of liquid oil
        
    %     Beta_f=(1+x_0*(p_atm./Pre))./(1+x_0*(p_atm./Pre.^2)*bulk_liq)*bulk_liq;
        eta=(x_0/(1-x_0))*(p_atm/Pre)^(1/1.4)*(Temp+273.15)/(T_atm+273.15);                           % Entrained air percentage
        Beta_f=(1.4*Pre*bulk_liq)/((1-eta)*1.4*Pre+eta*bulk_liq);
        rho=rho_f_0*(1./(1+eta))*exp((Pre-p_atm)/Beta_f-alpha_T*(Temp-T_atm));
        mu_f=mu_0.*(1+1.5*eta).*exp(2.3e-8*(Pre-p_atm)-4e-2*(Temp-T_atm));
    end


    function Q_leak=leak(P1,P2,mu,d_pist,d_tu,L_pist,u_p)
        v1=-1;
        v2=1;
        k=10000;

        ka=d_pist/d_tu;
        r2=d_tu/2;

        Q_leak=(pi*r2^4/8/mu/L_pist*abs(P1-P2)*(1-ka^4-(1-ka^2)^2/log(1/ka)) + ...
            pi*r2^2/2*((1-ka^2)/log(1/ka)-2*ka^2)*abs(u_p))*...
            ((v2-v1)*tanh(k*(P1-P2))/2 + (v1+v2)/2);
        % Q_leak=(pi*r2*Delta_p^3/6/mu/L_pist*abs(P1-P2)+pi*r2*Delta_p*u_p)*...
        % ((v2-v1)*tanh(k*(P1-P2))/2 + (v1+v2)/2);
    end

    function Q_ori=ori_reverse(P1,P2,d_ori,l_ori,cfmax,rho_P1,mu_P1)
    
    v1=-1;
    v2=1;
    k=10000;
    gamma_t=500;

    gamma=d_ori./mu_P1*sqrt(2*abs(P1-P2)*rho_P1);
    

    dP=P1-P2;
    c_f = cfmax.*tanh(1.8619*gamma./gamma_t);
    Q_ori=sum(c_f.*(pi./4*d_ori.^2)*sqrt(2*abs(dP)/rho_P1)*...
        ((v2-v1)*tanh(k*(P1-P2))/2 + (v1+v2)/2));
    
    end

    function Q_reflux=reflux(P1,P2,rho,mu,r_valve,l_valve,r_ext,A_ext,k_valve,x0_pre,stroke,l_spool,delta_spool,N_holes,d_ori,l_ori,cfmax_ori)
        if d_ori<0.9e-3
            alph=0.4;
        elseif d_ori<1e-3
            alph=0.6;
        else
            alph=0.99;
        end
        r_int=r_ext-l_valve;
        dP1=0.6*(1-(pi*r_int^2/A_tp))*rho/2*A_p^2*u^2/(pi*r_int^2)^2/9;
        dP2=0.6*(1-(2.5e-3^2/(2*r_int)^2))*rho/2*A_p^2*u^2/(pi/4*2.5e-3^2)^2/9;
        dP3=A_p^2*u^2*rho/2/(3*(pi/4*4e-3^2/(1+0.7071068*sqrt(1-(pi/4*4e-3^2)/A_p)-pi/4*4e-3^2/(pi/4*11e-3^2))))^2;     % Area Expansion from the valve to the piston side.

        Q_reflux=alph*ori_reverse(P1-dP1-dP2,P2+dP3,d_ori,l_ori,cfmax_ori,rho,mu);

    end

function [Q_valve,valve_o]=valve(P1,P2,rho,mu,r_valve,l_valve,r_ext,A_ext,k_valve,x0_pre,stroke,l_spool,delta_spool,N_holes,d_ori,l_ori,cfmax_ori)
  
    k=1e6;
    v1=-1;
    gamma_v=1200;
    r_int=r_ext-l_valve;
    if nargin==17
        
        dP1=((1-(pi*r_int^2/A_p))+0.7072*sqrt(1-(pi*r_int^2/A_p)))^2*rho/2*A_p^2*u^2/(pi*r_int^2)^2/9;      % Area contraction of flow through 2 valves

        dP3=A_p^2*u^2*rho/2/(6*(pi/4*3e-3^2/(1+0.7071068*sqrt(1-(pi/4*3e-3^2)/(pi/4*11e-3^2))-pi/4*3e-3^2/A_tp))+...
        3*(pi/4*4e-3^2/(1+0.7071068*sqrt(1-(pi/4*4e-3^2)/(pi/4*11e-3^2))-pi/4*4e-3^2/A_tp)))^2;    % Area Expansion from the valve to the piston side.

        F_m=rho*(A_p*u/N_ori)^2/(pi*r_int^2);
        x_p=(F_m+P1*(pi/4*8.5e-3^2)-dP1*(pi*r_int^2)-(P2+dP3)*(pi/4*9e-3^2))/k_valve-x0_pre;   %Static movement of the spring.

%         x_p=(P1*(pi/4*9e-3^2)-dP1*(pi*r_int^2)-(P2+dP3)*(pi/4*9e-3^2))/k_valve-x0_pre;   %Static movement of the spring.
	    % x_p=(F_m +(P1-dP1)*(pi/4*8.5e-3^2)-(P2+dP3)*(pi/4*9e-3^2))/k_valve-x0_pre;   %Static movement of the spring.

    %     x_p=((P1-dP1-P2-dP3)*(pi/4*9e-3^2))/k_valve-x0_pre;   %Static movement of the spring.


    if x_p<=stroke
        valve_o=0;
		Q_leak=pi*r_ext*delta_spool^3/6/mu/l_spool*(P1-P2);    %Leakage of the valve.
        
        dP2=0.22*(1-(2.5e-3^2/(2*r_int)^2))*rho/2*(A_p*abs(u)/N_ori-Q_leak)^2/(pi/4*2.5e-3^2)^2;      % Area contraction of flow through 3 valves
        Q_ori=ori_reverse(P1-dP1-dP2,P2+dP3,d_ori,l_ori,cfmax_ori,rho,mu);

        Q_valve=Q_leak+Q_ori;

    else
		[A_s,D_h]=area_valve(x_p-stroke,r_valve);
        valve_o=min(x_p*1e3-stroke*1e3,2*r_valve*1e3);

        Cfmax_v=0.34;

        dP=abs(P1-dP1-P2-dP3);

        Cfmax_v= Cfmax_v*tanh(2*(x_p-stroke)/(2*r_valve));
	    c_f=Cfmax_v*tanh(1.86*D_h/mu*sqrt(2*dP*rho)/gamma_v);

		Q_hole=c_f*N_holes*A_s*sqrt(2*dP/rho);       %Flow through the valve holes and leakage
            
        dP2=0.22*(1-(2.5e-3^2/(2*r_int)^2))*rho/2*(A_p*abs(u)/2-Q_hole)^2/(pi/4*2.5e-3^2)^2;      % Area contraction of flow through 5 valves
        Q_ori=ori_reverse(P1-dP1-dP2,P2+dP3,d_ori,l_ori,cfmax_ori,rho,mu);

        Q_valve=Q_hole+Q_ori;
    end
    else
        if ((P1)*(pi/4*8.5e-3^2)-(P2)*(pi/4*9e-3^2))<=k_valve*(x0_pre) %Static movement of the spring.
        
            valve_o=0;
		    Q_leak=pi*r_ext*delta_spool^3/6/mu/l_spool*(P1-P2);    %Leakage of the valve.
            Q_valve=Q_leak;
        else
        dP1=((1-(pi*r_int^2/A_p))+0.7072*sqrt(1-(pi*r_int^2/A_p)))^2*rho/2*A_p^2*u^2/(pi*r_int^2)^2/9;      % Area contraction of flow through 2 valves

        dP3=A_p^2*u^2*rho/2/(6*(pi/4*3e-3^2/(1+0.7071068*sqrt(1-(pi/4*3e-3^2)/(pi/4*11e-3^2))-pi/4*3e-3^2/A_tp))+...
        3*(pi/4*4e-3^2/(1+0.7071068*sqrt(1-(pi/4*4e-3^2)/(pi/4*11e-3^2))-pi/4*4e-3^2/A_tp)))^2;    % Area Expansion from the valve to the piston side.

        F_m=rho*(A_p*u/N_ori)^2/(pi*r_int^2);
        x_p=(F_m +P1*(pi/4*8.5e-3^2)-dP1*(pi*r_int^2)-(P2+dP3)*(pi/4*9e-3^2))/k_valve-x0_pre;   %Static movement of the spring.
%         x_p=((P1-dP1)*(pi/4*8.5e-3^2)-(P2+dP3)*(pi/4*9e-3^2))/k_valve-x0_pre;   %Static movement of the spring.
        if x_p<=stroke
		    Q_leak=pi*r_ext*delta_spool^3/6/mu/l_spool*(P1-P2);    %Leakage of the valve.
            Q_valve=Q_leak;
            valve_o=0;
        else

		    [A_s,D_h]=area_valve(x_p-stroke,r_valve);
            valve_o=min(x_p*1e3-stroke*1e3,2*r_valve*1e3);
    
            Cfmax_v=0.34 + 0.2/(1+exp(-100*((2*pi*amp*freq)-0.308)));            

            dP=abs(P1-dP1-P2-dP3);
    
            Cfmax_v= Cfmax_v*tanh(2*(x_p-stroke)/(2*r_valve));
	        c_f=Cfmax_v*tanh(1.86*D_h/mu*sqrt(2*dP*rho)/gamma_v);
    
		    Q_leak=pi*r_ext*delta_spool^3/6/mu/l_spool*(P1-P2);    %Leakage of the valve.
    
            Q_hole=c_f*N_holes*A_s*sqrt(2*dP/rho);       %Flow through the valve holes and leakage
            Q_valve=Q_hole+Q_leak;
        end

        end
    end
    end


  function [Q_valve,valve_o]=valve_emb(P1,P2,rho,mu,r_valve,l_valve,r_ext,A_ext,k_valve,x0_pre,stroke,l_spool,delta_spool,N_holes,d_ori,l_ori,cfmax_ori)
  
    k=1e6;
    v1=0;
    gamma_v=1200;
    r_int=r_ext-l_valve;
    if nargin==17
        dP1=((1-(pi*r_int^2/A_p))+0.7072*sqrt(1-(pi*r_int^2/A_p)))^2*rho/2*A_p^2*u^2/(pi*r_int^2)^2/9;      % Area contraction of flow through 2 valves

    %    x_p=(P1*A_ext-dP1*(pi*r_int^2)-(P2)*(pi/4*6e-3^2))/k_valve-x0_pre;   %Static movement of the spring.
        F_m=rho*(A_p*u/N_ori)^2/(pi*r_int^2);
        x_p=(F_m +P1*(pi/4*8.5e-3^2)-dP1*(pi*r_int^2)-(P2)*(pi/4*9e-3^2))/k_valve-x0_pre;   %Static movement of the spring.
%         x_p=(F_m+(P1-dP1)*(pi/4*8.5e-3^2)-(P2)*(pi/4*9e-3^2))/k_valve-x0_pre;   %Static movement of the spring.

    if x_p<=stroke
        valve_o=0;
		Q_leak=pi*r_ext*delta_spool^3/6/mu/l_spool*(P1-P2);    %Leakage of the valve.
        
        dP2=0.22*(1-(2.5e-3^2/(2*r_int)^2))*rho/2*(A_p*abs(u)/N_ori-Q_leak)^2/(pi/4*2.5e-3^2)^2;      % Area contraction of flow through 3 valves
        Q_ori=ori_reverse(P1-dP1-dP2,P2,d_ori,l_ori,cfmax_ori,rho,mu);

        Q_valve=Q_leak+Q_ori;

    else
        
	    [A_s,D_h]=area_valve(x_p-stroke,r_valve);
        valve_o=min(x_p*1e3-stroke*1e3,2*r_valve*1e3);
        
        Cfmax_v=0.42 + 0.3/(1+exp(-100*((2*pi*amp*freq)-0.308)));

        dP=abs(P1-dP1-P2);

        Cfmax_v= Cfmax_v*tanh(2*(x_p-stroke)/(2*r_valve));
	    c_f=Cfmax_v*tanh(1.8619*D_h/mu*sqrt(2*dP*rho)/gamma_v);

        Q_hole=c_f*N_holes*A_s*sqrt(2*dP/rho);       %Flow through the valve holes and leakage
        
        dP2=0.22*(1-(2.5e-3^2/(2*r_int)^2))*rho/2*(A_p*abs(u)/2-Q_hole)^2/(pi/4*2.5e-3^2)^2;      % Area contraction of flow through 5 valves
        Q_ori=ori_reverse(P1-dP1-dP2,P2,d_ori,l_ori,cfmax_ori,rho,mu);

        Q_valve=Q_hole+Q_ori;

    end
    else
	    if ((P1)*(pi/4*8.5e-3^2)-(P2)*(pi/4*9e-3^2))<=k_valve*(x0_pre)   %Static movement of the spring.

        valve_o=0;
		Q_leak=pi*r_ext*delta_spool^3/6/mu/l_spool*(P1-P2);    %Leakage of the valve.
        Q_valve=Q_leak;

    else
        dP1=((1-(pi*r_int^2/A_p))+0.7072*sqrt(1-(pi*r_int^2/A_p)))^2*rho/2*A_p^2*u^2/(pi*r_int^2)^2/9;      % Area contraction of flow through 2 valves

        F_m=rho*(A_p*u/N_ori)^2/(pi*r_int^2);
        x_p=(F_m +P1*(pi/4*8.5e-3^2)-dP1*(pi*r_int^2)-(P2)*(pi/4*9e-3^2))/k_valve-x0_pre;   %Static movement of the spring.
%         x_p=(F_m+(P1-dP1)*(pi/4*8.5e-3^2)-(P2)*(pi/4*9e-3^2))/k_valve-x0_pre;   %Static movement of the spring.
            if x_p<=stroke
    		    Q_leak=pi*r_ext*delta_spool^3/6/mu/l_spool*(P1-P2);    %Leakage of the valve.
                Q_valve=Q_leak;
                valve_o=0;
            else

	            [A_s,D_h]=area_valve(x_p-stroke,r_valve);
                valve_o=min(x_p*1e3-stroke*1e3,2*r_valve*1e3);
                
                Cfmax_v=0.42 + 0.3/(1+exp(-100*((2*pi*amp*freq)-0.308)));
    
                dP=abs(P1-dP1-P2);
        
                Cfmax_v= Cfmax_v*tanh(2*(x_p-stroke)/(2*r_valve));
	            c_f=Cfmax_v*tanh(1.8619*D_h/mu*sqrt(2*dP*rho)/gamma_v);
        
		        Q_leak=pi*r_ext*delta_spool^3/6/mu/l_spool*(P1-P2);    %Leakage of the valve.
        
                Q_hole=c_f*N_holes*A_s*sqrt(2*dP/rho);       %Flow through the valve holes and leakage
                Q_valve=Q_hole+Q_leak;
            end
    end

    end
    end

    function [A,D_h]=area_valve(x,r)
    x=min(x,2*r);
    A=r^2*acos(1-(x/r))-(r-x)*sqrt(2*x*r-x^2);
    Peri=2*r*acos(1-(x/r))+1e-10;
    D_h=4*A/Peri;

    end
end

