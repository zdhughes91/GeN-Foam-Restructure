%%Unkwons of the problem to find: Tin,Tout,Tmean and the power after the transient%%

Power_in=990000000;	%W
U_in=0.831669;	%m/s
Density=4125;	%kg/m^3
cp=1600;	%J/kg*m^3
Area=1;	%m^2
massflow_in=Area*Density*U_in;	%kg/s %initial massflow rate
		
Tin_in=963.97;	%K
Tout_in=1144.3;	%K 
T_hx=900;	%K   %Temperature of the hx
Tml_in=(Tout_in-Tin_in)/log((Tout_in-T_hx)/(Tin_in-T_hx)); %K  initial logarithmic mean temperature
Tmean_in=(Tin_in+Tout_in)/2;  %K
Condut=200*2*1.85*10^(4);	%W/K  %conductance of the hx  
		
U_red=0.297515;   %m/s	

%Compute the new massflow rate

massflow_fin=U_red*Area*Density; %kg/s
		
alpha_fuel=-3.54253*10^(-5);	%1/K  %feedback fuel coefficient
		
prec_ini=174.01*10^(-5); %-  %reactivity of the precursors on the initial steadystate
prec_fin=134.221*10^(-5); %- %reactivity of the precursors after the transient 

%Compute Tmean_fin

Tmean_fin=((prec_fin-prec_ini)/(alpha_fuel))+Tmean_in %K

%Compute Tin_fin 

A=0.5*(exp(Condut/(massflow_fin*cp))-1);  %coefficient found thanks to the power balance equations 
                                           %on the system and on the hx
                                         
Tin_fin=(Tmean_fin+A*T_hx)/(A+1)  %K  

%Compute Tout_fin

Tout_fin=2*Tmean_fin-Tin_fin  %K

Tml_fin=(Tout_fin-Tin_fin)/log((Tout_fin-T_hx)/(Tin_fin-T_hx)); %K  final logarithmic mean temperature

%Compute final power

Power_fin=massflow_fin*cp*(Tout_fin-Tin_fin)     %W