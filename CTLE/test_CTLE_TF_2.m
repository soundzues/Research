%descretprion: Place zero and pole atleast a decade apart for real values

%create transfer function for CTLE:
% H = (s+z)/((s+p1)(s+p2))
%put zero on location of poles 

zero_ctle = [1 5.732*10^7];

%zero_ctle = [1 10]

%chosing 2 poles at random locations > zero such that poles are at least 1
%decade apart from the zero
% p1 = 16*10^7
% p2 = 16.5*10^7
%(s+16*10^7)+(s+16.5*10^7) = (s^2 + 1.45*10^8s + 5.25*10^15)

%pole_ctle = [1 1.45*10^8 5.25*10^15];

%pole_ctle = 1;

%pole_ctle = [1 250 15000]

pole_ctle = [1 3.25*10^8 2.64*10^16]


sys_ctle = tf(zero_ctle,pole_ctle)

%{
  CTLE Transfer function:
        s + 5.732e07
  -------------------------
  s^2 + 1.45e08 s + 5.25e15
%}

%Plot Poles and Zeros plot
figure(1)
pzmap(sys_ctle)
grid on

%Plot Bode plot
figure(2)
%bode(sys_ctle)
margin(sys_ctle)
