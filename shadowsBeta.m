%   ------------------------------------------------------------/
%   uncertainty shadows
%
%   This is supplementary material for
%   https://arxiv.org/abs/1704.05840
%   coded by J Fuentes
%   jfuentes [at] fis.cinvestav.mx
%   ------------------------------------------------------------/

clc; clear all; close all;

%   files

pulse1      = fopen('beta1.dat','w');

trajectory1 = fopen('trajectory1.dat','w');

sha = fopen('shadow.dat','w');
ash = fopen('antiShadow.dat','w');

heisenberg  = fopen('heisenberg.dat','w');

%   interval

s   = 200;
d   = pi/s;
m   = 1;
k   = 1;

%   ini

q = 1;
p = 1;

x0 = [q;p]; 
x1 = [0;0];
x2 = [0,0];

l  = length(-pi/2:d:pi/2-d);

shadow      = zeros(1,2*l);
antishadow  = zeros(1,2*l);

%   --------------------------------
%   beta pulse 1
%   --------------------------------

U0 = eye(2);
A  = zeros(2);
n  = 0;

%   --------------------------------

for t = -pi/2:d:35*pi/32-d

    n  = n + 1;

    U  = U0;
    x1 = U*x0;
    
    shadow(n)      = sqrt(1/(2*k)*U(1,1)^2 + (k/2)*U(1,2)^2);
    antishadow(n)  = sqrt((k/2)*U(2,2)^2 + 1/(2*k)*U(2,1)^2);

    fprintf(sha,'%f \n',shadow(n));
    fprintf(ash,'%f \n',antishadow(n));
    
    fprintf(pulse1,'%d %f %f %f %f\n',n,U(1,:),U(2,:));
    fprintf(trajectory1,'%f %f\n',x1(1),x1(2));

    n1 = cos(t)^2*(31857611 - 32827504*cos(2*t) + 12927000*cos(4*t) + 6879216*cos(6*t) - 4452*cos(8*t) - 14768*cos(10*t) + 45736*cos(12*t) + 10800*cos(14*t) + 729*cos(16*t)); 
    d1 = 8*(1553 - 20*cos(2*t) - 20*cos(4*t) + 20*cos(6*t) + 3*cos(8*t))^2;           
    
    beta1 = n1/d1;
    
    A  = [0,1;-beta1,0];
    U0 = expm(d*A)*U;

end

%   --------------------------------
         
for n = 1:2*l            
    uncert  = shadow(n)*antishadow(n);
    fprintf(heisenberg,'%f \n',uncert);
end

%   --------------------------------
        
fclose(pulse1);

fclose(trajectory1);

fclose(sha);
fclose(ash);

fclose(heisenberg);

disp('done !')

%   eof
