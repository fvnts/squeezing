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
pulse2      = fopen('beta2.dat','w');

trajectory1 = fopen('trajectory1.dat','w');
trajectory2 = fopen('trajectory2.dat','w');

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

for t = -pi/2:d:pi/2-d

    n  = n + 1;

    U  = U0;
    x1 = U*x0;
    
    shadow(n)      = sqrt(1/(2*k)*U(1,1)^2 + (k/2)*U(1,2)^2);
    antishadow(n)  = sqrt((k/2)*U(2,2)^2 + 1/(2*k)*U(2,1)^2);

    fprintf(sha,'%f \n',shadow(n));
    fprintf(ash,'%f \n',antishadow(n));
    
    fprintf(pulse1,'%d %f %f %f %f\n',n,U(1,:),U(2,:));
    fprintf(trajectory1,'%f %f\n',x1(1),x1(2));

    n1 =  cos(t)^2 * (114714786 - 18293134 * cos(2*t) - 9511704 * cos(4*t) - 223035 * cos(6*t) + 15990*cos(8*t) + 1225 * cos(10*t));
    d1 =  2*(10126 + 677 * cos(2*t) - 46 * cos(4*t) - 5 * cos(6*t))^2;
            
    beta1 = n1/d1;
    
    A  = [0,1;-beta1,0];
    U0 = expm(d*A)*U;

end

%   --------------------------------
%   beta pulse 2
%   --------------------------------

U0  = eye(2);
A   = zeros(2);

n  = 0;
        
%   --------------------------------

for t = -pi/2+pi:d:pi/2+pi-d
            
    n  = n + 1;
    U  = U0;
  	x2 = U*x1;
    
    shadow(n+l)      = sqrt(1/(2*k)*U(1,1)^2 + (k/2)*U(1,2)^2);
    antishadow(n+l)  = sqrt((k/2)*U(2,2)^2 + 1/(2*k)*U(2,1)^2);
    
    fprintf(sha,'%f \n',shadow(n+l));
    fprintf(ash,'%f \n',antishadow(n+l));

    fprintf(pulse2,'%d %f %f %f %f\n',n,U(1,:),U(2,:));
    fprintf(trajectory2,'%f %f\n',x2(1),x2(2));
    
    n2 = (-38263 * cos(2*t) + 6*(4771 + 2743*cos(4*t) + 33*cos(6*t) - 26*cos(8*t)) + 49*cos(10*t))*sin(2*t)^2;
    d2 = 2*(195 * sin(t) - 2 * sin(5*t) + sin(7*t))^2;
            
    beta2 = n2/d2;
    
    A   =  [0,1;-beta2,0];
    U0  =  expm(d*A)*U;
        
end

%   --------------------------------
         
for n = 1:2*l            
    uncert  = shadow(n)*antishadow(n);
    fprintf(heisenberg,'%f \n',uncert);
end

%   --------------------------------
        
fclose(pulse1);
fclose(pulse2);

fclose(trajectory1);
fclose(trajectory2);

fclose(sha);
fclose(ash);

fclose(heisenberg);

disp('done !')

%   eof
