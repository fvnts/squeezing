%
%   use this code to compute
%   {beta1(t) + beta2(t) : beta0 = 0}
%   smooth evolution packets
%
%   This is supplementary material for
%   https://arxiv.org/abs/1704.05840
%   coded by J Fuentes
%   jfuentes [at] fis.cinvestav.mx
%   --------------

clc; clear; close;

%   --------------

%   file

f   =   fopen('beta1.dat','w');
g   =   fopen('beta2.dat','w');

fx1 =   fopen('x1.dat','w');
fx2 =   fopen('x2.dat','w');

fap1 =  fopen('ap1.dat','w');
fap2 =  fopen('ap2.dat','w');
fan1 =  fopen('an1.dat','w');
fan2 =  fopen('an2.dat','w');

fbp1 =  fopen('bp1.dat','w');
fbp2 =  fopen('bp2.dat','w');
fbn1 =  fopen('bn1.dat','w');
fbn2 =  fopen('bn2.dat','w');

fcp1 =  fopen('cp1.dat','w');
fcp2 =  fopen('cp2.dat','w');
fcn1 =  fopen('cn1.dat','w');
fcn2 =  fopen('cn2.dat','w');

fdp1 =  fopen('dp1.dat','w');
fdp2 =  fopen('dp2.dat','w');
fdn1 =  fopen('dn1.dat','w');
fdn2 =  fopen('dn2.dat','w');

fep1 =  fopen('ep1.dat','w');
fep2 =  fopen('ep2.dat','w');
fen1 =  fopen('en1.dat','w');
fen2 =  fopen('en2.dat','w');

%   interval

n   =   200;
d   =   pi/n;

x   =   [1;0];
x1  =   [0;0];
x2  =   [0,0];

%   --------------

ap   =   [1;2/6];
ap1  =   [0;0];
ap2  =   [0,0];

an   =   [1;-2/6];
an1  =   [0;0];
an2  =   [0,0];

%   --------------

bp   =   [1;3/6];
bp1  =   [0;0];
bp2  =   [0,0];

bn   =   [1;-3/6];
bn1  =   [0;0];
bn2  =   [0,0];

%   --------------

cp   =   [1;4/6];
cp1  =   [0;0];
cp2  =   [0,0];

cn   =   [1;-4/6];
cn1  =   [0;0];
cn2  =   [0,0];

%   --------------

dp   =   [1;5/6];
dp1  =   [0;0];
dp2  =   [0,0];

dn   =   [1;-5/6];
dn1  =   [0;0];
dn2  =   [0,0];

%   --------------

ep   =   [1;6/6];
ep1  =   [0;0];
ep2  =   [0,0];

en   =   [1;-6/6];
en1  =   [0;0];
en2  =   [0,0];


        %   ----------------------
        %   beta pulse 1
        %   ----------------------

        U0 = eye(2);
        A  = zeros(2);
        n  = 0;

        for t = -pi/2:d:pi/2-d
            
            n    = n + 1;

            U    = U0;

            x1   = U*x;

            ap1  = U*ap;
            an1  = U*an;

            bp1  = U*bp;
            bn1  = U*bn;

            cp1  = U*cp;
            cn1  = U*cn;

            dp1  = U*dp;
            dn1  = U*dn;

            ep1  = U*ep;
            en1  = U*en;

                   
            fprintf(f,'%d %f %f %f %f\n',n,U(1,:),U(2,:));

            fprintf(fx1,'%f %f\n',x1(1),x1(2));

            fprintf(fap1,'%f %f\n',ap1(1),ap1(2));
            fprintf(fan1,'%f %f\n',an1(1),an1(2));

            fprintf(fbp1,'%f %f\n',bp1(1),bp1(2));
            fprintf(fbn1,'%f %f\n',bn1(1),bn1(2));

            fprintf(fcp1,'%f %f\n',cp1(1),cp1(2));
            fprintf(fcn1,'%f %f\n',cn1(1),cn1(2));

            fprintf(fdp1,'%f %f\n',dp1(1),dp1(2));
            fprintf(fdn1,'%f %f\n',dn1(1),dn1(2));

            fprintf(fep1,'%f %f\n',ep1(1),ep1(2));
            fprintf(fen1,'%f %f\n',en1(1),en1(2));

            n1 =  cos(t)^2 * (114714786 - 18293134 * cos(2*t) - 9511704 * cos(4*t) - 223035 * cos(6*t) + 15990*cos(8*t) + 1225 * cos(10*t));
            d1 =  2*(10126 + 677 * cos(2*t) - 46 * cos(4*t) - 5 * cos(6*t))^2;

            beta1 = n1./d1;

            A   =   [0,1;-beta1,0];
            U0  =   expm(d*A)*U;

        end

        %   ----------------------
        %   beta pulse 2
        %   ----------------------

        U0 = eye(2);
        A  = zeros(2); 
        n  = 0;
        
        for t = -pi/2+pi:d:pi/2+pi-d
            
            n    = n + 1;
            
            U    = U0;

            x2   = U*x1;

            ap2  = U*ap1;
            an2  = U*an1;

            bp2  = U*bp1;
            bn2  = U*bn1;

            cp2  = U*cp1;
            cn2  = U*cn1;

            dp2  = U*dp1;
            dn2  = U*dn1;

            ep2  = U*ep1;
            en2  = U*en1;

                   
            fprintf(g,'%d %f %f %f %f\n',n,U(1,:),U(2,:));

            fprintf(fx2,'%f %f\n',x2(1),x2(2));

            fprintf(fap2,'%f %f\n',ap2(1),ap2(2));
            fprintf(fan2,'%f %f\n',an2(1),an2(2));

            fprintf(fbp2,'%f %f\n',bp2(1),bp2(2));
            fprintf(fbn2,'%f %f\n',bn2(1),bn2(2));

            fprintf(fcp2,'%f %f\n',cp2(1),cp2(2));
            fprintf(fcn2,'%f %f\n',cn2(1),cn2(2));

            fprintf(fdp2,'%f %f\n',dp2(1),dp2(2));
            fprintf(fdn2,'%f %f\n',dn2(1),dn2(2));

            fprintf(fep2,'%f %f\n',ep2(1),ep2(2));
            fprintf(fen2,'%f %f\n',en2(1),en2(2));

            n2 = (-38263 * cos(2*t) + 6*(4771 + 2743*cos(4*t) + 33*cos(6*t) - 26*cos(8*t)) + 49*cos(10*t))*sin(2*t)^2;
            d2 = 2*(195 * sin(t) - 2 * sin(5*t) + sin(7*t))^2;
            
            beta2 = n2./d2;

            A   =  [0,1;-beta2,0];
            U0  =  expm(d*A)*U;
        end
        
fclose(f);
fclose(g);

fclose(fx1);
fclose(fx2);

fclose(fap1);
fclose(fap2);
fclose(fan1);
fclose(fan2);

fclose(fbp1);
fclose(fbp2);
fclose(fbn1);
fclose(fbn2);

fclose(fcp1);
fclose(fcp2);
fclose(fcn1);
fclose(fcn2);

fclose(fdp1);
fclose(fdp2);
fclose(fdn1);
fclose(fdn2);

fclose(fep1);
fclose(fep2);
fclose(fen1);
fclose(fen2);

disp('done !')


%   eof
