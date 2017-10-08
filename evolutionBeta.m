%
%   use this code to compute
%   {beta1(t) + beta2(t) : beta0 = 0}
%   smooth evolution packets
%
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

fx =   fopen('x1.dat','w');

fap =  fopen('ap1.dat','w');
fan =  fopen('an1.dat','w');

fbp =  fopen('bp1.dat','w');
fbn =  fopen('bn1.dat','w');

fcp =  fopen('cp1.dat','w');
fcn =  fopen('cn1.dat','w');

fdp =  fopen('dp1.dat','w');
fdn =  fopen('dn1.dat','w');

fep =  fopen('ep1.dat','w');
fen =  fopen('en1.dat','w');

%   interval

n   =   200;
d   =   pi/n;

x   =   [1;0];
x1  =   [0;0];

%   --------------

ap   =   [1;2/6];
ap1  =   [0;0];

an   =   [1;-2/6];
an1  =   [0;0];


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
        %   beta pulse
        %   ----------------------

        U0 = eye(2);
        A  = zeros(2);
        n  = 0;

        for t = -pi/2:d:35*pi/32-d 
            
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

            fprintf(fx,'%f %f\n',x1(1),x1(2));

            fprintf(fap,'%f %f\n',ap1(1),ap1(2));
            fprintf(fan,'%f %f\n',an1(1),an1(2));

            fprintf(fbp,'%f %f\n',bp1(1),bp1(2));
            fprintf(fbn,'%f %f\n',bn1(1),bn1(2));

            fprintf(fcp,'%f %f\n',cp1(1),cp1(2));
            fprintf(fcn,'%f %f\n',cn1(1),cn1(2));

            fprintf(fdp,'%f %f\n',dp1(1),dp1(2));
            fprintf(fdn,'%f %f\n',dn1(1),dn1(2));

            fprintf(fep,'%f %f\n',ep1(1),ep1(2));
            fprintf(fen,'%f %f\n',en1(1),en1(2));

            n1 = cos(t)^2*(31857611 - 32827504*cos(2*t) + 12927000*cos(4*t) + 6879216*cos(6*t) - 4452*cos(8*t) - 14768*cos(10*t) + 45736*cos(12*t) + 10800*cos(14*t) + 729*cos(16*t)); 
            d1 = 8*(1553 - 20*cos(2*t) - 20*cos(4*t) + 20*cos(6*t) + 3*cos(8*t))^2;
            
            beta1 = n1./d1;

            A   =   [0,1;-beta1,0];
            U0  =   expm(d*A)*U;

        end


        
fclose(f);

fclose(fx);

fclose(fap);
fclose(fan);

fclose(fbp);
fclose(fbn);

fclose(fcp);
fclose(fcn);

fclose(fdp);
fclose(fdn);

fclose(fep);
fclose(fen);

disp('done !')


%   eof
