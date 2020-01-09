function x=porand(lambda);
% n=porand(m);
% Ger ett Poissonf�rdelat slumptal ur Po(m), d�r
% m  �r v�ntev�rdet.

% Finn Lindgren 1998-11-11
% FL            1999-11-24

% Fixa till numeriska problem f�r stora lambda
if lambda>80
 x=porand(lambda/2)+porand(lambda/2);
else
 u=rand(1,1);
 x=0;
 S=exp(-lambda);
 fakult=1;
 while u>S
   x=x+1;
   fakult=fakult*x;
   S=S+lambda^x/fakult*exp(-lambda);
 end
end
