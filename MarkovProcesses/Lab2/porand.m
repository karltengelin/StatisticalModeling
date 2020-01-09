function x=porand(lambda);
% n=porand(m);
% Ger ett Poissonfördelat slumptal ur Po(m), där
% m  är väntevärdet.

% Finn Lindgren 1998-11-11
% FL            1999-11-24

% Fixa till numeriska problem för stora lambda
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
