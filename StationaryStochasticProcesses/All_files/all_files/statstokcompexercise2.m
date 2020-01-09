C = 1;
a1 = 0.5;
A = [1 a1];

[H,w]=freqz(C,A);
R=abs(H).^2;
figure (1)
plot(w/2/pi,R)

H=freqz(C,A,512,'whole');
Rd=abs(H).^2;
r=ifft(Rd);
figure (2)
stem([0:49],r(1:50));

%%

m = 0;
n = 400;
sigma = 1;
e = normrnd(m, sigma, 1, n);
x = filter(C, A, e);
figure (3)
plot(1:400,x)

%%

A=[1 -1 0.5]; 
C=[1 1 0.5];

P = roots(A);
Z = roots(C);

figure(5)
zplane(Z,P)

A=poly(P);
C=poly(Z);

[H,w]=freqz(C,A);
R=abs(H).^2;
figure (6)
plot(w/2/pi,R)

%%

armagui

%%
figure (10)
plot((1:160),x_temp)
%%
A = arcov(x,20);
    
    % Generate AR process by filtering white noise
    a = [1, .1, -0.8];                      % AR coefficients
    v = 0.4;                                % noise variance
    w = sqrt(v)*randn(15000,1);             % white noise
    x = filter(1,a,w);                      % realization of AR process
    [ar,vr] = arcov(x,numel(a)-1)           % estimate AR model parameters 