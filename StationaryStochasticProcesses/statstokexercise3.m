load unknowndata
newdata = data-[3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3]';
figure(1)
plot(newdata)
%%
nd = newdata;

X=fft(nd);

n=length(nd);
Rhat=(X.*conj(X))/n;
f=[0:n-1]/n;
figure (2)
plot(f,Rhat)
%%
N = 1024;
X=fft(nd,N);
Rhat=(X.*conj(X))/n;
f=[0:N-1]/N;
figure (3)
plot(f,Rhat)
%%
rhat=ifft(Rhat);
figure (4)
plot(0:15,rhat(1:16))
%%
A =[1 -2.39 3.35 -2.34 0.96];
C =[1 0 1];
armagui
%%
e = randn(500,1); 
x = filter(modell.C, modell.A, e);

figure (5)
periodogram(x,[],4096);

figure (6)
periodogram(x,hanning(500),4096);
%%
a = length(x);
[Rhat,f]=periodogram(x,[],a,1);

figure (7)
plot(f,Rhat) % Linear scale
figure (30)
semilogy(f,Rhat) % Logarithmic scale

L = (2*a)/11;
figure (8)
pwelch(x,hanning(L),[],4096);
%%
Rhate=periodogram(e,[],4096);
Rhatew=pwelch(e,hanning(L),[],4096);


kvot = var(Rhate)/var(Rhatew);
%%
load threeprocessdata
L = (2*length(y1))/11;

figure(9)
periodogram(y1,[],4096);
figure (10)
periodogram(y1,hanning(500),4096);
figure(11)
pwelch(y1,hanning(L),[],4096);
% y1 = B
%%
L = (2*length(y2))/11;

figure(12)
periodogram(y2,[],4096);
figure (13)
periodogram(y2,hanning(500),4096);
figure(14)
pwelch(y2,hanning(L),[],4096);
% y2=C
%%
L = (2*length(y3))/11;

figure(15)
periodogram(y3,[],4096);
figure (16)
periodogram(y3,hanning(500),4096);
figure(17)
pwelch(y3,hanning(L),[],4096);
% y3=A
%%
figure (18)
mscohere(x1,y1,hanning(L),[],4096);
figure (19)
mscohere(x3,y1,hanning(L),[],4096);
%%
x1renamed = x3;

Rxy=cpsd(x1renamed,y1,hanning(L),[],4096);
Rxx=pwelch(x1renamed,hanning(L),[],4096);

figure (20)
Hf = abs(Rxy./Rxx);
plot(Hf)
figure(21)
H=tfestimate(x1renamed,y1,hanning(L),[],4096);
plot(0:2048,abs(H))
%%
H2=tfestimate(x2,y2,hanning(500),[],4096,'twosided');
h2=ifft(H2);
figure(22)
plot(0:4095,abs(h2))

H3=tfestimate(x1,y3,hanning(500),[],4096,'twosided');
h3=ifft(H3);
figure(23)
plot(0:4095,abs(h3))

h=ifft(H);
figure(24)
plot(0:2048,abs(h))
%%
y2new = y2(64:end);
x2new = x2(1:437);

H2new=tfestimate(x2new,y2new,hanning(437),[],4096,'twosided');
h2new=ifft(H2new);
figure(22)
plot(0:4095,abs(h2new))