function [] = myplotter(y,tau)
if nargin < 2
    tau = 50;
end
    figure
    subplot(311)
    acf(y,tau,0.05,1);
    title('ACF')
    subplot(312)
    pacf(y,tau,0.05,1);
    title('PACF')
    subplot(313)
    normplot(y);
    title('Normplot')

end