clear
close all
%{
    The section below is used to initialize all variables needed for the
    model. We used a time step of 1/100 of a millisecond, which may be
    overkill, but gives more accurate simulations of model neurons.
%}
gK= 36; %maximal K+ conductance
gNa=120; %maximal Na+ conductance
gL=.3; %Leak conductance
Ek=-12; % Reversal potential of Potassium
ENa=115;% Reversal potential of Sodium
El=10.6;% Reversal potential of the leak channels
Vrest=0;% resting potential
dt=.01; % time step for Eulers method
length=10000; % amount of time steps (100 ms / .01 time step)
time= .01:.01:100;
Vm= zeros(1,length); % Zero out the voltage array
Vm(1)=Vrest;%  start the voltage at rest
Cm = 1.0; %membrane capacitance
Icurrent=[00.*ones(1, length), zeros(1, length)];%current injection parameter
%{
    The section below is used to initialize the variables needed for
    euler's method to iterate.
%}
alphaM(1)=0.1.*((25-Vm(1))/(exp((25-Vm(1))/10)-1));%this is used to determine initial AlphaM to find original gating parameters
betaM(1)=4.*exp(-Vm(1)/18);

alphaN(1)=.01.*((10-Vm(1))/(exp((10-Vm(1))/10)-1));
betaN(1)= .125.*exp(-Vm(1)/80);

alphaH(1)=.07.*exp(-Vm(1)/20);
betaH(1)=1/(exp((30-Vm(1))/10)+1);

m(1)=alphaM(1)/(alphaM(1)+betaM(1));
n(1)=alphaN(1)/(alphaN(1)+betaN(1));
h(1)=alphaH(1)/(alphaH(1)+betaH(1));
%{
    Here begins the numerical solving of the differential equations that
    describe the HH model neuron. This iterates from 1 to the length of the
    simulation time points, which is 100/millisecond for 100 ms so 10000
    total iterations.
%}
for i= 1:length-1
    %% Finding alphas and betas, the opening and closing characteristics of the markov model basis of ion channels
    alphaM(i)= 0.1*((25-Vm(i))/(exp((25-Vm(i))/10)-1));%alphaM determined by membrane voltage
    betaM(i)= 4 * exp(-(Vm(i))/18);%betaM determined by membrane voltage only
    
    alphaN(i)= 0.01 * ((10-Vm(i))/(exp((10-Vm(i))/10)-1));%alphaN determined by membrane voltage only
    betaN(i)= 0.125*exp(-(Vm(i))/80);%betaN determined by membrane voltage only
    
    alphaH(i)= 0.07*exp(-(Vm(i))/20);%AlphaH determined by membrane voltage only
    betaH(i)= 1/(exp((30-Vm(i))/10)+1);%betaH determined by membrane voltage only
    %% Finding m, n, and h which determine the odds of a single Na+ activation gate, K+ activation gate
    %  and Na+ inactivation gate respectively being open. Based on the
    %  alphas and betas, and previous values for M, N, and H. Eulers method
    %  is used to numerically solve differential equations
    m(i+1) = m(i) + dt*(alphaM(i) *(1-m(i)) - betaM(i)*m(i));
    n(i+1) = n(i) + dt*(alphaN(i) *(1-n(i)) - betaN(i)*n(i));
    h(i+1) = h(i) + dt*(alphaH(i) *(1-h(i)) - betaH(i)*h(i));
    %% Sodium, potassium, and leak currents are determined using channel
    % characteristics discussed above. Models of Na+ channel include 3
    % activation gates and one inactivation gate, so cubing the activation
    % gate shows the odds that a single channel is open, which due to the
    % law of large numbers can describe total conductance as a percent of
    % max.
    gNaActual(i)=(m(i+1)^3)*gNa*h(i+1);
    gKActual(i)=((n(i+1))^4).*gK;
    INa(i)= (m(i+1)^3)*gNa*h(i+1)*(Vm(i)-ENa);
    IK(i)= ((n(i+1))^4) *gK*(Vm(i)-Ek);
    IL(i)= gL*(Vm(i)-El);
    %% Summing all membrane currents (without Na+/K+ pump which is not included in model)
    Im(i)= Icurrent(i) - IK(i) - INa(i)- IL(i);
    %% iterating to find next membrane voltage using Eulers method
    Vm(i+1)= Vm(i) + dt*(Im(i)/Cm);
end
Vm(:)=Vm(:)-70;
%% plots membrane voltage
figure
plot(time,Vm)

ylabel('Membrane Voltage (mV)')
xlabel('Time(ms)')
%}

%% Plotting Conductances
figure
plot(time(2:end), gNaActual, 'b', time(2:end), gKActual,'r')
legend('Sodium Conductance','Potassium Conductance')
xlabel('Time (ms)')
ylabel('Conductance (mS/cm2)')

