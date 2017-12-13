% intracellular_homeostasis.m
%
%   Tests how dual homeostasis works with a branching biochemical pathway
%   and an integration step either before or after the branch-point.
%
%   This code is used to produce Figure 6 in the paper,
%   "Combined mechanisms of neural firing rate homeostasis: a mathematical 
%   primer." 
%   by Paul Miller, Biological Cybernetics 2018.
%
%   The system of equations are such that
%   Firing rate depends on inputs, I, excitatory conductance, g, and
%   threshold, T.
%   Calcium, x1 depends on firing-rate, r
%   CaMKIV activation, x2, depends on calcium (x1)
%   mRNA for negative control of excitatory synaptic conductance, x3a, 
%   depends on CaMKIV (x2)
%   mRNA for channels that alter threshold, x3b, depends on CaMKIV (x2)
%   excitatory conductance, g, depends negatively on x3a
%   threshold, T, depends on x3b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

tmax = 3000000;         % Maximum time to simulate (sec)
dt = 0.1;               % Time-step for simulations (sec)
tvec = dt:dt:tmax;      % Vector of time-points
Nt = length(tvec);      % Number of time-points

r = zeros(size(tvec));  % Initialize firing-rate vector
%% Set up the input vector
I = zeros(size(tvec));

Imean = 2.0;        % Baseline mean input
Isigma = 2.0;       % Baseline standard deviation of input

% Input vector will be split into four sections, with abrupt changes in the
% distribution at distinct points in time
% First 1/4 of simulation with baseline inputs:
I(1:round(Nt/4)) = Imean + Isigma*randn(1,round(Nt/4));
% Second 1/4 of simulation with raised mean input, baseline variance
I(round(Nt/4)+1:round(Nt/2)) = 5*Imean + Isigma*randn(1,round(Nt/4));
% Third quarter of simulation with baseline mean, raised variance of input
I(round(Nt/2)+1:round(3*Nt/4)) = Imean + 5*Isigma*randn(1,round(Nt/4));
% Final quarter of simulation with baseline mean, reduced variance of input
I(round(3*Nt/4)+1:round(Nt)) = Imean + 0.2*Isigma*randn(1,round(Nt/4));

%% Define all variables
x1 = zeros(size(tvec));     % x1 represents calcium
x2 = zeros(size(tvec));     % x2 represents CaMKIV
x3a = zeros(size(tvec));    % x3a represents mRNA for excitatory synaptic conductance
x3b = zeros(size(tvec));    % x3b represents mRNA for channels that affect threshold (eg somatic potassium channels)

g = zeros(size(tvec));      % excitatory synaptic conductance
rthresh = zeros(size(tvec));    % firing rate threshold

%% Set up the parameters
% First the set of time constants for each reaction
tau1 = 0.1;             % rate to calcium
tau2 = 2;               % calcium to CaMKIV activation
tau3a = 1800;           % CaMKIV to mRNA (g)
tau3b = 1800;           % CaMKIV to mRNA (T)
tau_g = 300;            % mRNA (g) to synaptic conductance
tau_th = 30;            % mRNA (T) to threshold change

% Next set up the relative rate constants for degradation of each variable
% as a first-order process.
% Any rate constant set to zero here means the corresponding variable
% integrates its inputs.
alpha2 = 1;         % CaMKIV activation
alpha3a = 0.0;      % mRNA for g
alpha3b = 0;        % mRNA for T

% For the dual controller, g_exponent = 2 ensures the control variable for 
% synaptic conductance integrates a quadratic function of firing rate
g_exponent = 2;     
gmax = 50;          % Maximum synaptic conductance
g_init = 2;         % Initial synaptic conductance

%% The next three quantities define the input-output curve of the neuron. 
%   These are the paramters of the logistic function:
%   r = rmax/(1+exp( (rthresh-input)/rsigma ) )
rthresh0 = 20;      % Baseline (initial) input for 1/2-maximum firing
rsigma = 20;        % Defines steepness of f-I curve (constant)
rmax = 100;         % Maximum firing rate of neuron

x2_loss = 2;        % Baseline constant rate of loss of active CaMKII
x3a_loss = 400;     % Baseline constant rate of loss of x3a
x3b_loss = 18;
x1bar = x3b_loss;
    
% x3a_scale is used to scale down x3a when converting to conductance
% change.
x3a_scale = 1000;

% Do two simulations, the first with functioning dual homeostasis, the
% second with a single controller/integrator (CaMKIV = x2) before the 
% branch-point.
for simulation_number = 1:2
    if ( simulation_number == 2 )
        x2_loss = 20;       % constant rate of loss of x2
        x3a_loss = 0;       % remove constant rate of loss of x3a
        x3b_loss = 0;       % remove constant rate of loss of x3b
        alpha2 = 0;         % remove proportionality to x2 in loss of x2
        alpha3a = 4;        % proportionality to x3a in rate of loss of x3a
        alpha3b = 4;        % proportionality to x3b in rate of loss of x3b
        g_exponent = 1;     % no need for quadratic dependence with one controller
        tau2 = 30;          % dx2/dt time constant (sec)
        tau3a = 1800;       % dx3a/dt time constant (sec)
        tau3b = 1800;       % dx3b/dt time constant (sec)
        
    end
    
    g(1) = g_init;              % Initial synaptic conductance
    rthresh(1) = rthresh0;      % Initial input for 1/2-max firing rate
    x3a(1) = x3a_scale-g(1);    % Initial value for x3a
    x3b(1) = rthresh(1);        % Initial value for x3b 
    x2(1) = r(1)-x2_loss;       % Initial value for x2 (CaMKII)
    
    
    %% Simulate through all time-points
    t_out = round(Nt/100);      % Used to calculate progress every 1% of time points
    for i = 2:Nt;
        
        % The next few lines simply display the progress as % of total time
        if ( mod(i,t_out) == 0 )
            progress = i/t_out;
            disp(strcat(num2str(progress),'% complete'));
        end
        
        % Immediate update of firing rate with logistic f-I curve
        r(i) = rmax/(1+exp(-(g(i-1)*I(i-1)-rthresh(i-1))/rsigma));
        
        % x1 is e.g. somatic calcium, a linear filter of firing rate
        x1(i) = x1(i-1) + dt/tau1*(r(i-1) - x1(i-1));
        x1(i) = max(x1(i),0);
        
        % x2 is e.g. CaMKII activation.
        % In method-1 a linear filter of calcium
        % In method-2 an integration of calcium (alpha2=0)
        x2(i) = x2(i-1) + dt/tau2*(x1bar*(x1(i-1)/x1bar)^3 - alpha2*x2(i-1) - x2_loss);
        x2(i) = max(x2(i),0);
        
        % x3a is a factor which decreases synaptic conductance
        % In method-1 it is an integrator of activated CaMKII (alpha3a=0)
        % In method-2 it is a linear filter of activated CaMKII
        x3a(i) = x3a(i-1) + dt/tau3a*(x2(i-1)^g_exponent - alpha3a*x3a(i-1) - x3a_loss);
        x3a(i) = max(x3a(i),0);
        
        % x3b is a factor which decreases intrinsic excitability 
        % In method-1 it is an integrator of activated CaMKII (alpha3a=0)
        % In method-2 it is a linear filter of activated CaMKII      
        x3b(i) = x3b(i-1) + dt/tau3b*(x2(i-1) - alpha3b*x3b(i-1) - x3b_loss);
        x3b(i) = max(x3b(i),0);
        
        % Synaptic conductance depends negatively on x3a
        g(i) = g(i-1) + dt*( gmax*(1 - x3a(i-1)/x3a_scale) - g(i-1) )/tau_g;
        if ( g(i) < 0 )
            g(i) = 0;
        end
        
        % Threshold (input for 1/2-maximum rate) increases with x3b
        rthresh(i) = rthresh(i-1) + dt*(x3b(i-1) - rthresh(i-1) )/tau_th;
        
    end
    
    %% Set default styles for the plot
    set(0,'DefaultLineLineWidth',2,...
        'DefaultLineMarkerSize',8, ...
        'DefaultAxesLineWidth',2, ...
        'DefaultAxesFontSize',14,...
        'DefaultAxesFontWeight','Bold');
    
    
    %% Display some results indicating whether homeostasis was achieved
    % In the second half of each period of constant input current distribution,
    % find the mean rate.
    m1 = mean(r(round(Nt/8)+1:round(Nt/4)))
    m2 = mean(r(round(3*Nt/8)+1:round(Nt/2)))
    m3 = mean(r(round(5*Nt/8)+1:round(3*Nt/4)))
    m4 = mean(r(round(7*Nt/8)+1:round(Nt)))
    
    % In the second half of each period of constant input current distribution,
    % find the standard deviation of the rate.
    s1 = std(r(round(Nt/8)+1:round(Nt/4)))
    s2 = std(r(round(3*Nt/8)+1:round(Nt/2)))
    s3 = std(r(round(5*Nt/8)+1:round(3*Nt/4)))
    s4 = std(r(round(7*Nt/8)+1:round(Nt)))
    
    % In the second half of each period of constant input current distribution,
    % find the correlation between input and response.
    corrtest1 = corr(I(round(Nt/8)+1:round(Nt/4)-1)',r(round(Nt/8)+2:round(Nt/4))')
    corrtest2 = corr(I(round(3*Nt/8)+1:round(Nt/2)-1)',r(round(3*Nt/8)+2:round(Nt/2))')
    corrtest3 = corr(I(round(5*Nt/8)+1:round(3*Nt/4)-1)',r(round(5*Nt/8)+2:round(3*Nt/4))')
    corrtest4 = corr(I(round(7*Nt/8)+1:round(Nt)-1)',r(round(7*Nt/8)+2:round(Nt))')
    
    %% Now plot all the results
    
    % Set default styles for the plot
    set(0,'DefaultLineLineWidth',2,...
        'DefaultLineMarkerSize',8, ...
        'DefaultAxesLineWidth',2, ...
        'DefaultAxesFontSize',14,...
        'DefaultAxesFontWeight','Bold');
    
    % Calculate smoothed versions of the variables so that outliers do not
    % dominate the visibility of the fluctuations on a graph.
    step = 600;
    rs = r(step:step:end);
    rs1 = smooth(rs,10);
    Is = I(step:step:end);
    Is1 = smooth(Is,10);
    gs = g(step:step:end);
    gs1 = smooth(gs,10);
    ts = rthresh(step:step:end);
    ts1 = smooth(ts,10);
    x2s = x2(step:step:end);
    x2s1 = smooth(x2s,10);
    
    % Divide time into block of size averange = 5mins (300s)
    averange = 300/dt;
    mr = zeros(1,round(Nt/averange));   % Will contain mean rate in each block
    mr3 = zeros(1,round(Nt/averange));   % Will contain mean rate in each block
    vr = zeros(1,round(Nt/averange));   % Will contain variance of rate in each block
    % Next loop through each block of time and record both mean and variance of
    % rate.
    for i = 1:round(Nt/averange)
        mr(i) = mean(r((i-1)*averange+1:i*averange));
        mr3(i) = mean(r((i-1)*averange+1:i*averange).^3);
        vr(i) = var(r((i-1)*averange+1:i*averange));
    end
    
    % Plot the final figure with multiple panels.
    column_posn = 0.12 + 0.5*(simulation_number-1);
    
    figure(4)
    subplot('Position',[column_posn 0.855 0.35 0.12])
    plot(Is1(10:end-10),'k')
    xlabel('Time')
    set(gca,'XTick',[])
    ylabel('Applied current')
    subplot('Position',[column_posn 0.69 0.35 0.12])
    plot(mr,'k')
    xlabel('Time')
    set(gca,'XTick',[])
    ylabel('Mean rate (Hz)')
    subplot('Position',[column_posn 0.525 0.35 0.12])
    plot(vr,'k')
    xlabel('Time')
    set(gca,'XTick',[])
    ylabel('Var(rate) (Hz^{2})')
    subplot('Position',[column_posn 0.36 0.35 0.12])
    plot(x2s1(10:end-10),'k')
    xlabel('Time')
    set(gca,'XTick',[])
    ylabel('Active Kinase')
    subplot('Position',[column_posn 0.195 0.35 0.12])
    plot(gs1(10:end-10),'k')
    xlabel('Time')
    set(gca,'XTick',[])
    ylabel('Synaptic g')
    subplot('Position',[column_posn 0.03 0.35 0.12])
    plot(ts1(10:end-10),'k')
    xlabel('Time')
    set(gca,'XTick',[])
    ylabel('Threshold, T')
    
end

% Label all panels
annotation('textbox',[0 0.96 0.04 0.04],'String','A1','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.795 0.04 0.04],'String','A2','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.63 0.04 0.04],'String','A3','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.465 0.04 0.04],'String','A4','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.30 0.04 0.04],'String','A5','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0 0.135 0.04 0.04],'String','A6','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.5 0.96 0.04 0.04],'String','B1','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.5 0.795 0.04 0.04],'String','B2','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.5 0.63 0.04 0.04],'String','B3','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.5 0.465 0.04 0.04],'String','B4','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.5 0.30 0.04 0.04],'String','B5','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.5 0.135 0.04 0.04],'String','B6','LineStyle','none','FontSize',16,'FontWeight','Bold')

