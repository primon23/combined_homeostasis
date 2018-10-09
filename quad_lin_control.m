% quad_lin_control.m
%
% Code to produce a dual controller one quadratic, one linear, for
% threshold and excitatory synaptic conductance.
%
%   This code is used to produce Figure 3 in the paper,
%   "Combined mechanisms of neural firing rate homeostasis: a mathematical 
%   primer." 
%   by Paul Miller, Biological Cybernetics 2018.
%
%   It is important to note that the failure of dual homeostasis depicted
%   in Figure 3B is the result of switching the response functions, not the
%   goal rates between the threshold feedback and the conductance feedback 
%   functions -- i.e. the caption is incorrect. If goal rates are switched,
%   dual homeostasis still fails (conductance becomes zero and rate becomes
%   constant at the threshold-goal) but that is not the plotted result in
%   the paper.
%
% Tests dual homeostasis in response to changes of input current.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

r_t_goal = 10;      % Goal rate for threshold control mechanism
r_g_goal = 15;      % Goal rate for synaptic conductance control mechanism
r_t_power = 1;      % The rate is raised to this power in the threshold-control function
r_g_power = 2;      % The rate is raised to this power in the conductance-control function

epsilon_g = 0.002;  % Inverse time constant for homeostatic control of g
epsilon_t = 0.002;  % Inverse time constant for homeostatic control of T

% Set of parameters for the firing rate model sigmoid function follow where
% r = rmax/(1+exp( (rthresh-g.I)/rsigma ) )
rthresh0 = 50;      % Initial value of rthresh
rsigma = 20;        % Input range that causes rate variation
rmax = 100;         % Maximum rate

tmax = 3000000;     % Maximum simulation time
dt = 0.1;           % Time step in simulations
tvec = dt:dt:tmax;  % Vector of time points
Nt = length(tvec);  % Number of time points

r = zeros(size(tvec));  % Vector to store firing rates

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

%% Set up the vectors for conductance and threshold, which change in time
g = zeros(size(tvec));
g(1) = 1;
rthresh = zeros(size(tvec));
rthresh(1) = rthresh0;

% Do two simulations, the first with functioning dual homeostasis, the
% second with a switching of the exponents of rate in the feedback functions.
for simulation_number = 1:2
    if ( simulation_number == 2 )
    % Important note, the following two lines cause panel B of the figure
    % to be the result of a switching of the exponent/power of the
    % rate-dependence of the threshold feedback and the conductance
    % feedback, such that threshold feedback becomes more supralinear than
    % the conductance feedback. This causes "wind-up" and is the case
    % plotted in Figure 3B in the paper (the caption is incorrect).
         r_t_power = 2;
         r_g_power = 1;
         
         % To switch goal rates, use the following lines instead. In this
         % case, the conductance becomes zero and the variance of firing
         % rate also becomes zero, with rate at a constant of r_t_goal =
         % 15Hz. Clearly such a result is a failure of dual homeostasis, if
         % not exactly wind-up because conductance is constrained to remain
         % non-negative.
        %r_t_goal = 15;
        %r_g_goal = 10;
    end
    
    %% Simulate through all time-points
    for i = 2:Nt
        % Firing rate is a sigmoidal function of input scaled by g and shifted
        % by the threshold.
        r(i) = rmax/(1+exp(-(g(i-1)*I(i-1)-rthresh(i-1))/rsigma));
        
        % Now update g based on the new firing rate using the slow homeostatic
        % control function
        g(i) = g(i-1) - epsilon_g*( (r(i)/r_g_goal)^r_g_power - 1);
        if ( g(i) < 0 )
            g(i) = 0;       % Conductance can not be negative
        end
        
        % Now update rthresh based on the new firing rate using the slow homeostatic
        % control function
        rthresh(i) = rthresh(i-1) + epsilon_t*((r(i)/r_t_goal)^r_t_power - 1);
        
    end
    
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
    step = 60;
    rs = r(step:step:end);
    rs1 = smooth(rs,10);
    Is = I(step:step:end);
    Is1 = smooth(Is,10);
    gs = g(step:step:end);
    gs1 = smooth(gs,10);
    ts = rthresh(step:step:end);
    ts1 = smooth(ts,10);
    
    % Divide time into block of size averange = 5mins (300s)
    averange = 300/dt;
    mr = zeros(1,round(Nt/averange));   % Will contain mean rate in each block
    vr = zeros(1,round(Nt/averange));   % Will contain variance of rate in each block
    % Next loop through each block of time and record both mean and variance of
    % rate.
    for i = 1:round(Nt/averange)
        mr(i) = mean(r((i-1)*averange+1:i*averange));
        vr(i) = var(r((i-1)*averange+1:i*averange));
    end
    
    % Plot the final figure with multiple panels.
    column_posn = 0.12 + 0.5*(simulation_number-1);
    figure(4)
    % Top panel plot input current
    subplot('Position',[column_posn 0.83 0.35 0.15])
    plot(Is1(10:end-10),'k')
    xlabel('Time')
    set(gca,'XTick',[])
    ylabel('Applied current')
    
    % Second panel plot mean rate
    subplot('Position',[column_posn 0.635 0.35 0.15])
    plot(mr,'k')
    xlabel('Time')
    set(gca,'XTick',[])
    ylabel('Mean rate (Hz)')
    % Third panel plot variance in rate
    subplot('Position',[column_posn 0.44 0.35 0.15])
    plot(vr,'k')
    xlabel('Time')
    set(gca,'XTick',[])
    ylabel('Var(rate) (Hz^{2})')
    
    % Fourth panel plot synaptic conductance
    subplot('Position',[column_posn 0.245 0.35 0.15])
    plot(gs1(10:end-10),'k')
    xlabel('Time')
    set(gca,'XTick',[])
    ylabel('Conductance, g')
    
    % Fifth panel plot threshold
    subplot('Position',[column_posn 0.05 0.35 0.15])
    plot(ts1(10:end-10),'k')
    xlabel('Time')
    set(gca,'XTick',[])
    ylabel('Threshold, T')
end

% Label panels
annotation('textbox',[0.0 0.96 0.04 0.04],'String','A1','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.0 0.77 0.04 0.04],'String','A2','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.0 0.58 0.04 0.04],'String','A3','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.0 0.39 0.04 0.04],'String','A4','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.0 0.20 0.04 0.04],'String','A5','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.5 0.96 0.04 0.04],'String','B1','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.5 0.77 0.04 0.04],'String','B2','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.5 0.58 0.04 0.04],'String','B3','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.5 0.39 0.04 0.04],'String','B4','LineStyle','none','FontSize',16,'FontWeight','Bold')
annotation('textbox',[0.5 0.20 0.04 0.04],'String','B5','LineStyle','none','FontSize',16,'FontWeight','Bold')
