function precoder_reflectors(varargin)

if isempty(varargin)
    disp('using default simulation settings and parameters...');
    par.runId       = 0;        % simulation ID (used to reproduce results)
    par.plot        = false;     % plot results (true/false)
    par.save        = false;     % save results (true/false)
    par.BS          = 4;       % Antanna of BaseStation
    par.UE          = 4;       % Number of User Equipment
    par.order       = 4;        % PSK modulation order
    par.P           = 64;       % Number of reflectors
    par.P_list      = 512;
    par.BS_list     = 8:8:64;
    par.trials      = 100;      % number of Monte-Carlo trials (transmissions)
    par.precoder    = {'GD'};%'GD','AO-GD','AO-CVX','CI','CI-R'
    % select precoding scheme(s) to be evaluated
    par.s_fac       =  1;       % strong factor multiply on G
    par.kmax        =  5;       % max iteration of AO
    par.tol         =  1e-1;    % tollerence of iteration error
    par.rep_time    =  1e4;     % how many times of diferrent noise to on trial
else
    disp('use custom simulation settings and parameters...')
    par = varargin{1}; % load custom simulation parameters
end
BS_LIST_MOD = false;
P_LIST_MOD  = false;
for i=1:length(par.precoder)
    if(strcmp(par.precoder,'CI'))
        BS_LIST_MOD = true;
    else
        P_LIST_MOD  = true;
    end
end
if(BS_LIST_MOD && (length(par.precoder > 1)))
    error('BS_LIST_MOD is only supported for CI method')
elseif(BS_LIST_MOD && P_LIST_MOD)
    error('CI and other methods can not be used in one sim')
end

if(BS_LIST_MOD)
    L= length(par.BS_list);
else
    L= length(par.P_list);
end
% -- initialization
par.simName = ['TIME_',num2str(par.s_fac),'strong',num2str(par.UE),'Users_',num2str(par.BS),'antennas_and',num2str(par.P),'reflectors'...
    num2str(par.order),'-th_order_',num2str(par.runId),'_',datestr(clock,30)];
res.time = zeros(length(par.precoder),L);

% -- start simulation
profile on 
% track simulation time
for ll = 1:L
    if (BS_LIST_MOD)
        par.BS= par.BS_list(ll);
    else
        par.P  = par.P_list(ll);
    end
    for tt=1:par.trials
        u = randi([0 par.order-1],par.UE,1);
        Hd = sqrt(0.5)*(randn(par.UE,par.BS)+1i*randn(par.UE,par.BS));
        Hr = sqrt(0.5)*(randn(par.UE,par.P )+1i*randn(par.UE,par.P ));
        G  = sqrt(0.5*par.s_fac)*(randn(par.P ,par.BS)+1i*randn(par.P ,par.BS));
        for pp=1:length(par.precoder)
            % noise-dependent precoders
            switch (par.precoder{pp})
                case {'GD'}
                    %Soft max gradient descent method
                    V =eye(par.P);
                    x =randn(par.UE,1)+1j*randn(par.UE,1);
                    tic
                    [x,V]=soft_max_gradient_projection(Hr,V,G,Hd,x,u,par.order,1,200);
                    res.time(pp,ll)=res.time(pp,ll)+toc;
                case {'AO-CVX'}
                    %Soft max gradient descent method
                    V =eye(par.P);
                    x =randn(par.UE,1)+1j*randn(par.UE,1);
                    tic
                    [x,V] = A2_DAO(Hr,Hd,G,u,par.kmax,par.tol,x,par.order);
                    res.time(pp,ll)=res.time(pp,ll)+toc;
                case {'AO-GD'}
                    %Soft max gradient descent method
                    V =eye(par.P);
                    x =randn(par.UE,1)+1j*randn(par.UE,1);
                    tic
                    [x,V] = A2_DAO3(Hr,Hd,G,u,par.kmax,par.tol,x,par.order);
                    res.time(pp,ll)=res.time(pp,ll)+toc;
                case{'CI'}
                    tic
                    [x,V] = CI(Hr,Hd,G,u,par.order);
                    res.time(pp,ll)=res.time(pp,ll)+toc;
                case{'CI-R'}
                    tic
                    [x,V] = CI_R(Hr,Hd,G,u,par.order);
                    res.time(pp,ll)=res.time(pp,ll)+toc;
                otherwise
                    error('precoder not supported!')
            end
        end % algorithm loop
        
        % keep track of simulation time
    end % trials loop
end
    fprintf('\nsimulation has finished!\n\n');
    profile viewer
    profile off
    % -- save results
    
    if par.save
        save([par.simName],'par','res');
        %         save ('mx_test2.mat','par','res');
        if par.plot
            print('-depsc',[par.simName '.eps' ]);
        end
    end
end
