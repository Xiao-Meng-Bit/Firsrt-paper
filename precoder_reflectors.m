function precoder_reflectors(varargin)

if isempty(varargin)
    disp('using default simulation settings and parameters...');
    par.runId       = 0;        % simulation ID (used to reproduce results)
    par.plot        = false;     % plot results (true/false)
    par.save        = true;    % save results (true/false)
    par.BS          = 4;       % Antanna of BaseStation
    par.UE          = 4;       % Number of User Equipment
    par.order       = 4;        % PSK modulation order
    par.P           = 64;       % Number of reflectors
    par.trials      = 1500;      % number of Monte-Carlo trials (transmissions)
    par.SNRdB_list  = -40:1:20;    % list of SNR [dB] values to be simulated
    par.precoder    = {'CI-R'};%'GD','AO-GD','AO-CVX','CI','CI-R'
    % select precoding scheme(s) to be evaluated
    par.s_fac       =  1;       % strong factor multiply on G
    par.kmax        =  5;       % max iteration of AO
    par.tol         =  1e-1;    % tollerence of iteration error
    par.rep_time    =  1e4;     % how many times of diferrent noise to on trial
else
    disp('use custom simulation settings and parameters...')
    par = varargin{1}; % load custom simulation parameters
end

% -- initialization
par.simName = ['BER_',num2str(par.s_fac),'strong',num2str(par.UE),'Users_',num2str(par.BS),'antennas_and',num2str(par.P),'reflectors'...
    num2str(par.order),'-th_order_',num2str(par.runId),'_',datestr(clock,30)];
% initialize result arrays
res.VER = zeros(length(par.precoder),length(par.SNRdB_list));
res.SER = zeros(length(par.precoder),length(par.SNRdB_list));
res.BER  = zeros(length(par.precoder),length(par.SNRdB_list));


% Tx and Rx power (average and max)
res.TxAvgPower = zeros(length(par.precoder),length(par.SNRdB_list));
res.RxAvgPower = zeros(length(par.precoder),length(par.SNRdB_list));
res.TxMaxPower = zeros(length(par.precoder),length(par.SNRdB_list));
res.RxMaxPower = zeros(length(par.precoder),length(par.SNRdB_list));

% -- start simulation

% track simulation time
time_elapsed = 0; tic;

for tt=1:par.trials
    u = randi([0 par.order-1],par.UE,1);
    u_rep = repmat(u,1,par.rep_time);
    Hd = sqrt(0.5)*(randn(par.UE,par.BS)+1i*randn(par.UE,par.BS));
    Hr = sqrt(0.5)*(randn(par.UE,par.P )+1i*randn(par.UE,par.P ));
    G  = sqrt(0.5*par.s_fac)*(randn(par.P ,par.BS)+1i*randn(par.P ,par.BS));
    for pp=1:length(par.precoder)
        % SNR loop
        
        % noise-dependent precoders
        switch (par.precoder{pp})
            case {'GD'}
                %Soft max gradient descent method
                V =eye(par.P);
                x =randn(par.UE,1)+1j*randn(par.UE,1);
                [x,V]=soft_max_gradient_projection(Hr,V,G,Hd,x,u,par.order,1);
            case {'AO-CVX'}
                %Soft max gradient descent method
                V =eye(par.P);
                x =randn(par.UE,1)+1j*randn(par.UE,1);
                [x,V] = A2_DAO(Hr,Hd,G,u,par.kmax,par.tol,x,par.order);
            case {'AO-GD'}
                %Soft max gradient descent method
                V =eye(par.P);
                x =randn(par.UE,1)+1j*randn(par.UE,1);
                [x,V] = A2_DAO3(Hr,Hd,G,u,par.kmax,par.tol,x,par.order);
            case{'CI'}
                [x,V] = CI(Hr,Hd,G,u,par.order);
            case{'CI-R'}
                [x,V] = CI_R(Hr,Hd,G,u,par.order);
            otherwise
                error('precoder not supported!')
        end
        Hx = Hr*V*G*x+Hd*x;
        for k=1:length(par.SNRdB_list)
            % set noise variance
            N0 = 10.^(-par.SNRdB_list(k)/10);
            n = sqrt(0.5)*(randn(par.UE,par.rep_time)+1i*randn(par.UE,par.rep_time));
            Hx_rep=repmat(Hx,1,par.rep_time);
            % transmit data over noisy channel
            y = Hx_rep(:) + sqrt(N0)*n(:);

            % extract maximum instantaneous transmitted and received power
            res.TxMaxPower(pp,k) = max(res.TxMaxPower(pp,k), sum(abs(x).^2));
            res.RxMaxPower(pp,k) = max(res.RxMaxPower(pp,k), sum(abs(Hx).^2)/par.UE);
            
            % extract average transmitted and received power
            res.TxAvgPower(pp,k) = res.TxAvgPower(pp,k) + sum(abs(x).^2);
            res.RxAvgPower(pp,k) = res.RxAvgPower(pp,k) + sum(abs(Hx).^2)/par.UE;
            
            % UE-side nearest-neighbor detection
            r = pskdemod(y,par.order,pi/par.order);
            u_bi = de2bi(u,par.order);
            u_bi_rep = repmat(u_bi,1,par.rep_time);
            r_bi = de2bi(r,par.order);
            
            % -- compute error metrics
            err = (u_rep(:)~=r); % check for symbol errors
            res.VER(pp,k) = res.VER(pp,k) + any(err)/par.rep_time; % vector error rate
            res.SER(pp,k) = res.SER(pp,k) + sum(err)/par.UE/par.rep_time; % symbol error rate
            res.BER(pp,k) = res.BER(pp,k) + ...
                sum(sum(u_bi_rep(:)~=r_bi(:)))/(par.UE*log2(par.order)*par.rep_time); % bit error rate
            
        end % SNR loop
    end % algorithm loop

% keep track of simulation time
if toc>10
    
    time=toc;
    time_elapsed = time_elapsed + time;
    fprintf('estimated remaining simulation time: %3.0f min.\n',time_elapsed*(par.trials/tt-1)/60);
    tic;
    
end
end % trials loop
% normalize results
res.VER = res.VER/par.trials;
res.SER = res.SER/par.trials;
res.BER = res.BER/par.trials;
res.TxAvgPower = res.TxAvgPower/par.trials;
res.RxAvgPower = res.RxAvgPower/par.trials;
res.time_elapsed = time_elapsed;

fprintf('\nsimulation has finished!\n\n');

% -- save results

if par.save
    save([par.simName],'par','res');
    %         save ('mx_test2.mat','par','res');
    if par.plot
        print('-depsc',[par.simName '.eps' ]);
    end
end
end