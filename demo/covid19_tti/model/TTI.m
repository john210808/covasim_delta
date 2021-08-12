classdef TTI
    properties       
        params = []
        p
%         kmax;    
%         t, tmax, tmin, x0;
%         str = cell(1,kmax);
%         t = ti(1):dt:tf(end);
%         Xs = cell(kmax,1);
%         R0s = cell(kmax,1);
%         Nhat = cell(kmax,1);
%         Nhatobs = cell(kmax,1);
        N_hat = [];
        N_sum = [];
    end
    
    
    methods   
       function o = TTI()
           dt = 0.005;
           
           p1 = TTI.baseParams;
           p1.id = 1;
           p1.t = -40:dt:0;
           p1.eta = 0;
           p1.lambda_s = 0;
           
           p2 = TTI.baseParams;
           p2.id = 2;
           p2.t = 0:dt:120;
           p2.eta = 0.66;
           p2.lambda_s = 0.1;
           
           o.params = [p1 p2];
           
%             o.tmin = -20;
%             o.tmax = tf(end);
        
%             o.xi = o.xi + (1-o.xi)*o.varphi; 
%             o.rt = [o.Rt o.Rt];
%             o.kmax = length(ti);
        end
        
        % solve differential equations
        function x = doSolve(o)
            [~,x] = ode45(@(t, x) o.doSolve0(t, x), o.p.t, o.p.x0);
            x = x';
        end
        
        % autonomous differential equation, t is not used.
        function F = doSolve0(o, t, x)
            p = o.p;
            T = x(1); H = x(2); Hs = x(3);
            
            % tracing
            ne = p.eta*p.Rt*(p.lambda_s*Hs + p.lambda_r*H); 
            if ne >= p.nmax
                ne = p.nmax;
            end
            
            % differential equations, 
            dT = p.Gamma*(p.nu*p.Rt-1)*T + p.lambda_s*Hs + p.lambda_r*H + ne;
            dH = p.Gamma*(p.Rt-1)*H - (p.lambda_s*Hs + p.lambda_r*H) - ne + p.Gamma*p.epsilon*p.Rt*T + p.Phi;
            dHs = p.xim*p.Gamma*p.Rt*H - p.Gamma*Hs - (p.lambda_s+p.lambda_r)*Hs + p.xim*(-ne + p.Gamma*p.epsilon*p.Rt*T + p.Phi);
            
            F = zeros(3,1); 
            F = [dT; dH; dHs];
        end
        
        function r = calc(o)
            X = [];
            for j = 1:length(o.params)
                if j>1
                    x0 = X(:,end)'; X(:,end)=[];
                end
                o.p = o.params(j);
                x = o.doSolve();
%                 if j>1
%                     o.N_hat = [o.N_hat(1:end-1) ; o.newInfections(x)];
%                     o.N_sum = [o.N_sum(1:end-1) ; o.newInfections_Total(x)];
%                 else
%                     o.N_hat = o.newInfections(x);
%                     o.N_sum = o.newInfections_Total(x);
%                 end
                X = [X x];
            end
            r = X;    
        end
        
        function N = newInfections(o, x)
            p = o.p;
            % obs X = [TSum HSum Hs];
            T = x(1,:);
            Ht = x(2,:);
            Hs = x(3,:);
            nbt0 = p.eta(j)*(p.lambda_s(j)*p.Rt*Hs + p.lambda_r*p.Rt*Ht);
            ne = zeros(size(T));
            for i = 1:length(nbt0)
                ne(i) = min(p.nmax,nbt0(i));
            end
            N = p.nu*p.Gamma*p.Rt*T + p.lambda_s*Hs + p.lambda_r*Ht + ne;
            N = N';
        end
        
        function N = newInfections_Total(o, x)
            p = o.p;
            % obs X = [TSum HSum Hs];
            T = x(1,:);
            Ht = x(2,:);
            if length(p.Phi)>1
                N = (p.nu+p.epsilon)*p.Gamma*p.Rt*T + p.Gamma*p.Rt*Ht + p.Phi;
            else
                N = (p.nu+p.epsilon)*p.Gamma*p.Rt*T + p.Gamma*o.Rt*Ht + p.Phi;
            end
            N = N';
        end
        
        function r = doPlot()
            %%% Figure 2 %%%
            clear all
            close all
            clc
            
            %% Discretizacion

            idx = t==floor(t);
            T = t(idx);
            N_hat = N_hat(idx);
            N_sum = N_sum(idx);

            %% Calculo de los casos nuevos por d?a

            N_hat_obs = EstimDelay(N_hat,4,1,0.95);
            tau = 4;
            Rt_hat_obs = N_hat_obs./[NaN(tau,1) ; N_hat_obs(1:end-tau)];
            tau = 4;
            N_T = EstimDelay(N_sum,4,1,0.95);
            Rt_hat = N_T./[NaN(tau,1) ; N_T(1:end-tau)];

            %% Visualization

            load('Colores.mat')
            lt = '-';
            ls = '-.';
            lh = ':';

            fact_axis = 1;
            fact_label = 1.5;
            fact_curva = 1;
            siz = 15;
            green = [0 1 0]; 
            blue = [0 0 1];

            ylimsupI = 11000;
            ylimsupN = 1200;
            loc = 'northeast';
            %% Fig a

            figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
            ax = subplot(1,1,1);
            ax.Position = [0.25 0.25 0.65 0.65];
            ax.ActivePositionProperty = 'position';

            plot(t,x(1,:),lt,'Color',bl1,'LineWidth',3*fact_curva)
            hold on
            plot(t,x(1,:) + x(2,:),ls,'Color',bl2,'LineWidth',3*fact_curva)
            hold on
            plot(t,x(2,:),lh,'Color',bl3,'LineWidth',3*fact_curva)
            hold on
            legend({'Traced $T$','Total $T+H$','Hidden $H$'},'interpreter','latex','FontSize',15*fact_axis,'Location',loc);
            set(gca,'FontSize',15*fact_axis)
            hold on
            xlabel('days since TTI','interpreter','latex','FontSize',15*fact_label)
            ylabel('Active cases','interpreter','latex','FontSize',15*fact_label)
            ax.TickLabelInterpreter='latex';
            xlim([tmin tmax])
            ylim([0 ylimsupI])
            hold on
        end
        
    end
    
    
    methods (Static)
        function params=baseParams
            params.id = 0;
            % init values: T, H, Hs
            params.T = 0; 
            params.H = 20;
            params.Hs = 0;
            params.x0 = [params.T params.H params.Hs];

            params.Rt = 1.8;          % Reproduction number (hidden)
            params.Rtlim = 2;

            params.Gamma = 0.1;       % Recovery rate
            params.Phi = 15;          % Influx rate (hidden)
            params.xi= 0.15;          % Asymptomatic ratio
            params.xim = 1 - params.xi;

            params.varphi = 0.2;      % Fraction skipping testing
            params.lambda_s = 0;      % symptom-driven testing rate
            params.lambda_r = 0;      % random-testing rate

            params.nmax = 1e9;        % Limitless tracing
        %     params.eta_sc = 0.66;     % Tracing efficiency
            params.eta = 0;
            params.nu = 0.1;          % Isolation factor (traced)
            params.epsilon = 0.1;     % Missed contacts (traced) 

            params.ne = 0;
        end
    end 
end