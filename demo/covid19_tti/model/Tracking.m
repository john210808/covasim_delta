classdef Tracking
    properties
        %% Timeframe and initial conditions
        ti = [-40 0];
        tf = [0 120];
        tmax;
        tmin = -20;
        x0 = [0 20 0];
        dt = 0.005;
        
        % Default 
        xi= 0.15; 
        varphi = 0.2;
        eta_sc = 0.66;
        lambda_s = [0 0.1];
        Rt = 1.8;
        Rtlim = 2;
        eta;
        % xi; 
        xim;
        Gamma = 0.1;
        nu = 0.1;
        epsilon = 0.1;
        lambda_r = 0;
        Phi = 15;
        rt;
        nmax = 1e9;                     % Limitless tracing
        kmax;
        
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
        function o = Tracking()
            o.tmax = o.tf(end);
            o.eta = [0 o.eta_sc];
            o.xi = o.xi + (1-o.xi)*o.varphi; 
            o.xim = 1-o.xi;
            o.rt = [o.Rt o.Rt];
            o.kmax = length(o.ti);
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
        
        function r = calc(o)
            X = [];
            for j = 1:length(o.ti)
                if j>1
                    x0 = X(:,end)'; X(:,end)=[];
                end
                x = Solver.solver_por_tramos(o, j);
                if j>1
                    o.N_hat = [o.N_hat(1:end-1) ; o.newInfections(x, j)];
                    o.N_sum = [o.N_sum(1:end-1) ; o.newInfections_Total(x, j)];
                else
                    o.N_hat = o.newInfections(x, j);
                    o.N_sum = o.newInfections_Total(x, j);
                end
                X = [X x];
            end
            r = X;    
        end
        
        function N = newInfections(o, x, j)
            % obs X = [TSum HSum Hs];
            T = x(1,:);
            Ht = x(2,:);
            Hs = x(3,:);
            nbt0 = o.eta(j)*(o.lambda_s(j)*o.Rt*Hs + o.lambda_r*o.Rt*Ht);
            ne = zeros(size(T));
            for i = 1:length(nbt0)
                ne(i) = min(o.nmax,nbt0(i));
            end
            N = o.nu*o.Gamma*o.Rt*T + o.lambda_s(j)*Hs + o.lambda_r*Ht + ne;
            N = N';
        end
        
        function N = newInfections_Total(o, x, j)
            % obs X = [TSum HSum Hs];
            T = x(1,:);
            Ht = x(2,:);
            if length(o.Phi)>1
                N = (o.nu+o.epsilon)*o.Gamma*o.Rt*T + o.Gamma*o.Rt*Ht + o.Phi;
            else
                N = (o.nu+o.epsilon)*o.Gamma*o.Rt*T + o.Gamma*o.Rt*Ht + o.Phi;
            end
            N = N';
        end
    end
end