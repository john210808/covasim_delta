classdef TTI < handle
    properties       
        params = [];
        p;
        X = [];
        N_hat = [];
        N_sum = [];
    end
    
    methods (Static)
        function params=baseParams
            params.id = 0;
            % init values: T, H, Hs
            params.x0 = [0 20 0];

            params.Rt = 1.8;          % Reproduction number (hidden)
            params.Rtlim = 2;

            params.Gamma = 0.1;       % Recovery rate
            params.Phi = 15;          % Influx rate (hidden)
            params.xi= 0.15;          % Asymptomatic ratio
            
            params.varphi = 0.2;      % Fraction skipping testing
            params.lambda_s = 0;      % symptom-driven testing rate
            params.lambda_r = 0;      % random-testing rates
            
            params.xi =  params.xi + (1-params.xi)*params.varphi;
            params.xim = 1 - params.xi;
            
            params.nmax = 1e9;        % Limitless tracing
        %     params.eta_sc = 0.66;     % Tracing efficiency
            params.eta = 0;
            params.nu = 0.1;          % Isolation factor (traced)
            params.epsilon = 0.1;     % Missed contacts (traced) 
        end
        
        function r = Run(eta)
            TTI(eta).calc().doPlot();
        end
    end 
    
    methods   
       function this = TTI(eta)
            dt = 0.005;

            p1 = TTI.baseParams;
            p1.id = 1;
            p1.t = -40:dt:0;
            p1.eta = 0;
            p1.lambda_s = 0;

            p2 = TTI.baseParams;
            p2.id = 2;
            p2.t = 0:dt:120;
            p2.eta = eta; % 0.66 or 0.33
            p2.lambda_s = 0.1;
           
            this.params = [p1 p2];
            this.p;
            this.X = [];
            this.N_hat = [];
            this.N_sum = [];
        end
        
        % solve differential equations
        function x = doSolve(this, p)
            [~,x] = ode45(@(t, x) this.doSolve0(t, x), p.t, p.x0);
            x = x';
        end
        
        % autonomous differential equation, t is not used.
        function x = doSolve0(this, t, x)
            p = this.p;
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
            
            x = zeros(3,1); 
            x = [dT; dH; dHs];
        end
        
        function this = calc(this)
            X = [];
            for j = 1:length(this.params)
                this.p = this.params(j);
                if j>1
                    this.p.x0 = X(:,end)'; X(:,end)=[];
                end
                x = this.doSolve(this.p);
                if j>1
                    this.N_hat = [this.N_hat(1:end-1) ; this.newInfections(x)];
                    this.N_sum = [this.N_sum(1:end-1) ; this.newInfections_Total(x)];
                else
                    this.N_hat = this.newInfections(x);
                    this.N_sum = this.newInfections_Total(x);
                end
                X = [X x];
            end
            this.X = X;    
        end
        
        function N = newInfections(this, x)
            p = this.p;
            % obs X = [TSum HSum Hs];
            T = x(1,:);
            Ht = x(2,:);
            Hs = x(3,:);
            nbt0 = p.eta*(p.lambda_s*p.Rt*Hs + p.lambda_r*p.Rt*Ht);
            ne = zeros(size(T));
            for i = 1:length(nbt0)
                ne(i) = min(p.nmax,nbt0(i));
            end
            N = p.nu*p.Gamma*p.Rt*T + p.lambda_s*Hs + p.lambda_r*Ht + ne;
            N = N';
        end
        
        function N = newInfections_Total(this, x)
            p = this.p;
            % obs X = [TSum HSum Hs];
            T = x(1,:);
            Ht = x(2,:);
            if length(p.Phi)>1
                N = (p.nu+p.epsilon)*p.Gamma*p.Rt*T + p.Gamma*p.Rt*Ht + p.Phi;
            else
                N = (p.nu+p.epsilon)*p.Gamma*p.Rt*T + p.Gamma*p.Rt*Ht + p.Phi;
            end
            N = N';
        end
        
        function r = doPlot(this)
            %% Discretizacion
            t = -40:0.005:120;
            idx = t==floor(t);
            T = t(idx);
            this.N_hat = this.N_hat(idx);
            this.N_sum = this.N_sum(idx);
%             r = 0

            %% Calculo de los casos nuevos por d?a

            N_hat_obs = EstimDelay(this.N_hat,4,1,0.95);
            tau = 4;
            Rt_hat_obs = N_hat_obs./[NaN(tau,1) ; N_hat_obs(1:end-tau)];
            tau = 4;
            N_T = EstimDelay(this.N_sum,4,1,0.95);
            Rt_hat = N_T./[NaN(tau,1) ; N_T(1:end-tau)];

            %% Visualization

            load('Colores.mat')
            lt = '-';
            ls = '-.';
            lh = ':';

            fact_axis = 1;
            fact_label = 1;
            fact_curva = 1;
            siz = 5;
            green = [0 1 0]; 
            blue = [0 0 1];

            ylimsupI = 11000;
            ylimsupN = 1200;
            loc = 'northeast';
            %% Fig a

            tmin = -20;
            tmax = 120;
            x = this.X;
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
            
            
            
            %% Fig b

            figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
            ax = subplot(1,1,1);
            ax.Position = [0.25 0.25 0.65 0.65];
            ax.ActivePositionProperty = 'position';

            plot(T,N_hat_obs,lt,'Color',or1,'LineWidth',3*fact_curva)
            hold on
            plot(T,this.N_sum,ls,'Color',or2,'LineWidth',3*fact_curva)
            hold on
            legend({'$\hat{N}^{\mbox{obs}}$','$N$'},'interpreter','latex','FontSize',15*fact_axis);
            set(gca,'FontSize',15*fact_axis)
            hold on
            xlabel('days since TTI','interpreter','latex','FontSize',15*fact_label)
            ylabel('New infections','interpreter','latex','FontSize',15*fact_label)
            hold on
            xlim([tmin tmax])
            ylim([0 ylimsupN])
            ax.TickLabelInterpreter='latex';


            Rtlim = 2;
            %% Subfig c

            figure('units','normalized','position',[0.3 0.3 0.5 0.6]);
            ax = subplot(1,1,1);
            ax.Position = [0.25 0.25 0.65 0.65];
            ax.ActivePositionProperty = 'position';


            plot(T,Rt_hat_obs,lt,'Color',red1,'LineWidth',3*fact_curva)
            hold on
            plot(T,Rt_hat,ls,'Color',red2,'LineWidth',3*fact_curva)
            hold on
            plot(tmin:1:tmax,ones(size(tmin:1:tmax)),'k--','LineWidth',1*fact_curva,'HandleVisibility','off');
            hold on
            legend({'$\hat{R}_t^{\mbox{obs}}$','$\hat{R}_t$'},'interpreter','latex','FontSize',15*fact_axis);
            set(gca,'FontSize',15*fact_axis)
            hold on
            xlabel('days since TTI','interpreter','latex','FontSize',15*fact_label)
            ylabel('Rep. number','interpreter','latex','FontSize',15*fact_label)
            hold on
            ylim([0.75 Rtlim])
            xlim([tmin tmax])

            ax.TickLabelInterpreter='latex';
        end
    end
end