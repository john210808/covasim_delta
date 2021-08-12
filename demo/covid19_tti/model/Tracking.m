classdef Tracking
    properties
        params
    end
    
    methods
        function o = Tracking()
            o.params.Rt = 1.8;          % Reproduction number (hidden)
            o.params.Rtlim = 2;

            o.params.Gamma = 0.1;       % Recovery rate
            o.params.Phi = 15;          % Influx rate (hidden)
            o.params.xi= 0.15;          % Asymptomatic ratio

            o.params.varphi = 0.2;      % Fraction skipping testing
            o.params.lambda_s = 0;% symptom-driven testing rate
            o.params.lambda_r = 0;      % random-testing rate

            o.params.nmax = 1e9;        % Limitless tracing
            o.params.eta_sc = 0.66;     % Tracing efficiency
%             p.eta;
            o.params.nu = 0.1;          % Isolation factor (traced)
            o.params.epsilon = 0.1;     % Missed contacts (traced) 
        end
    end
end