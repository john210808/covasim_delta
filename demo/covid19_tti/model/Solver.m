
classdef Solver
%     properties
% %         t
%     end

    methods (Static)
        function x = solver_por_tramos(o, j)
            t = o.ti:o.dt:o.tf;
            [~,x] = ode45(@(t,x) Solver.Pools_solver(o, j, t, x), t, o.x0);
            x = x';
            o.lambda_s(j)
        end
        
        function F = Pools_solver(o, j, t, x)
            % 1->T^sum; 2->H^sum; 3->H^s
            F = zeros(3,1); 
            xim = 1-o.xi;
            nbt0 = o.eta(j)*(o.lambda_r*o.Rt*x(2) + o.lambda_s(j)*o.Rt*x(3)); 
            if nbt0 >=o.nmax
                ne = o.nmax;
            else
                ne = nbt0;
            end
            F(1) =   o.Gamma*(o.nu*o.Rt-1)*x(1)         + o.lambda_r*x(2)                 + o.lambda_s(j)*x(3)         + ne;
            F(2) = o.Gamma*o.epsilon*o.Rt*x(1)        + (o.Gamma*(o.Rt-1)-o.lambda_r)*x(2)	- o.lambda_s(j)*x(3)         - ne + o.Phi;
            F(3) = o.xim*o.Gamma*o.epsilon*o.Rt*x(1)    + o.xim*o.Gamma*o.Rt*x(2)             -(o.lambda_s(j)+o.lambda_r+o.Gamma)*x(3)  + o.xim*(o.Phi-ne);
        end
    end
end