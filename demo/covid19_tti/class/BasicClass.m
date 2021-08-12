classdef BasicClass
   properties
      Value
      ti = [-40 0];
   end
   methods
        function obj = BasicClass(val)
            if nargin == 1
                obj.Value = val;
            end
        end
        function r = roundOff(obj)
         r = round([obj.Value],2);
        end
        function r = multiplyBy(obj,n)
         r = [obj.Value] * n;
        end
        function r = tt(obj)
            r = obj.roundOff() + 1;
        end
   end
   methods (Static)
       function r = testStatic(v)
           v + 1
       end
   end
end