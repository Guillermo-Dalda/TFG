function val = ridge(x)
   val = 0;
   for i=2:length(x)
       val = val + x(i)^2;
   end
   val = x(1) + sqrt(val);
end