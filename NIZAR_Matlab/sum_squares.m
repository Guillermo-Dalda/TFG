function val = sum_squares(x)
   val = 0;
   for i=1:length(x)
       val = val + i * x(i)^2;
   end
end