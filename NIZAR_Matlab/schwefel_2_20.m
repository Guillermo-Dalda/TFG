function val = schwefel_2_20(x)
   val = 0;
   for i=1:length(x)
       val = val + abs(x(i));
   end
end