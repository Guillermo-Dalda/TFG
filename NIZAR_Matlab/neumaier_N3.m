function val = neumaier_N3(x)
   sum2 = 0;
   sum1 = (x(1) - 1) * (x(1) - 1);
   for i=2:length(x)
        sum1 = sum1 + (x(i) - 1) * (x(i) - 1);
        sum2 = sum2 + x(i) * x(i-1);
   end
   val = sum1 - sum2;
end