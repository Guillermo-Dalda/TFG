function val = powell_sum(x)
   val = 0;
   for i=1:length(x)
       val = val + abs(x(i))^(i+1);
   end
end