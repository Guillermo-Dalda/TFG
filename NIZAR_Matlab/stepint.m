function val = stepint(x)
   val = 0;
   for i=1:length(x)
       val = val + floor(x(i));
   end
   val = val + 25;
end