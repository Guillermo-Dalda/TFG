function val = quartic(x)
   val = 0;
   for i=1:length(x)
       val = val + x(i)^4;
   end
end