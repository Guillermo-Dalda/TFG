function val = sphere(x)
   val = 0;
   for i=1:length(x)
       val = val + x(i)^2;
   end
end