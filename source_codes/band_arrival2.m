function ts=band_arrival2(lambda_m,T,delta)
t=0;
i=2;
ts(1)=delta;
while t<T
   s=exprnd(1/lambda_m);  % poisson process, interval time
   t=s+t;
   if s+t<T
       ts(i,1)=t;
       i=i+1;
   else
   end
end

end
