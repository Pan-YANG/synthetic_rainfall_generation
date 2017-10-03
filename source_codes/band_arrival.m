function ts=band_arrival(lambda_m,T,delta,cell_dur)
% the start time should be one cell duration earlier than the first arrival
% of the rain band
% the defaul value of cell_dur is 0, which will allow the function to be
% used in other applications, e.g., the time at which crowd observations
% are taken
if nargin==3
    cell_dur=0;
end
i=2;
ts(1)=delta+ceil(cell_dur/delta)*delta;
t=ts(1);
while t<T
   s=exprnd(1/lambda_m);  % poisson process, interval time
   s=round(s/delta)*delta;
   t=s+t;
   if s+t<T
       ts(i,1)=t;
       i=i+1;
   else
   end
end

end
