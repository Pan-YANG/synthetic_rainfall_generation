function CPC=cluster_potential_larger_border(rou_L,X,Y,Ub,alpha,D,sigma)
% cluster potential function that produces CPC in a larger domain
% when generating the CPC in the middle point, it considers the wind speed
X_wind=Ub(1)/alpha;
Y_wind=Ub(2)/alpha;
lambda=rou_L*(X+7*X_wind+4*D+4*sigma(1))*(Y+7*Y_wind+4*D+4*sigma(2));
npoints=poissrnd(lambda);
if npoints==1
    CPC(1,1)=X/2+Ub(1)*1/alpha; CPC(1,2)=Y/2+Ub(2)*1/alpha;
else
  CPC(:,1)= (X+7*X_wind+4*D+4*sigma(1))*...
      (rand(npoints,1)-((4*X_wind+2*D+2*sigma(1)))/(X+7*X_wind+4*D+4*sigma(1)));
  CPC(:,2)= (Y+7*Y_wind+4*D+4*sigma(2))*...
      (rand(npoints,1)-((4*Y_wind+2*D+2*sigma(2)))/(Y+7*Y_wind+4*D+4*sigma(2)));
  CPC(1,1)=X/2+Ub(1)*1/alpha; CPC(1,2)=Y/2+Ub(2)*1/alpha;
end
end