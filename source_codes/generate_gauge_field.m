function gauge_field=generate_gauge_field(raingauge_obs,raingauge_pos,X_range_sam,Y_range_sam,x_step,y_step,t_step)
[N_t,~]=size(raingauge_obs); % calculate the number of time steps
gauge_field_obs=zeros(X_range_sam/x_step,Y_range_sam/y_step,N_t);
for i=1:X_range_sam/x_step
   for j=1:Y_range_sam/y_step 
      [~,ind_theissen]=min((raingauge_pos(:,1)-(i-0.5)*x_step).^2+(raingauge_pos(:,2)-(j-0.5)*x_step).^2);
% identify the closet gauge of a grid point, ind_theissen provides the index
% of the gauge
      gauge_field_obs(i,j,:)=raingauge_obs(:,ind_theissen);
% the gauge field at the specified grid point is then assigned the value of 
% selected raingauge
   end
end
% calculate the temperoal average if the time step (t_step) to report is not the
% same as the time step to generate the rainfield
if t_step~=1
    gauge_field=zeros(X_range_sam/x_step,Y_range_sam/y_step,N_t/t_step);
for i=1:N_t/t_step
    gauge_field(:,:,i)=mean(gauge_field_obs(:,:,(i-1)*t_step+1:i*t_step),3);
end
else
    gauge_field=gauge_field_obs;
end

end