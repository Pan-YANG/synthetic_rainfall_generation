function gauge_field=generate_gauge_field_inter(raingauge_obs,raingauge_pos,X_range_sam,Y_range_sam,x_step,y_step,t_step,method)
[N_t,N_gauge]=size(raingauge_obs); % calculate the number of time steps
gauge_field_obs=zeros(X_range_sam/x_step,Y_range_sam/y_step,N_t);
[xq,yq]=meshgrid(x_step/2:x_step:X_range_sam-x_step/2,y_step/2:y_step:Y_range_sam-y_step/2);
for i=1:N_t
    F=scatteredInterpolant(raingauge_pos(:,1),raingauge_pos(:,2),raingauge_obs(i,:)',method,'linear');
    gauge_field_obs(:,:,i)=F(xq,yq)';
end

for i=1:N_gauge
    gauge_field_obs(ceil(raingauge_pos(i,1)/x_step),ceil(raingauge_pos(i,2)/y_step),:)=raingauge_obs(:,i);
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