function [raingauge_obs,raingauge_pos]=get_raingauge_obs(rain_field,X_range_sam,Y_range_sam,delta_x,delta_y,N_gauge,smooth_ratio)

x=X_range_sam*rand(N_gauge,1); % generate raingauge positions
y=Y_range_sam*rand(N_gauge,1);  

x_grid=ceil(x/delta_x); y_grid=ceil(y/delta_y); % specify which grid is the raingauge
raingauge_pos=[x y];
[~,~,N_t]=size(rain_field);
gauge_rain=zeros(N_t,N_gauge);
for i=1:N_gauge
gauge_rain(:,i)=squeeze(rain_field(x_grid(i),y_grid(i),:)); % get the rain data for each site
end
raingauge_obs=zeros(size(gauge_rain)); % define the output variable
N_t_smooth=N_t/smooth_ratio; % calculate the number of time steps after smooth

for i=1:N_t_smooth % generate the smoothed data
    temp=mean(gauge_rain(smooth_ratio*(i-1)+1:smooth_ratio*i,:),1);
    raingauge_obs(smooth_ratio*(i-1)+1:smooth_ratio*i,:)=repmat(temp,smooth_ratio,1);
end

end