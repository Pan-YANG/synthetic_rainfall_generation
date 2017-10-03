function [crowd_field,N_obs_c,rep_bias]=generate_crowd_field_withfixedpoint_2dinter(rain_field,bias,lambda_c,pop_den,fixed_point_den,X_range_sam,Y_range_sam,delta_x,delta_y,sigma_x,sigma_y,T,delta,sigma_obs,t_step,x_step,y_step)
% bias is calculated as a ratio of the ground truth rainfall
% fixed_point_den is the density of fixed crowd monitoring facilities 
[N_x_grid,N_y_grid,N_t]=size(rain_field);
lambda=lambda_c*pop_den*X_range_sam*Y_range_sam; % calculate the poisson rate 
% lambda based on the participation rate lambda_c, pop_den and spatial
% extents X_range_sam and Y_range_sam
N_fixed_point=floor(fixed_point_den*X_range_sam*Y_range_sam); % calculate the number of fixed crowd moniotring facilities
x_fiex=X_range_sam*rand(N_fixed_point,1); % generate fixed positions
y_fixed=Y_range_sam*rand(N_fixed_point,1);  

x_grid_fixed=ceil(x_fiex/delta_x); y_grid_fixed=ceil(y_fixed/delta_y); % specify which grid is the fixed point

t_c=band_arrival2(lambda,T,delta); % generate the times when scattered observations
% occur, through a poisson process

N_obs_c=length(t_c); % number of scattered crowd observations
N_obs_c_fixed=N_t*N_fixed_point; % number of crowd observations from fixed point
x_obs_c=X_range_sam*rand(N_obs_c,1); % x axis values of observations
y_obs_c=Y_range_sam*rand(N_obs_c,1); % y axis values of observations
% specify the grids of observations
t_obs_c_grid=ceil(t_c/delta); x_obs_c_grid=ceil(x_obs_c/delta_x); y_obs_c_grid=ceil(y_obs_c/delta_y);
% get the rain data, crowd_true = scatterd + fixed
crowd_true=zeros(N_obs_c+N_obs_c_fixed,1);
crowd_rep=zeros(N_obs_c+N_obs_c_fixed,1);
for i=1:N_obs_c % first N_obs_c points are scatted observation
   crowd_true(i)=rain_field(x_obs_c_grid(i),y_obs_c_grid(i),t_obs_c_grid(i));
   % generate the level of observation errors
   sd_level=sigma_obs;
   crowd_rep(i)=normrnd(crowd_true(i),sd_level*crowd_true(i))+bias*crowd_true(i);
end

for i=1:N_fixed_point
    crowd_true(N_obs_c+(i-1)*N_t+1:N_obs_c+i*N_t)=squeeze(rain_field(x_grid_fixed(i),y_grid_fixed(i),:));
       % generate the level of observation errors
   sd_level=sigma_obs;
   crowd_rep(N_obs_c+(i-1)*N_t+1:N_obs_c+i*N_t)=normrnd(crowd_true(N_obs_c+(i-1)*N_t+1:N_obs_c+i*N_t),...
       crowd_true(N_obs_c+(i-1)*N_t+1:N_obs_c+i*N_t).*repmat(sd_level,N_t,1))+...
       bias*crowd_true(N_obs_c+(i-1)*N_t+1:N_obs_c+i*N_t);
end

% calculate the reported position of scattered observations, by adding
% errors of GPS
x_rep_c_grid=ceil(normrnd(x_obs_c,repmat(sigma_x,N_obs_c,1))/delta_x);
y_rep_c_grid=ceil(normrnd(y_obs_c,repmat(sigma_y,N_obs_c,1))/delta_y);
for i=1:N_fixed_point % adding the reported position of fixed observations, which are assumed to have no location error
    x_rep_c_grid(N_obs_c+(i-1)*N_t+1:N_obs_c+i*N_t)=repmat(x_grid_fixed(i),N_t,1);
    y_rep_c_grid(N_obs_c+(i-1)*N_t+1:N_obs_c+i*N_t)=repmat(y_grid_fixed(i),N_t,1);
    t_obs_c_grid(N_obs_c+(i-1)*N_t+1:N_obs_c+i*N_t)=[1:N_t]';
end

Positive_ID=crowd_rep>0;
crowd_rep=crowd_rep(Positive_ID);
crowd_true=crowd_true(Positive_ID);
x_rep_c_grid=x_rep_c_grid(Positive_ID);
y_rep_c_grid=y_rep_c_grid(Positive_ID);
t_obs_c_grid=t_obs_c_grid(Positive_ID);

rep_bias=mean(mean((crowd_rep-crowd_true)./crowd_true));
% the actual measurement bias as a ratio of the average rainfall

% calculate the grid where the reported position is
x_ratio=x_step/delta_x; y_ratio=y_step/delta_y; t_ratio=t_step/delta; 

crowd_field=-2000*ones(N_x_grid/x_ratio,N_y_grid/y_ratio,N_t/t_ratio); 
% the variable to store crowd field obsverations
for i=1:N_t/t_ratio
    for j=1:N_x_grid/x_ratio
        for k=1:N_y_grid/y_ratio
            
            ind=(ismember(t_obs_c_grid,(i-1)*t_ratio+1:i*t_ratio))...
                &(ismember(x_rep_c_grid,(j-1)*x_ratio+1:j*x_ratio))...
                &(ismember(y_rep_c_grid,(k-1)*y_ratio+1:k*y_ratio));
% get the indices of observations occured in the specified position            
            if sum(ind)~=0
                crowd_field(j,k,i)=mean(crowd_rep(ind)); % otherwise set the value
% as the average of all observations in the specified position               
            else
            end
            
        end
    end
end

[x_no_obs,y_no_obs,t_no_obs]=ind2sub(size(crowd_field),find(crowd_field==-2000));
% identify the grid points where there is no observation, the x,y and t
% grid indices are stored in x_no_obs, y_no_obs and t_no_obs
[x_yes_obs,y_yes_obs,t_yes_obs]=ind2sub(size(crowd_field),find(crowd_field~=-2000));
% identify the grid points where there are observations,the x,y and t
% grid indices are stored in x_yes_obs, y_yes_obs and t_yes_obs
if ~isempty(x_no_obs)

    for i=1:N_t % for every timestep, set the values of the closet grid point as the value of grid point with no observation
        % the potentially sellectable grid points are the grid points
        % happened at timestep i
        
        index_NO=t_no_obs==i;
        x_NO_index=x_no_obs(index_NO);
        y_NO_index=y_no_obs(index_NO);
         % the grid points that have no observation at time i
        
        index_YES=t_yes_obs==i;
        x_YES_index=x_yes_obs(index_YES);
        y_YES_index=y_yes_obs(index_YES);
         % all the grid points that have observation at time i
        if isempty(x_YES_index)&&i==1
            crowd_field(:,:,i)=0;
        elseif isempty(x_YES_index)
            crowd_field(:,:,i)=crowd_field(:,:,i-1);
        else            
        for j=1:length(x_NO_index) % for all the grid points without observation at time i
            
            [~,ind_min]=min((x_YES_index-x_NO_index(j)).^2+(y_YES_index-y_NO_index(j)).^2);    
            % identify the closet grid point to the grid point that have no observation
            crowd_field(x_NO_index(j),y_NO_index(j),i)=crowd_field(x_YES_index(ind_min),y_YES_index(ind_min),i);
            % use the observation value reported by the grid point identified before as
            % the reported value
        end
        end
        
    end
    
else
end

end