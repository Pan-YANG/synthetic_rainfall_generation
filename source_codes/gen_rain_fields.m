function [results]=gen_rain_fields(rainfield_param,study_area_pos)
%% get the parameters for generating rainfield
lambda_m=rainfield_param.lambda_m; T=rainfield_param.T; 
delta=rainfield_param.delta; rou_L=rainfield_param.rou_L;
X_Lim=rainfield_param.X_Lim; Y_Lim=rainfield_param.Y_Lim; 
delta_x=rainfield_param.delta_x; delta_y=rainfield_param.delta_y;
E_v=rainfield_param.E_v; alpha=rainfield_param.alpha; 
sigma=rainfield_param.sigma; beta=rainfield_param.beta; 
Ub=rainfield_param.Ub; i0=rainfield_param.i0; D=rainfield_param.D; 
% get the parameters for generating gauge field
N_gauge=rainfield_param.N_gauge; smooth_ratio=rainfield_param.smooth_ratio; 
t_step=rainfield_param.t_step; x_step=rainfield_param.x_step; 
y_step=rainfield_param.y_step;
% get the parameters for generating crowd field
lambda_c=rainfield_param.lambda_c; pop_den=rainfield_param.pop_den; 
X_range_sam=rainfield_param.X_range_sam; Y_range_sam=rainfield_param.Y_range_sam; 
sigma_x=rainfield_param.sigma_x; sigma_y=rainfield_param.sigma_y;
sigma_obs=rainfield_param.sigma_obs; 
fixed_point_den=rainfield_param.fixed_point_den;
bias=rainfield_param.bias;

%% ground truth rainfall field generation
% First start with generating the rainfield
    % generate the rainfield on a larger space than the study area
   rain_field_large=generate_rain_WRG_large_domain_slow...
       (lambda_m,T,delta,rou_L,X_Lim,Y_Lim,delta_x,delta_y,E_v,alpha,sigma,beta,Ub,i0,D); 
   % then get the rain field for study area
   rain_field=rain_field_large...
       (ceil(study_area_pos(1,1)/delta_x)+1:ceil(study_area_pos(2,1)/delta_x)...
       ,ceil(study_area_pos(1,2)/delta_y)+1:ceil(study_area_pos(2,2)/delta_y),:); 
   
   %% rain gauge rainfall field generation
    % generate rain gauge observations    
        [raingauge_obs,raingauge_pos]=get_raingauge_obs(rain_field_large,...
        X_Lim,Y_Lim,delta_x,delta_y,N_gauge,smooth_ratio);
    % generate rain field using theissen
    gauge_field_large=generate_gauge_field...
        (raingauge_obs,raingauge_pos,...
        X_Lim,Y_Lim,x_step,y_step,t_step/delta);
    gauge_field=gauge_field_large...
        (ceil(study_area_pos(1,1)/x_step)+1:ceil(study_area_pos(2,1)/x_step),...
        ceil(study_area_pos(1,2)/y_step)+1:ceil(study_area_pos(2,2)/y_step),:);
    % generate rain field using linear interpolation
    gauge_field_large=generate_gauge_field_inter...
        (raingauge_obs,raingauge_pos,...
        X_Lim,Y_Lim,x_step,y_step,t_step/delta,'linear');
    gauge_field_linear=gauge_field_large...
        (ceil(study_area_pos(1,1)/x_step)+1:ceil(study_area_pos(2,1)/x_step),...
        ceil(study_area_pos(1,2)/y_step)+1:ceil(study_area_pos(2,2)/y_step),:);
    % generate rain field using natural interpolation
    gauge_field_large=generate_gauge_field_inter...
            (raingauge_obs,raingauge_pos,...
        X_Lim,Y_Lim,x_step,y_step,t_step/delta,'natural');
    gauge_field_natural=gauge_field_large...
        (ceil(study_area_pos(1,1)/x_step)+1:ceil(study_area_pos(2,1)/x_step),...
        ceil(study_area_pos(1,2)/y_step)+1:ceil(study_area_pos(2,2)/y_step),:);
    
    %% crowd field generation
    [crowd_field,~]=generate_crowd_field_withfixedpoint_2dinter(rain_field,bias,...
        lambda_c,pop_den,fixed_point_den,X_range_sam,Y_range_sam,...
        delta_x,delta_y,sigma_x,sigma_y,T,delta,sigma_obs,...
        t_step,x_step,y_step);
    
    %% outputing results
    results.rain_field=rain_field;
    results.gauge_field=gauge_field;
    results.gauge_field_linear=gauge_field_linear;
    results.gauge_field_natural=gauge_field_natural;
    results.crowd_field=crowd_field;
    