function rain_field=generate_rain_WRG_large_domain_slow(lambda_m,T,delta,rou_L,X_Lim,Y_Lim,delta_x,delta_y,E_v,alpha,sigma,beta,Ub,i0,D)
sample_t(1,1,:)=delta:delta:T; % delta is the time step, T is the duration
N_t=length(sample_t); % N_t is the number of sample timeperiods
% sample_t are the sampling time periods
sample_point_X(:,1,1)=delta_x/2:delta_x:X_Lim-delta_x/2; % X_Lim is the X horizon range, x_step is the x stepsize
sample_point_Y(1,:,1)=delta_y/2:delta_y:Y_Lim-delta_y/2; % Y_Lim is the Y horizon range, y_step is the y stepsize
N_X=length(sample_point_X); N_Y=length(sample_point_Y);
% N_X and N_Y are the number of sample points in X and Y horizon
% sample_point_X and sample_point_Y are the sampling points
sample_t=repmat(sample_t,[N_X,N_Y,1]);
sample_point_X=repmat(sample_point_X,[1,N_Y,N_t]);
sample_point_Y=repmat(sample_point_Y,[N_X,1,N_t]);
% sample_t, sample_point_X and sample_point_Y are the repeated to 3D for
% matrix calculation

ts=band_arrival(lambda_m,T,delta,1/alpha); % lambda_m is the lamdbaM statisc
% ts are the time periods that rainbands arrives
N_band=length(ts); % N_band is the number of rainbands arrived during T
rain_i_band=zeros(N_X,N_Y,N_t,N_band); % all the rain intensities in each rainband

for i=1:N_band % sum up all the rain intensities in all rain bands
    CPC=cluster_potential_larger_border(rou_L,X_Lim,Y_Lim,Ub,alpha,D,sigma); % rou_L is the static rouL
    % CPC are the cluster potential centers in a rainband, first and second rows are the X and Y positions 
    [N_CPC,~]=size(CPC); % N_CPC is the number of CPCs in a rainband
    rain_i_CPC=zeros(N_X,N_Y,N_t);
    j=1;
    
    while j<=N_CPC % sum up all the rain intensities in all CPCs
        N_cell=poissrnd(E_v); % E_v is the statistic Ev, N_cell is the number of cells in a CPC
        rain_i_cell=zeros(N_X,N_Y,N_t,N_cell); % all the rain intensities in each cell
        [tao,cell_center]=cell_born(ts(i),1/alpha,CPC(j,:),sigma,beta,Ub,N_cell);
        
        for k=1:N_cell % sum up all the rain intensities in all cells
            rain_i_cell(:,:,:,k)=i0*(sample_t>tao(k)).*... % =0 if sample_t<tao (s-T in Valdes 1985)
                exp(-(sample_t-tao(k))*alpha).*... % g1(a) in Valdes 1985
                exp(-((sample_point_X-(cell_center(k,1)+Ub(1)*(sample_t-tao(k)))).^2+...
                (sample_point_Y-(cell_center(k,2)+Ub(2)*(sample_t-tao(k)))).^2)/(2*D^2));
            % g2(r) in Valdes 1985
            
        end
        
        rain_i_CPC=rain_i_CPC+sum(rain_i_cell,4);
        j=j+1;
    end
    
    rain_i_band(:,:,:,i)=sum(rain_i_CPC,4);
end

rain_field=sum(rain_i_band,4);
end