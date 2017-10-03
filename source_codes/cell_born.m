function [tao,x]=cell_born(ts,T,CPC,sigma,beta,Ub,N_cell)
tao=exprnd(1/beta,N_cell,1)+ts-T; % time of arrival of the band - T
% fs(1)t in Valdes 1985
x=mvnrnd([0,0],[sigma(1),0;0,sigma(2)],N_cell)+...
    repmat(CPC,[N_cell,1])+...
    repmat(Ub,[N_cell,1]).*repmat((tao-ts),[1,2]);
% f(2)x in Valdes 1985
end