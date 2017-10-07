function [R, t, TData] = BdICP(R0, t0, MoveStep)
% input:
%          R0: initial rotation matrix;
%          t0: initial translation vector;
%          MoveStep: maximum iteration 
% output  
%          R: rotation matrix
%          t: translation vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref:
% Jihua Zhu, Di Wang, Xiuxiu Bai, Huimin Lu, Congcong Jin, Zhongyu Li.
% Registration of Point Clouds Based on the Ratio of Bidirectional Distances. 
% The International Conference on 3D Vision(3DV), 2016: 102-107.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jihua Zhu, mailzhujh@xjtu.edu.cn


global Model Data
CurrStep = 0;
Mns= createns(Model');                                   % build k-d tree for NN search
Dns= createns(Data');
phiMd= 100;                                           
minMd= 10;
R= R0;
t= t0;
while ((CurrStep < MoveStep)&(abs(minMd-phiMd)>10^(-6)))
    TData = transform_to_global(Data, R, t);            
    phiMd= minMd;
    [corr1,TD1] = knnsearch(Mns,TData');   
    minMd= mean(TD1);
    [C,IA,IC]= unique(corr1);
    TModel = transform_to_global(Model(:,C), inv(R),-inv(R)*t);
    [corr2,TD2] = knnsearch(Dns,TModel');
    TD21= TD2(IC);
    k= (TD1+10^-6)./(TD21+10^-6);
  % p= exp(-(k.^4-1));                                     % the  probability can be calculated by different ways             
    p= exp(-6*(k-1));        
    corr1(:,2) = [1 : length(corr1)]';                   
    [R, t] = reg(corr1, p);
    CurrStep= CurrStep+1;
end

% Solution for the weighted least square problem
function [R1, t1] = reg(corr,p)
global Model Data
M = Model(:,corr(:,1)); 
S = Data(:,corr(:,2));
sumP= sum(p);
p1= [p';p';p'];
mm= sum((M.*p1)')./sumP;
ms= sum((S.*p1)')./sumP;
Sshifted = [S(1,:)-ms(1); S(2,:)-ms(2); S(3,:)-ms(3)];
Mshifted = [M(1,:)-mm(1); M(2,:)-mm(2); M(3,:)-mm(3)];
K= Sshifted.*p1*Mshifted';
[U A V] = svd(K);
R1 = V*U';
if det(R1)<0
    B = eye(3);
    B(3,3) = det(V*U');
    R1 = V*B*U';
end
t1= sum(((M-R1*S).*p1)')./sumP;
t1= t1';