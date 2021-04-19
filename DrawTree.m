function values = DrawTree(chromosome,x0,y0,ss,variableRangeMatrix)
%%
% n,rR, RL,rSigma,phiR, phiL,phiSigma,yb,colorTreeR,colorTreeG,colorTreeB
[n nNorm] = DecodeChromosome(chromosome(1,:),variableRangeMatrix(1,:));
n = round(n);
[rR rRnorm] = DecodeChromosome(chromosome(2,:),variableRangeMatrix(2,:));
[rL rLnorm] = DecodeChromosome(chromosome(3,:),variableRangeMatrix(3,:));
r = [rL rR];
[rSigma rSigmaNorm] = DecodeChromosome(chromosome(4,:),variableRangeMatrix(4,:));
[phiR phiRnorm] = DecodeChromosome(chromosome(5,:),variableRangeMatrix(5,:));
[phiL phiLnorm]= DecodeChromosome(chromosome(6,:),variableRangeMatrix(6,:));
phi = [phiL phiR];
[phiSigma phiSigmaNorm]= DecodeChromosome(chromosome(7,:),variableRangeMatrix(7,:));
[x1 x1norm] = DecodeChromosome(chromosome(8,:),variableRangeMatrix(8,:));
[y1 y1norm] = DecodeChromosome(chromosome(9,:),variableRangeMatrix(9,:));
colorTreeR = DecodeChromosome(chromosome(10,:),variableRangeMatrix(10,:));
colorTreeG = DecodeChromosome(chromosome(11,:),variableRangeMatrix(11,:));
colorTreeB = DecodeChromosome(chromosome(12,:),variableRangeMatrix(12,:));
xb = [x0 x1]; yb = [y0 y1];
values = [nNorm rRnorm rLnorm rSigmaNorm phiRnorm phiLnorm phiSigmaNorm x1norm y1norm colorTreeR colorTreeG colorTreeB];
% This function allow you to bild  fractal trees
%  by using modified algorithms based on the so-called Kantor`s array
%  and method of inverse trace
%  These methods allow you to economise time and computer memory
%  considerably
% It`s arguments:
%  Fraktal(n,r,phi,xb,yb)
%     n - number of iterations
%     r - scale factor
%       in case of nonuniform fractals where different scale-factors take
%       place, r may be vector (in this case the lengths of r and phi must be equal)
%     phi - vector of angles in fractal generator that are calculated
%     relative to the horizontal axis connecting end points of generator
%     xb and yb - coordinates of trunk
  
% axes(ss);

mN=length(phi);
if length(r)==1
    rM=ones(1,mN)*r;
elseif length(r)==mN
    rM=r;
else
    warndlg('The sizes of scale vector and vector of angles in fraktal`s generator don`t equal');
end

% The adress to the system to ignore warnings as divide by zero
warning off MATLAB:divideByZero

% ksi - angle of trunk with vertical axes
% it`s calculated on the basis of coordinates xb and yb of trunk
%  by using the line-equation:
%  y=k*x+d,   here
%      k=(yb(2)-yb(1))/(xb(2)-xb(1));
%      d=(yb(1)*xb(2)-xb(1)*yb(2))/(xb(2)-xb(1));
a=sqrt((xb(1)-xb(2))^2+(yb(1)-yb(2))^2);   % scale factor
b=1;
c=sqrt((1-xb(2))^2+(yb(2))^2);
alpha=acos((a^2+b^2-c^2)/(2*a*b));
k=(yb(2)-yb(1))/(xb(2)-xb(1));                       % define the tangent
d=(yb(1)*xb(2)-xb(1)*yb(2))/(xb(2)-xb(1));           % define the free member
if yb(2)>=yb(1)
    ksi=alpha;
else
    ksi=-alpha;
end

%    auxiliary vectors for theta and rD (see below their definition)
psi=zeros(n,mN^n);
ralt=ones(n,mN^n);
%   psi and ralt have Kantor array`s structur !!!

for i=1:1:n
    z=1;
    for j=1:1:mN^(i-1)
        for k=1:1:mN
            for m=1:1:mN^(n-i)
%                 psi(i,z)=phi(k);               % define psi on the base of phi
                psi(i,z)= normrnd(phi(k),phiSigma);
%                 ralt(i,z)=rM(k);               % define ralt on the base of rM
                ralt(i,z)= normrnd(rM(k),rSigma);               % define ralt on the base of rM
                z=z+1;
            end
        end
    end
end

theta=zeros(n,mN^n);  % vector of angles between each branch of fractal and vertical axes
rD=ones(n,mN^n);      % lengths of branches

for i=1:1:mN^n
    for j=1:1:n
        for k=1:1:j
            theta(j,i)=theta(j,i)+psi(k,i);        % define theta on the base of psi
            rD(j,i)=rD(j,i)*ralt(k,i);             % define rD on the base of rD
        end
    end
end

theta=theta+ksi;


% ----Matrix for coordinates-------
A=ones(n+1,mN^n,2);
% initial coordinates
A(1,:,1)=ones(1,mN^n)*xb(2);   % x-coordinate
A(1,:,2)=ones(1,mN^n)*yb(2);   % y-coordinate
for j=1:1:mN^n
    for i=1:1:n
        
        % define following coordinates
        A(i+1,j,1)=A(i,j,1)+a*rD(i,j)*cos(theta(i,j));
        A(i+1,j,2)=A(i,j,2)+a*rD(i,j)*sin(theta(i,j));
    end
end

% By visualisation the method of inverse trace is used
tau=1;
for i=1:mN:mN^n
    z=1;
    for k=1:1:mN
        for j=z:1:n
            line([A(j,i,1),A(j+1,i,1)],[A(j,i,2),A(j+1,i,2)], ...
'Color',[colorTreeR colorTreeG colorTreeB 0.75],'LineWidth',2-2*(j-1)./(1+abs(j-0.5)));
        end
        z=z+1;
    end
end

% -------------- Trunk -----------------
line([xb(1),xb(2)],[yb(1),yb(2)],'Color',[colorTreeR colorTreeG colorTreeB],'LineWidth',2);
set(gca,'xticklabels','','yticklabels','')
