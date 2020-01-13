function [k,TRCS,TC] = RWG_Sphere_v2(tesselation,iradius,ifreq,itheta,iphi,plotting,display,savepic)
%%
close all

% load custom colormap
% load 'magma.mat';

%EM parameters (f=3e8 means that f=300 MHz) 
f           = ifreq;  
epsilon_    = 8.85418782e-012;
mu_         = 1.25663706e-006;
%Speed of light
c_= 1/sqrt(epsilon_*mu_);
%Free-space impedance 
eta_= sqrt(mu_/epsilon_);
% Normal incidence
%Example: d=[0 0 -1] means that the incident signal
% is in the -z direction. 
theta   =itheta;
phi     =iphi;
d       =[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];
Pol     =[0 -cos(theta) sin(theta)*sin(phi)];
Pol     =Pol*sqrt(1/(norm(Pol)^2));
% radius of sphere
radius = iradius;

[pplot,tplot] = icosphere(2);
[po,t] = icosphere(tesselation);
p = po*radius;
pplot = pplot*radius;
t(:,4) = 1;
p = p';
t = t';

%% RWG1 Geometry calculations
Nf = length(t);
Nv = length(p);
Ne = Nf + Nv - 2;
TR = triangulation(t(1:3,:)',p');
Area = zeros(1,Nf);
Center = incenter(TR)';
%Find areas of separate triangles
for m=1:Nf
   N=t(1:3,m);
   Vec1=p(:,N(1))-p(:,N(2));
   Vec2=p(:,N(3))-p(:,N(2));
   Area(m) =norm(cross(Vec1,Vec2))/2;
end

%Find all edge elements "Edge_" with at least two 
%adjacent triangles
Edge_ = edges(TR)';
eA = edgeAttachments(TR,Edge_')';
eA = cell2mat(eA')';
TrianglePlus = eA(1,:);
TriangleMinus = eA(2,:);

%All structures of this chapter have EdgeIndicator=2
%EdgeIndicator=t(4,TrianglePlus)+t(4,TriangleMinus);

%Find edge length
EdgeLength = zeros(1,Ne);
for m=1:Ne
   EdgeLength(m)=norm(p(:,Edge_(1,m))-p(:,Edge_(2,m)));
end

%% RWG2 Geometry calculations
[RHO_Plus,RHO__Plus,RHO_Minus,RHO__Minus,Center_] = subtri(TrianglePlus,TriangleMinus,Ne,Nf,Center,Edge_,p,t);

%% RWG3 Calculates the impedance matrix using function IMPMET
%Contemporary variables - introduced to speed up 
%the impedance matrix calculation
omega       =2*pi*f;                                            
k           =omega/c_;
K           =1j*k;
Constant1   =mu_/(4*pi);
Constant2   =1/(1j*4*pi*omega*epsilon_);
Factor      =1/9;    

FactorA     =Factor*(1j*omega*EdgeLength/4)*Constant1;
FactorFi    =Factor*EdgeLength*Constant2;

RHO_P = zeros(3,9,Ne);
RHO_M = zeros(3,9,Ne);
for m=1:Ne
    RHO_P(:,:,m)=repmat(RHO_Plus(:,m),[1 9]);   %[3 9 EdgesTotal]
    RHO_M(:,:,m)=repmat(RHO_Minus(:,m),[1 9]);  %[3 9 EdgesTotal]
end
FactorA=FactorA.';
FactorFi=FactorFi.';


Z=  impmet( Ne,Nf,...
            EdgeLength,K,...
            Center,Center_,...
            TrianglePlus,TriangleMinus,...
            RHO_P,RHO_M,...
            RHO__Plus,RHO__Minus,...
            FactorA,FactorFi);   
%PLOT PLOT PLOT PLOT PLOT
% figure(1);
% imagesc(abs(Z)); colorbar(); pbaspect([1 1 1]);

%% RWG4 Solves MoM equations for the scattering problem

kv=k*d;

V = zeros(1,Ne);
for m=1:Ne    
   ScalarProduct=sum(kv.*Center(:,TrianglePlus(m))');
   EmPlus =Pol.'*exp(-1j*ScalarProduct);      
   ScalarProduct=sum(kv.*Center(:,TriangleMinus(m))');
   EmMinus=Pol.'*exp(-1j*ScalarProduct);      
   ScalarPlus =sum(EmPlus.* RHO_Plus(:,m));
   ScalarMinus=sum(EmMinus.*RHO_Minus(:,m));
   V(m)=EdgeLength(m)*(ScalarPlus/2+ScalarMinus/2);   
end
I=Z\V.';

%% RWG5 Visualizes the surface current magnitude

%Find the current density for every triangle
CurrentNorm = zeros(1,Nf);
for kk=1:Nf
    i=[0 0 0]';
    for m=1:Ne
        IE=I(m)*EdgeLength(m);
        if(TrianglePlus(m)==kk)
            i=i+IE*RHO_Plus(:,m)/(2*Area(TrianglePlus(m)));
        end
        if(TriangleMinus(m)==kk)
            i=i+IE*RHO_Minus(:,m)/(2*Area(TriangleMinus(m)));
        end
    end
    CurrentNorm(kk)=abs(norm(i));
end

%Jmax=max(CurrentNorm);
%MaxCurrent = strcat(num2str(Jmax),'[A/m]');
%CurrentNorm1=CurrentNorm/max(CurrentNorm);


%% RCS 
DipoleCenter = zeros(3,Ne);
DipoleMoment = zeros(3,Ne);

p_observation = po';
for m=1:Ne
    Point1=Center(:,TrianglePlus(m));
    Point2=Center(:,TriangleMinus(m));
    DipoleCenter(:,m)=0.5*(Point1+Point2);
    DipoleMoment(:,m)=EdgeLength(m)*I(m)*(-Point1+Point2); 
end

TotalPower=0;
Poynting = zeros(3,Nf);
Area__ = zeros(1,Nf);
U = zeros(1,Nf);
for m=1:Nf
    N=t(1:3,m);
    ObservationPoint=1/3*sum(p_observation(:,N),2);
    [E,H]=point(ObservationPoint,eta_,K,DipoleMoment,DipoleCenter);
    ET=sum(E,2); HT=sum(H,2);
    Poynting(:,m)=0.5*real(cross(ET,conj(HT)));
    U(m)=(norm(ObservationPoint))^2*norm(Poynting(:,m)); 
    Vector1=p_observation(:,N(1))-p_observation(:,N(2));
    Vector2=p_observation(:,N(3))-p_observation(:,N(2));
    Area__(m) =0.5*norm(cross(Vector1,Vector2)); 
    TotalPower=TotalPower+norm(Poynting(:,m))*Area__(m);
end

%GainLogarithmic     =10*log10(4*pi*max(U)/TotalPower);
%GainLinear          =4*pi*max(U)/TotalPower;
%UNorm=U/norm(U);

TRCS = TotalPower*2*eta_;
TC = sum(CurrentNorm.*Area);




if plotting
    close all;
    if ~display
        figure('visible','off');
    else
        figure(1);
    end
    fig1=patch('Faces',t(1:3,:)','Vertices',p',...
        'VertexNormals',p',...
        'LineStyle','none',...%'LineWidth',0.1,&...
        'FaceLighting','phong',...
        'BackFaceLighting','unlit',...
        'AmbientStrength',0.3,'DiffuseStrength',0.6,...
        'SpecularExponent',10,'SpecularStrength',0.9,...
        'FaceColor','flat','CData',CurrentNorm,...
        'Tag','Icosphere');
    cb2 = colorbar();
    title('Magnitude of current along surface')
    ylabel(cb2,'Absolute Current [A/m]');
    pbaspect([1 1 1]);
    xlabel('X[m]');
    ylabel('Y[m]');
    zlabel('Z[m]');
    rotate3d
    view(45,45) 
%     colormap(magma)
    if savepic
        saveas(gcf,['J',num2str(radius/1e-2),'cm.png'])
    end
    
    if ~display
        figure('visible','off');
    else
        figure(2);
    end
    fig2=patch('Faces',tplot','Vertices',pplot);
    hold on;
    fig3 = patch('Faces',t(1:3,:)','Vertices',p_observation',...
        'VertexNormals',p_observation',...
        'LineStyle','none',...%'LineWidth',0.1,...
        'FaceLighting','phong',...
        'BackFaceLighting','unlit',...
        'AmbientStrength',0.3,'DiffuseStrength',0.6,...
        'SpecularExponent',10,'SpecularStrength',0.9,...
        'FaceColor','flat','CData',U,...
        'Tag','Icosphere');
    set(fig3,'FaceAlpha',0.75);
    cb3 = colorbar();
    title(['R=',num2str(radius/1e-2),'cm, Scattering observed from a r=1m sphere'])
    ylabel(cb3,'RCS');
    pbaspect([1 1 1]);
    xlabel('X[m]');
    ylabel('Y[m]');
    zlabel('Z[m]');
    rotate3d
    view(45,45) 
%     colormap(magma)
    if savepic
        saveas(gcf,['RCS',num2str(radius/1e-2),'cm.png'])
    end
end

end