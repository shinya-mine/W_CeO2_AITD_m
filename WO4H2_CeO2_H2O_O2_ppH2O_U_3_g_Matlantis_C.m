%% clean working space
clear all
close all
format long  

%% control

name_head = 'WO4H2_CeO2_H2O_O2_ppH2O_g_Matlantis';


% aWO4H2 + CeO2 + bH2O + cO2 → (WxOyHz)/CeO2;

% Chemical potential: O2(g) & H2O(g)

% coefficient
% a=x, b=-x+(1/2)z, c=(-3/2)x+(1/2)y-(1/4)z


%%
% WxOyHz/CeO2
energy =[
    -1052.75690642081 % W1O4H2 (x=1, y=4, z=2)
    -1041.86750219616 % W1O3H0 (x=1, y=3, z=0)
    -1063.50950571173 % W1O5H4 (x=1, y=5, z=4)
    -1077.66199368742 % W2O7H2 (x=2, y=7, z=2)
    -1088.61691763514 % W2O8H4 (x=2, y=8, z=4)
    -1099.49582012334 % W2O9H6 (x=2, y=9, z=6)
    -1102.54019725891 % W3O10H2 (x=3, y=10, z=2)
    -1113.82203562399 % W3O11H4 (x=3, y=11, z=4)
    -1125.08005754843 % W3O12H6 (x=3, y=12, z=6)
    -1127.13776118487 % W4O13H2 (x=4, y=13, z=2)
    -1138.26184916674 % W4O14H4 (x=4, y=14, z=4)
    -1148.98382020161 % W4O15H6 (x=4, y=15, z=6)
    -1159.56310169897 % W4O16H8 (x=4, y=16, z=8)
    -1169.85726554221 % W4O17H10 (x=4, y=17, z=10)
    ];

%%
coefficient = [
    1   0   0
    1   -1   0
    1   1   0
    2	-1	0
    2	0	0
    2	1	0
    3	-2	0
    3	-1	0
    3	0	0
    4	-3	0
    4	-2	0
    4	-1	0
    4	0	0
    4	1	0
];

p_O2 = 0.1; %atm for fixed partial pressure
p_H2O = 0.00; %atm for fixed partial pressure

%%
energy_O2 = -6.035181266;
energy_H2O =  -10.03120878;   
energy_WO3 = -38.60191691; % W2O6: -77.20383382
energy_CeO2 = -1017.20013926025; % Ce48O96 slab
energy_WO6H6 = -78.83463101; % aW(OH)6 + CeO2 +bO2 + cH2O → (WxOyHz)/CeO2  
                             % coefficient: a=x, b=(-3/2)x+(1/2)y-(1/4)z, c=(1/2)z-3x
energy_WO4H2 = -35.556305737631; % aWO2H4 + CeO2 + bO2 + cH2O → (WxOyHz)/CeO2
                                 % coefficient: a=x, b=(-3/2)x+(1/2)y-(1/4)z, c=-x+(1/2)z

%%

%% O2 chemical potential

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluate the chemical potential of the O2
% http://janaf.nist.gov/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Oxygen (O2)	O2(ref)
%T(K)	S	H-H(Tr)
table_O2 = [
0	0	-8.683
100	173.307	-5.779
200	193.485	-2.868
250	199.99	-1.41
300	205.329	0.054
350	209.88	1.531
400	213.871	3.025
450	217.445	4.543
500	220.693	6.084
600	226.451	9.244
700	231.466	12.499
800	235.921	15.835
900	239.931	19.241
1000	243.578	22.703
1100	246.922	26.212
1200	250.01	29.761
1300	252.878	33.344
];

for Temp = 1:1:14

p0 = 0.986923169314 ;% in atm 

T_O2 = table_O2(Temp,1); 
R  = 8.31451 * 0.00001036427; % eV/K
H0_O2 = table_O2(Temp,3);
S0_O2 = table_O2(Temp,2); 

%a_O2 = [-8:1:5]; % look at the corresponding mu interval!*
%p_O2 = 10.^(a_O2);

delta_mu_O2_p0 = 0.0103642 * ( H0_O2 - table_O2(1,3) - (T_O2 .* S0_O2) ./ 1000) ;
delta_mu_O2 = delta_mu_O2_p0 + (R .* T_O2) * log ( p_O2 ./ p0);

cp_O2(:,Temp) = delta_mu_O2;


end


for ii = 1:1:14
    cp_O2(ii,:) = cp_O2(1,:);
end

%%

%% H2o chemical potential

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluate the chemical potential of the H2O
% http://kinetics.nist.gov/janaf/html/H-064.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Water (H2O)           H2O1(g)

% Enthalpy Reference Temperature = Tr = 298.15 K  
% Standard State Pressure = p° = 0.1 MPa = 0.986923169314 atm
% T              S°               H-H°(Tr)    
% K             J·mol-1·deg-1     kJ·mol-1        
table_H2O = [
0	0	-9.904
100	152.388	-6.615
200	175.485	-3.282
300	189.042	0.062
400	198.788	3.452
500	206.534	6.925
600	213.052	10.501
700	218.739	14.192
800	223.825	18.002
900	228.459	21.938
1000	232.738	26
1100	236.731	30.191
1200	240.485	34.506
1300	244.035	38.942
];

for Temp = 1:1:14

p0 = 0.986923169314 ;% in atm 



T_H2O = table_H2O(Temp,1); 
R  = 8.31451 * 0.00001036427; % eV/K
H0_H2O = table_H2O(Temp,3); 
S0_H2O = table_H2O(Temp,2); 

a_H2O = [-8:1:5]; % look at the corresponding mu interval!*
p_H2O = 10.^(a_H2O);

delta_mu_H2O_p0 = 0.0103642 * ( H0_H2O - table_H2O(1,3) - (T_H2O .* S0_H2O) ./ 1000) ; % delta_mu_H2O_p0 for T_H2O = 300 K

delta_mu_H2O = delta_mu_H2O_p0 + (R .* T_H2O) * log ( p_H2O ./ p0);


cp_H2O(:,Temp) = delta_mu_H2O;

end


pp_H2O = zeros(11,14);
pp_H2O = a_H2O.';

%check matlix

for iii = 1:1:14
    temperature(iii,:) = [0:100:1300];
    pp_H2O(:,iii) = pp_H2O(:,1);
end

%%


for i = 1:1:14
    E_O2(i,:) = energy_O2;
    E_H2O(i,:) = energy_H2O;
    E_WO3(i,:) = energy_WO3;
    E_CeO2(i,:) = energy_CeO2;
    E_WO6H6(i,:) = energy_WO6H6;
    E_WO4H2(i,:) = energy_WO4H2;
    
    E_O2(:,i) = energy_O2;
    E_H2O(:,i) = energy_H2O;
    E_WO3(:,i) = energy_WO3;
    E_CeO2(:,i) = energy_CeO2;
    E_WO6H6(:,i) = energy_WO6H6;
    E_WO4H2(:,i) = energy_WO4H2;
    
end

for j = 1:1:size(energy)
    for i = 1:1:14
        
       % if i < 12 
        E(i,:,j) = energy(j,1);
        E(:,i,j) = energy(j,1);

    end
    
    Z(:,:,j) = E(:,:,j) - E_CeO2 - coefficient(j,1)*E_WO4H2 - coefficient(j,2)*E_H2O - coefficient(j,3)*E_O2;
    %G(:,:,j) = Z(:,:,j) + (-coefficient(j,2)*cp_H2O-coefficient(j,3)*cp_O2);
    G(:,:,j) = Z(:,:,j) + (-coefficient(j,2)*cp_H2O); 
    
end

%% PLOT

colors = [
    128./255 0./255 128./255   ; ... % W1O4H2 (x=1, y=4, z=2)       
    153./255 50./255 204./255   ;...   % W1O3H0 (x=1, y=3, z=0)  
    186./255 85./255 211./255  ; ... % W1O5H4 (x=1, y=5, z=4)      
    65./255 105./255 225./255   ; ... % W2O7H2 (x=2, y=7, z=2)   
    135./255 206./255 250./255 ; ... % W2O8H4 (x=2, y=8, z=4)                  
    0./255 206./255 209./255 ; ... % W2O9H6 (x=2, y=9, z=6)
    50./255 205./255 50./255   ; ... % W3O10H2 (x=3, y=10, z=2) 
    0./255 128./255 0./255 ; ... % W3O11H4 (x=3, y=11, z=4)
    102./255 205./255 170./255 ; ... % W3O12H6 (x=3, y=12, z=6)***
    250./255 128./255 114./255 ; ... % W4O13H2 (x=4, y=13, z=2)***
    255./255 105./255 180./255 ; ... % W4O14H4 (x=4, y=14, z=4)
    255./255 192./255 203./255 ; ... % W4O15H6 (x=4, y=15, z=6)
    255./255 165./255 0./255 ; ... % W4O16H8 (x=4, y=16, z=8)
    255./255 69./255 0./255 ; ... % W4O17H10 (x=4, y=17, z=10)
           ];

figure(1);

hold on

for i=1:size(energy)
    hSurface = surf(temperature-273,pp_H2O,G(:,:,i),'EdgeColor','none');
    set(hSurface,'FaceColor',colors(i,:));
    xlabel('Temperature (℃)','FontSize',15,'FontName','Arial')
    ylabel('log(p_{H_{2}O})','FontSize',15,'FontName','Arial')
    zlabel('\DeltaG (eV/{A}^2)','FontSize',15,'FontName','Arial')
    title('')
end

legend('W_{1}O_{4}H_{2}','W_{1}O_{3}','W_{1}O_{5}H_{4}',...
    'W_{2}O_{7}H_{2}','W_{2}O_{8}H_{4}','W_{2}O_{9}H_{6}',...
    'W_{3}O_{10}H_{2}','W_{3}O_{11}H_{4}','W_{3}O_{12}H_{6}',...
    'W_{4}O_{13}H_{2}','W_{4}O_{14}H_{4}','W_{4}O_{15}H_{6}',...
    'W_{4}O_{16}H_{8}','W_{4}O_{17}H_{10}',...
    'location','northeastoutside','FontSize',12,'FontName','Arial')
%legend(comp,'location','eastoutside','FontSize',legand_fontsize,'FontName','Arial')

% turn the grid on
grid on;

% set the view
view(3);
ax1 = gca;
ax1.XTick = 0:100:900;
set(ax1, 'XLim',[0 800],'TickDir','out', 'box','on',...
    'YLim',[-8 5],'ZLim',[-8 2],'FontName','Arial','Fontsize',15);
pbaspect([5 5 3])

 print_name_3d = strcat(name_head,'-3d.jpg');
 print(print_name_3d,'-djpeg','-r600')


%%

figure(2)
scrsz = get(0,'ScreenSize');
hf2=figure(2);
set(hf2,'Position',[0.8,1.2,0.8,1.2].*scrsz); %full screensize
set(gcf,'Renderer','ZBuffer','RendererMode','manual');
hold on;

% turn the grid on
%grid on;

% set the view
view([0 0 -10]);

% draw the planes in their respective colors
for i=1:size(energy)
    hSurface = surf(temperature-273,pp_H2O,G(:,:,i),'EdgeColor','none');
    set(hSurface,'FaceColor',colors(i,:));
end

% title
%title({'\fontsize{12} \bf FAU (T = 600K)\rm ';'';''});
title({' ';'';''});

% set text for the default x/y axes
ax1 = gca;
set(ax1,'XColor','k','YColor','k','TickDir','out','box','on',...
    'Position',[0.2,0.2,0.6,0.6],'FontName','Arial','Fontsize',30);
ax1.XTick = 0:100:900;
xlabel('Temperature (℃)','Fontsize',30,'FontName','Arial');
ylabel('log(p_{H_{2}O})','Fontsize',30,'FontName','Arial');
set(ax1,'XLim',[0 800],...
        'YLim',[-6 1],...
        'LineWidth',1.5,...
        'XAxisLocation','bottom',...
        'Layer','Top',...
        'Box','off');
pbaspect([1.25 1 1])
% 
% % create another pair of axes
% ax2 = axes('Position',get(ax1,'Position'),...
%            'XAxisLocation','top',...
%            'YAxisLocation','right',...
%            'Color','none',...
%            'XColor','k','YColor','k','YDir','reverse');
%        
%        
% x2 = [a_O2];
% y2 = [a_H2O];
% set(ax2,'FontSize',20);
% xlabel('Temperature (K)','Fontsize',30,'FontName','Arial'); 
% ylabel('log(p_{H_2})','Fontsize',30,'FontName','Arial');
% x2s=size(x2,2);
% y2s=size(y2,2);
% set(ax2,'XLim',[0 1300],'YLim',[-8 5],'Layer','Bottom','FontName','Arial','Fontsize',30);
% 


 print_name_2d = strcat(name_head,'-2d.jpg');
 print(print_name_2d,'-djpeg','-r600')
 