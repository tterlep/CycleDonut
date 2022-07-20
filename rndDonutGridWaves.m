clear


%declare weights and the mode of operation
a1 =(2-sqrt(2))/2;
a2 = (3-sqrt(5))/2;
b1 = sqrt(2) - 1;
b2 = sqrt(5) - 2;

% a1 = 0.1;
% a2 = 0.2;
% b1 = 0.5;
% b2 = 0.2;


mode = 1;


% a1 = 4/10;
% a2 = 4/10;
% b1 = 4/10;
% b2 = 0.8;

Ag = [...
    0 a1 b1 a1 a2 a2 0 0;...
    a1 0 a1 b1 0 a2 a2 0;...
    b1 a1 0 a1 0 0 a2 a2;...
    a1 b1 a1 0 a2 0 0 a2;...
    a2 0 0 a2 0 b2 0 b2;...
    a2 a2 0 0 b2 0 b2 0;...
    0 a2 a2 0 0 b2 0 b2;...
    0 0 a2 a2 b2 0 b2 0];

% F = Ag;
% F(F>0) = 1;
%
% a1 = 1*a1;
% a2 = 1*a2;
% b1 = 1*b1;
% b2 = 1*b2;
%
% Ag = [...
%     0 a1*rand b1*rand a1*rand a2*rand a2*rand 0 0;...
%     a1*rand 0 a1*rand b1*rand 0 a2*rand a2*rand 0;...
%     b1*rand a1*rand 0 a1*rand 0 0 a2*rand a2*rand;...
%     a1*rand b1*rand a1*rand 0 a2*rand 0 0 a2*rand;...
%     a2*rand 0 0 a2*rand 0 b2*rand 0 b2*rand;...
%     a2*rand a2*rand 0 0 b2*rand 0 b2*rand 0;...
%     0 a2*rand a2*rand 0 0 b2*rand 0 b2*rand;...
%     0 0 a2*rand a2*rand b2*rand 0 b2*rand 0];
%
% Ag = Ag+0.3*F;
% Ag = triu(Ag,1);
% Ag = Ag+Ag';

G = graph(Ag);
[cycles,edgecycles] = allcycles(G);

xycyc = [...
    0 0;...
    1 0;...
    0 0;...
    0 0;...
    0 0;...
    0 0;...
    0 1;...
    0 0;...
    0 0;...
    0 0;...
    0 0;...
    0 0;...
    0 0;...
    0 0;...
    0 -1;...
    1 0;...
    1 0;...
    0 1];


%remove all non-modal cycles
for k1 = size(cycles,1):-1:1
    index = find([cycles{k1,1}] == mode);
    if isempty(index)==1
        cycles(k1,:) = [];
        edgecycles(k1,:) = [];
    else
    end
end


xy = zeros(size(edgecycles,1),2);
gdist = zeros(size(edgecycles,1),1);

for k1 = 1:size(cycles,1)
    for k2 = 1:size(edgecycles{k1,1},2)
        if abs(xycyc((edgecycles{k1,1}(k2)),1)) == 1
            if k2 < size(edgecycles{k1,1},2)
                if cycles{k1,1}(k2) < cycles{k1,1}(k2+1)
                    sgn = 1;
                else
                    sgn = -1;
                end
            else
                if cycles{k1,1}(k2) < cycles{k1,1}(1)
                    sgn = 1;
                else
                    sgn = -1;
                end
            end
            xy(k1,1) = xy(k1,1) + sgn*xycyc(edgecycles{k1,1}(k2),1);
        end
        if abs(xycyc((edgecycles{k1,1}(k2)),2)) == 1
            if k2 < size(edgecycles{k1,1},2)
                if cycles{k1,1}(k2) < cycles{k1,1}(k2+1)
                    sgn = 1;
                else
                    sgn = -1;
                end
            else
                if cycles{k1,1}(k2) < cycles{k1,1}(1)
                    sgn = 1;
                else
                    sgn = -1;
                end
            end
            xy(k1,2) = xy(k1,2) + sgn*xycyc(edgecycles{k1,1}(k2),2);
        end
        gdist(k1,1) = gdist(k1,1) + G.Edges.Weight(edgecycles{k1,1}(k2));
    end
end

%Remove all of the zero-motion cycles
idx = find(xy(:,1)==0);
idy = find(xy(:,2)==0);
id = intersect(idx,idy);
xy(id,:) =[];
gdist(id,:) =[];
cycles(id,:) =[];
edgecycles(id,:) =[];


edist = xy(:,1).^2 + xy(:,2).^2;
edist = edist.^0.5;

gxy = [(1./gdist).*xy(:,1) (1./gdist).*xy(:,2)];
gxy = gxy';

nxy = [(1./edist).*xy(:,1) (1./edist).*xy(:,2)];
nxy = nxy';

U = zeros(1,size(gxy,2));
V = zeros(1,size(gxy,2));

xy = xy';
gdist = gdist';
xy = [xy [-xy(1,:); -xy(2,:)]];
gdist = [gdist; gdist];

nxy = [nxy [-nxy(1,:); -nxy(2,:)]];
gxy = [gxy [-gxy(1,:); -gxy(2,:)]];

nxy(nxy(:,:)==inf) = 0;
nxy(isnan(nxy(:,:))==1) = 0;

gxy(gxy(:,:)==inf) = 0;
gxy(isnan(gxy(:,:))==1) = 0;

U = [U U];
V = [V V];


% h = figure(1);
% q = quiver(U,V,nxy(1,:),nxy(2,:),0);
% q.ShowArrowHead = 'off';
% xlim([-1.2 1.2]);
% ylim([-1.2 1.2]);
% 
% hold on
% scatter(nxy(1,:),nxy(2,:), 40, 'MarkerEdgeColor',[0.5 0 0],...
%     'MarkerFaceColor',[0.7 0 0],...
%     'LineWidth',1.5)
% [k,av] = convhull(nxy(1,:),nxy(2,:));
% plot(nxy(1,k),nxy(2,k));
% 
% q = quiver(U,V,gxy(1,:),gxy(2,:),0);
% q.ShowArrowHead = 'off';
% hold on
% scatter(gxy(1,:),gxy(2,:),40,'MarkerEdgeColor',[0 .5 .5],...
%     'MarkerFaceColor',[0 .7 .7],...
%     'LineWidth',1.5)
% [k,av] = convhull(gxy(1,:),gxy(2,:));
% plot(gxy(1,k),gxy(2,k));
% hold off
% 
% 
% 
% 
% 
% maxsamples = 30000;
% vx = zeros(2,maxsamples);
% px = zeros(2,maxsamples);
% ex = zeros(2,maxsamples);
% for samples = 1:maxsamples
%     for steps = 1:5
%         spin = randi(496);
%         vx(1,samples) = vx(1,samples) + gxy(1,spin);
%         vx(2,samples) = vx(2,samples) + gxy(2,spin);
%         spin = randi(4);
%         if spin ==1
%             px(1,samples) = px(1,samples) + 1;
%         elseif spin ==2
%             px(1,samples) = px(1,samples) - 1;
%         elseif spin ==3
%             px(2,samples) = px(2,samples) + 1;
%         elseif spin ==4
%             px(2,samples) = px(2,samples) - 1;
%         end
%         spin = randi(16);
%         if spin ==1
%             ex(1,samples) = ex(1,samples) + 1;
%         elseif spin ==2
%             ex(1,samples) = ex(1,samples) - 1;
%         elseif spin ==3
%             ex(2,samples) = ex(2,samples) + 1;
%         elseif spin ==4
%             ex(2,samples) = ex(2,samples) - 1;
%         elseif spin ==5
%             ex(2,samples) = ex(2,samples) + 1/sqrt(2);
%             ex(1,samples) = ex(1,samples) + 1/sqrt(2);
%         elseif spin ==6
%             ex(2,samples) = ex(2,samples) + 1/sqrt(2);
%             ex(1,samples) = ex(1,samples) - 1/sqrt(2);
%         elseif spin ==7
%             ex(2,samples) = ex(2,samples) - 1/sqrt(2);
%             ex(1,samples) = ex(1,samples) + 1/sqrt(2);
%         elseif spin ==8
%             ex(2,samples) = ex(2,samples) - 1/sqrt(2);
%             ex(1,samples) = ex(1,samples) - 1/sqrt(2);
%         elseif spin ==9
%             ex(2,samples) = ex(2,samples) + 2/sqrt(5);
%             ex(1,samples) = ex(1,samples) + 1/sqrt(5);
%         elseif spin ==10
%             ex(2,samples) = ex(2,samples) + 2/sqrt(5);
%             ex(1,samples) = ex(1,samples) - 1/sqrt(5);
%         elseif spin ==11
%             ex(2,samples) = ex(2,samples) - 2/sqrt(5);
%             ex(1,samples) = ex(1,samples) + 1/sqrt(5);
%         elseif spin ==12
%             ex(2,samples) = ex(2,samples) - 2/sqrt(5);
%             ex(1,samples) = ex(1,samples) - 1/sqrt(5);
%         elseif spin ==13
%             ex(2,samples) = ex(2,samples) + 1/sqrt(5);
%             ex(1,samples) = ex(1,samples) + 2/sqrt(5);
%         elseif spin ==14
%             ex(2,samples) = ex(2,samples) + 1/sqrt(5);
%             ex(1,samples) = ex(1,samples) - 2/sqrt(5);
%         elseif spin ==15
%             ex(2,samples) = ex(2,samples) - 1/sqrt(5);
%             ex(1,samples) = ex(1,samples) + 2/sqrt(5);
%         elseif spin ==16
%             ex(2,samples) = ex(2,samples) - 1/sqrt(5);
%             ex(1,samples) = ex(1,samples) - 2/sqrt(5);
%         end
%     end
% end
% 
% figure(3)
% scatter(vx(1,:),vx(2,:),1, 'MarkerEdgeColor',[0.5 0 0],...
%     'MarkerFaceColor',[0.7 0 0],...
%     'LineWidth',1.5)
% 
% 
% 
% figure(5)
% hold on
% scatter(px(1,:),px(2,:),40, 'MarkerEdgeColor',[0 0.5 0],...
%     'MarkerFaceColor',[0 0.7 0],...
%     'LineWidth',1.5)
% scatter(ex(1,:),ex(2,:),0.5, 'MarkerEdgeColor',[0 0 0.5],...
%     'MarkerFaceColor',[0 0 0.7],...
%     'LineWidth',1.5)
% hold off
% 
% 
% figure(4)
% histogram(ex(1,:),'FaceColor',[0 0 0.5]...
%     )
% hold on
% histogram(px(1,:),'FaceColor',[0 0.5 0]...
%     )
% histogram(vx(1,:),'FaceColor',[0.5 0 0]...
%     )
% legend
% hold off
% 

%

maxsamples = 10000;
vx = zeros(2,maxsamples);
dx = zeros(1,maxsamples);

px = zeros(2,maxsamples);
dpx = zeros(1,maxsamples);


samplearea = 50;
surveyP = zeros(samplearea,samplearea);
surveyC = ones(samplearea,samplearea);
surveyD = zeros(samplearea,samplearea);
surveyPZ = zeros(samplearea,samplearea);

surveycount = 1;
freq = 1/600;
% freq = 1/4255.3;
numPhi = 0;

for phi = 3*pi/4:pi/8:5*pi/4
    numPhi = numPhi+1;
for samples = 1:maxsamples
%     randsteps = randi(6000);
    for steps = 1:1000
        
        %NSWE spinner
%         spin = randi(4);
%         if spin ==1
%             vx(1,samples) = vx(1,samples) + 1;
%             dx(1,samples) = dx(1,samples)+1;
%         elseif spin ==2
%             vx(1,samples) = vx(1,samples) - 1;
%             dx(1,samples) = dx(1,samples)+1;
%         elseif spin ==3
%             vx(2,samples) = vx(2,samples) + 1;
%             dx(1,samples) = dx(1,samples)+1;
%         elseif spin ==4
%             vx(2,samples) = vx(2,samples) - 1;
%             dx(1,samples) = dx(1,samples)+1;
%         end
        
% %         16 way spinner
%         spin = randi(16);
%         if spin ==1
%             vx(1,samples) = vx(1,samples) + 1;
%             dx(1,samples) = dx(1,samples)+1; 
%         elseif spin ==2
%             vx(1,samples) = vx(1,samples) - 1;
%             dx(1,samples) = dx(1,samples)+1; 
%         elseif spin ==3
%             vx(2,samples) = vx(2,samples) + 1;
%             dx(1,samples) = dx(1,samples)+1; 
%         elseif spin ==4
%             vx(2,samples) = vx(2,samples) - 1;
%             dx(1,samples) = dx(1,samples)+1; 
%         elseif spin ==5
%             vx(2,samples) = vx(2,samples) + 1;
%             vx(1,samples) = vx(1,samples) + 1;
%             dx(1,samples) = dx(1,samples) + 1/sqrt(2); 
%         elseif spin ==6
%             vx(2,samples) = vx(2,samples) + 1;
%             vx(1,samples) = vx(1,samples) - 1;
%             dx(1,samples) = dx(1,samples) + 1/sqrt(2); 
%         elseif spin ==7
%             vx(2,samples) = vx(2,samples) - 1;
%             vx(1,samples) = vx(1,samples) + 1;
%             dx(1,samples) = dx(1,samples) + 1/sqrt(2); 
%         elseif spin ==8
%             vx(2,samples) = vx(2,samples) - 1;
%             vx(1,samples) = vx(1,samples) - 1;
%             dx(1,samples) = dx(1,samples) + 1/sqrt(2); 
% %         elseif spin ==9
% %             vx(2,samples) = vx(2,samples) + 2;
% %             vx(1,samples) = vx(1,samples) + 1;
% %             dx(1,samples) = dx(1,samples) + 1/sqrt(5); 
% %         elseif spin ==10
% %             vx(2,samples) = vx(2,samples) + 2;
% %             vx(1,samples) = vx(1,samples) - 1;
% %             dx(1,samples) = dx(1,samples) + 1/sqrt(5); 
% %         elseif spin ==11
% %             vx(2,samples) = vx(2,samples) - 2;
% %             vx(1,samples) = vx(1,samples) + 1;
% %             dx(1,samples) = dx(1,samples) + 1/sqrt(5); 
% %         elseif spin ==12
% %             vx(2,samples) = vx(2,samples) - 2;
% %             vx(1,samples) = vx(1,samples) - 1;
% %             dx(1,samples) = dx(1,samples) + 1/sqrt(5); 
% %         elseif spin ==13
% %             vx(2,samples) = vx(2,samples) + 1;
% %             vx(1,samples) = vx(1,samples) + 2;
% %             dx(1,samples) = dx(1,samples) + 1/sqrt(5); 
% %         elseif spin ==14
% %             vx(2,samples) = vx(2,samples) + 1;
% %             vx(1,samples) = vx(1,samples) - 2;
% %             dx(1,samples) = dx(1,samples) + 1/sqrt(5); 
% %         elseif spin ==15
% %             vx(2,samples) = vx(2,samples) - 1;
% %             vx(1,samples) = vx(1,samples) + 2;
% %             dx(1,samples) = dx(1,samples) + 1/sqrt(5); 
% %         elseif spin ==16
% %             vx(2,samples) = vx(2,samples) - 1;
% %             vx(1,samples) = vx(1,samples) - 2;
% %             dx(1,samples) = dx(1,samples) + 1/sqrt(5); 
%         end
%            
%         %Supergraph spinner
        spin = randi(496);
        vx(1,samples) = vx(1,samples) + xy(1,spin);
        vx(2,samples) = vx(2,samples) + xy(2,spin);
        dx(1,samples) = dx(1,samples)+ gdist(spin);
% 

        if abs(vx(1,samples))<samplearea && abs(vx(2,samples))<samplearea
            xcoord = abs(vx(1,samples))+1;
            ycoord = abs(vx(2,samples))+1;
            for k4 = 1:size(freq)
                wv = exp(-freq(k4)*dx(1,samples))*sin(phi+2*pi*freq(k4)*dx(1,samples));
            surveyP(xcoord,ycoord) = surveyP(xcoord,ycoord)+ wv;
            end
            surveyC(xcoord,ycoord) = surveyC(xcoord,ycoord)+1;
            surveyD(xcoord,ycoord) = surveyD(xcoord,ycoord)+ dx(1,samples);
        end
    end
end
surveyPZ = surveyPZ+surveyP;
end
surveyPZ = surveyPZ/numPhi;

surveyC(:,:) = max(surveyC(:,:),1);

surveyP(2:end,2:end) = surveyP(2:end,2:end)/4;
surveyD(2:end,2:end) = surveyD(2:end,2:end)/4;
surveyC(2:end,2:end) = surveyC(2:end,2:end)/4;
surveyP(2:end,1) = surveyP(2:end,1)/2;
surveyD(2:end,1) = surveyD(2:end,1)/2;
surveyC(2:end,1) = surveyC(2:end,1)/2;
surveyP(1,2:end) = surveyP(1,2:end)/2;
surveyD(1,2:end) = surveyD(1,2:end)/2;
surveyC(1,2:end) = surveyC(1,2:end)/2;

surveyP2 = surveyP;
surveyP2 =surveyP2./(surveyC);
% surveyP2 = abs(surveyP2);


figure(10)
imagesc(surveyP2)
colormap('jet')
colorbar

figure(11)
surf(surveyP2)

figure(20)
imagesc(surveyD)
colormap('jet')
colorbar

figure(21)
surf(surveyD)


