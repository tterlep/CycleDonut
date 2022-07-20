clear


%declare weights and the mode of operation
a1 =(2-sqrt(2))/2;
a2 = (3-sqrt(5))/2;
b1 = sqrt(2) - 1;
b2 = sqrt(5) - 2;

% a1 = 1;
% a2 = 0.1;
% b1 = 0.2;
% b2 = 0.3;


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

xy = [xy [-xy(1,:); -xy(2,:)]];
nxy = [nxy [-nxy(1,:); -nxy(2,:)]];
gxy = [gxy [-gxy(1,:); -gxy(2,:)]];

nxy(nxy(:,:)==inf) = 0;
nxy(isnan(nxy(:,:))==1) = 0;

gxy(gxy(:,:)==inf) = 0;
gxy(isnan(gxy(:,:))==1) = 0;

U = [U U];
V = [V V];


h = figure(1);
q = quiver(U,V,nxy(1,:),nxy(2,:),0);
q.ShowArrowHead = 'off';
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);

hold on
scatter(nxy(1,:),nxy(2,:), 40, 'MarkerEdgeColor',[0.5 0 0],...
    'MarkerFaceColor',[0.7 0 0],...
    'LineWidth',1.5)
[k,av] = convhull(nxy(1,:),nxy(2,:));
plot(nxy(1,k),nxy(2,k));

q = quiver(U,V,gxy(1,:),gxy(2,:),0);
q.ShowArrowHead = 'off';
hold on
scatter(gxy(1,:),gxy(2,:),40,'MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[0 .7 .7],...
    'LineWidth',1.5)
[k,av] = convhull(gxy(1,:),gxy(2,:));
plot(gxy(1,k),gxy(2,k));
hold off
% 
% maxsamples = 10000;
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

% figure(3)
% scatter(vx(1,:),vx(2,:),1, 'MarkerEdgeColor',[0.5 0 0],...
%     'MarkerFaceColor',[0.7 0 0],...
%     'LineWidth',1.5)
% 
% 
% figure(5)
% close
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

maxsamples = 100000;
vx = zeros(2,maxsamples);
dx = zeros(1,maxsamples);

px = zeros(2,maxsamples);
dpx = zeros(1,maxsamples);

surveyX = 60;
surveyY = 2;
surveyP = [];
surveycount = 1;
for samples = 1:maxsamples
    for steps = 1:5000
        spin = randi(496);
        vx(1,samples) = vx(1,samples) + xy(1,spin);
        vx(2,samples) = vx(2,samples) + xy(2,spin);
        dx(1,samples) = dx(1,samples)+ sqrt(gxy(1,spin).^2+gxy(2,spin).^2);
        if (vx(1,samples) == surveyX) && (vx(2,samples) == surveyY)
            surveyP(surveycount) = dx(1,samples);
            surveycount = surveycount+1;
        end
        spin = randi(4);
        if spin ==1
            px(1,samples) = px(1,samples) + 1;
            dpx(1,samples) = dpx(1,samples) + 1;
        elseif spin ==2
            px(1,samples) = px(1,samples) - 1;
            dpx(1,samples) = dpx(1,samples) + 1;
        elseif spin ==3
            px(2,samples) = px(2,samples) + 1;
            dpx(1,samples) = dpx(1,samples) + 1;
        elseif spin ==4
            px(2,samples) = px(2,samples) - 1;
            dpx(1,samples) = dpx(1,samples) + 1;
        end
    end
end

figure(3)
scatter(px(1,:),px(2,:),10, 'MarkerEdgeColor',[0 0 0.5],...
    'MarkerFaceColor',[0 0 0.7],...
    'LineWidth',1.5)
hold on
scatter(vx(1,:),vx(2,:),1, 'MarkerEdgeColor',[0.5 0 0],...
    'MarkerFaceColor',[0.7 0 0],...
    'LineWidth',1.5)
hold off

figure(4)
histogram(dx)
hold on
histogram(dpx)
hold off

figure(5)
hold on
histogram(surveyP,40)


