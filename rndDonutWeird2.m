clear

%declare weights and the mode of operation
a1 =(2-sqrt(2))/2;
a2 = (3-sqrt(5))/2;
b1 = sqrt(2) - 1;
b2 = sqrt(5) - 2;

% a1 = rand;
% a2 = rand;
% b1 = rand;
% b2 = rand;


mode = 2;


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

nxy = [nxy [-nxy(1,:); -nxy(2,:)]];
gxy = [gxy [-gxy(1,:); -gxy(2,:)]];

nxy(nxy(:,:)==inf) = 0;
nxy(isnan(nxy(:,:))==1) = 0;

gxy(gxy(:,:)==inf) = 0.1;
gxy(isnan(gxy(:,:))==1) = 0.1;

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

maxsamples = 20000;
vx = zeros(2,maxsamples);
px = zeros(2,maxsamples);
for samples = 1:maxsamples
    for steps = 1:1000
        spin = randi(496);
        vx(1,samples) = vx(1,samples) + gxy(1,spin);
        vx(2,samples) = vx(2,samples) + gxy(2,spin);
        spin = randi(4);
        if spin ==1
            px(1,samples) = px(1,samples) + 1/pi;
        elseif spin ==2
            px(1,samples) = px(1,samples) - 1/pi;
        elseif spin ==3
            px(2,samples) = px(2,samples) + 1/pi;
        elseif spin ==4
            px(2,samples) = px(2,samples) - 1/pi;
        end
    end
end

figure(3)
xlim([-15 15]);
ylim([-15 15]);
scatter(vx(1,:),vx(2,:))
hold on
scatter(px(1,:),px(2,:))
legend
hold off

figure(4)
histogram(vx())
hold on
histogram(px)
hold off






