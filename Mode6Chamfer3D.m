clear
tic
a1 = (2-sqrt(2))/2;
a2 = 1;
A = [a2 0.00001 0.00001 0.00001];
A2 =[a2 0.00001 0.00001 0.00001];
b1 = (sqrt(2)-1)/2;
b1 = 0.00001;
b2 = sqrt(2)*b1;

d1 = [0.5 0.00001 0.00001 0.00001];
d2 = [sqrt(2)-0.5 0.00001 0.00001 0.00001];
d3 = [sqrt(2)-1 0.00001 0.00001 0.00001];


Bx = [b1 0.5 0.00001 0.00001];
By = [b1 0.00001 0.5 0.00001];
Bz = [b1 0.00001 0.00001 0.5];

% Bx = [b1 1 0.00001 0.00001];
% By = [sqrt(2)-1 0.5 sqrt(3)/2 0.00001];
% Bz = [sqrt(2)-1 -0.5 sqrt(3)/2 0.00001 1];

% Bxy = [b1 0.00001 0.00001 1];
% Byz = [b1 0.00001 0.00001 1];
% Bzx = [b1 0.00001 0.00001 1];

% 
% Cx = [b1/2 1 0.00001 0.00001];
% Cy = [b1/2 0.00001 1 0.00001];
% Cz = [b1/2 0.00001 0.00001 1];


mode = 3;



k = 1;
Adj = [ ...
    0 d1(k) d1(k) A(k) d1(k) d1(k) Bx(k) 0 0; ...
    d1(k) 0 A(k) d2(k) A(k) A(k) 0 By(k) 0; ...
    d1(k) A(k) 0 d2(k) d3(k) A(k) 0 0 Bz(k); ...
    A(k) d2(k) d2(k) 0 d2(k) d2(k) Bx(k) 0 0; ...
    d1(k) A(k) d3(k) d2(k) 0 d3(k) 0 By(k) 0; ...
    d1(k) A(k) A(k) d2(k) d3(k) 0 0 0 Bz(k); ...
    Bx(k) 0 0 Bx(k) 0 0 0 0 0; ...
    0 By(k) 0 0 By(k) 0 0 0 0; ...
    0 0 Bz(k) 0 0 Bz(k) 0 0 0];

% Adj = [ ...
%     0 A(k) A(k) Bx(k) A(k) A(k); ...
%     A(k) 0 A(k) A(k) By(k) A(k); ...
%     A(k) A(k) 0 A(k) A(k) Bz(k); ...
%     Bx(k) A(k) A(k) 0 A(k) A(k); ...
%     A(k) By(k) A(k) A(k) 0 A(k); ...
%     A(k) A(k) Bz(k) A(k) A(k) 0];


G = graph(Adj);
[c1,e1] = allcycles(G);
G1 = G.Edges.Weight;
toc

k = 2;
Adj = [ ...
    0 d1(k) d1(k) A(k) d1(k) d1(k) Bx(k) 0 0; ...
    d1(k) 0 A(k) d2(k) A(k) A(k) 0 By(k) 0; ...
    d1(k) A(k) 0 d2(k) d3(k) A(k) 0 0 Bz(k); ...
    A(k) d2(k) d2(k) 0 d2(k) d2(k) Bx(k) 0 0; ...
    d1(k) A(k) d3(k) d2(k) 0 d3(k) 0 By(k) 0; ...
    d1(k) A(k) A(k) d2(k) d3(k) 0 0 0 Bz(k); ...
    Bx(k) 0 0 Bx(k) 0 0 0 0 0; ...
    0 By(k) 0 0 By(k) 0 0 0 0; ...
    0 0 Bz(k) 0 0 Bz(k) 0 0 0];
G = graph(Adj);
[c2,e2] = allcycles(G);
G2 = G.Edges.Weight;

k = 3;
Adj = [ ...
    0 d1(k) d1(k) A(k) d1(k) d1(k) Bx(k) 0 0; ...
    d1(k) 0 A(k) d2(k) A(k) A(k) 0 By(k) 0; ...
    d1(k) A(k) 0 d2(k) d3(k) A(k) 0 0 Bz(k); ...
    A(k) d2(k) d2(k) 0 d2(k) d2(k) Bx(k) 0 0; ...
    d1(k) A(k) d3(k) d2(k) 0 d3(k) 0 By(k) 0; ...
    d1(k) A(k) A(k) d2(k) d3(k) 0 0 0 Bz(k); ...
    Bx(k) 0 0 Bx(k) 0 0 0 0 0; ...
    0 By(k) 0 0 By(k) 0 0 0 0; ...
    0 0 Bz(k) 0 0 Bz(k) 0 0 0];
G = graph(Adj);
[c3,e3] = allcycles(G);
G3 = G.Edges.Weight;

k = 4;
Adj = [ ...
    0 d1(k) d1(k) A(k) d1(k) d1(k) Bx(k) 0 0; ...
    d1(k) 0 A(k) d2(k) A(k) A(k) 0 By(k) 0; ...
    d1(k) A(k) 0 d2(k) d3(k) A(k) 0 0 Bz(k); ...
    A(k) d2(k) d2(k) 0 d2(k) d2(k) Bx(k) 0 0; ...
    d1(k) A(k) d3(k) d2(k) 0 d3(k) 0 By(k) 0; ...
    d1(k) A(k) A(k) d2(k) d3(k) 0 0 0 Bz(k); ...
    Bx(k) 0 0 Bx(k) 0 0 0 0 0; ...
    0 By(k) 0 0 By(k) 0 0 0 0; ...
    0 0 Bz(k) 0 0 Bz(k) 0 0 0];
G = graph(Adj);
[c4,e4] = allcycles(G);
G4 = G.Edges.Weight;

toc

F = [ ...
    0 0 0 0 0 0 1 0 0; ...
    0 0 0 0 0 0 0 1 0; ...
    0 0 0 0 0 0 0 0 1; ...
    0 0 0 0 0 0 -1 0 0; ...
    0 0 0 0 0 0 0 -1 0; ...
    0 0 0 0 0 0 0 0 -1; ...
    -1 0 0 1 0 0 0 0 0; ...
    0 -1 0 0 1 0 0 0 0; ...
    0 0 -1 0 0 1 0 0 0; ...
    ];


for k1 = size(c1,1):-1:1
    index = find([c1{k1,1}] == mode);
    if isempty(index)==1
        c1(k1,:) = [];
        e1(k1,:) = [];
        c2(k1,:) = [];
        e2(k1,:) = [];
        c3(k1,:) = [];
        e3(k1,:) = [];
        c4(k1,:) = [];
        e4(k1,:) = [];
    else
    end
end

toc
for k1 = 1:size(e1,1)
    clength(k1) = 0;
    xlength(k1) = 0;
    ylength(k1) = 0;
    zlength(k1) = 0;
    ec = e1{k1};
    ex = e2{k1};
    ey = e3{k1};
    ez = e4{k1};
    cx = c2{k1};
    cy = c3{k1};
    cz = c4{k1};
    for k2 = 1:size(e1{k1},2)
        xdir = 1;
        ydir = 1;
        zdir = 1;
        if k2<size(e1{k1},2)
            xdir = F(cx(k2),cx(k2+1));
            ydir = F(cy(k2),cy(k2+1));
            zdir = F(cz(k2),cz(k2+1));
        else
            xdir = F(cx(k2),cx(1));
            ydir = F(cy(k2),cy(1));
            zdir = F(cz(k2),cz(1));
        end
        clength(k1) = clength(k1)+G1(ec(k2));
        xlength(k1) = xlength(k1)+xdir*G2(ex(k2));
        ylength(k1) = ylength(k1)+ydir*G3(ey(k2));
        zlength(k1) = zlength(k1)+zdir*G4(ez(k2));
    end
end
toc

clength(k1) = clength(k1)';
xlength(k1) = xlength(k1)';
ylength(k1) = ylength(k1)';
zlength(k1) = zlength(k1)';

nlength = sqrt(xlength.^2+ylength.^2+zlength.^2);
gxy = [(1./clength).*xlength; (1./clength).*ylength; (1./clength).*zlength];
gxy = [gxy -gxy];
U = zeros(1,size(gxy,2));
V = zeros(1,size(gxy,2));
W = zeros(1,size(gxy,2));
gxy = gxy';
U = U';
V = V';
W = W';

P = gxy;
% P = P*sqrt(2)/sqrt(3);
[k,vol] = convhulln(P);

figure(1)
[x y z] = sphere(128);
h = surfl(x, y, z);
set(h, 'FaceAlpha', 0.5)
shading flat
hold on
scatter3(P(:,1),P(:,2),P(:,3),70,'filled')
trisurf(k,P(:,1),P(:,2),P(:,3),'FaceColor','cyan')
hold off
% q = quiver(U,V,gxy(:,2),gxy(:,3),0);
% q.ShowArrowHead = 'off';
% hold on
% scatter3(gxy(:,1),gxy(:,2),40,'MarkerEdgeColor',[0 .5 .5],...
%     'MarkerFaceColor',[0 .7 .7],...
%     'LineWidth',1.5)
%
% figure(3)
% trisurf(k,P(:,1),P(:,2),P(:,3),'FaceColor','cyan')
