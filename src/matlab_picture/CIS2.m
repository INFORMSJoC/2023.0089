clear
clc
set(0,'defaultfigurecolor','w');
filename = '.\self.txt50.txt';
[nu,x,y] = textread(filename,'%n%n%n') ;
Ncircles = nu(1) ;
R = x(1) ; 
Circle  = zeros(Ncircles,2) ;
contact = zeros(Ncircles,1);
for i=1:Ncircles
    Circle(i,1) = x(i+1);
    Circle(i,2) = y(i+1); 
end 
figure;
hold on;
filename = '.\self.txt';
[p,q] = textread(filename,'%n%n') ;
nenda = 0.1;

numP = q(1);
 i=1:numP;
 x = p(i+1);
 y = q(i+1);
 fill(x,y,'w');

n = numP+3;
while p(n) ~= 404
    numP = q(n);
    i=n+1:n+numP;
    x = p(i);
    y = q(i);
    fill(x,y,[0.7,0.7,0.7]);
    n = n+numP+2;
end


numP = q(1);
for i=1:numP
    line([p(i+1),p(i+2)],[q(i+1),q(i+2)], 'linewidth',0.5, 'linestyle','-','color','k');
end 
n = numP+3;
while p(n) ~= 404
    numP = q(n);
    for i=n+1:n+numP
        line([p(i),p(i+1)],[q(i),q(i+1)], 'linewidth',0.5, 'linestyle','-','color','k');
    end
    n = n+numP+2;
end

axis equal
axis off
epslon = 1.0e-8; 

for i = 1:Ncircles
    x_i = Circle(i,1);
    y_i = Circle(i,2);
    numP = q(1);
    for j=2:numP+1
        dist = sqrt( (x_i-p(j))*(x_i-p(j)) + (y_i - q(j))*(y_i - q(j)) );
        if dist < R*nenda + epslon
            contact(i) =  contact(i)+1;
        end
        %圆和容器顶点
        if p(j+1) == p(j) %线段方程为Ax+By+C=0 �?x-x1=0
            a = 1.0;
            b = 0;
            c = -p(j);
        else                              %线段方程为Ax+By+C=0 �?kx-y+y1-kx1=0
            a = (q(j+1) - q(j)) / (p(j+1) - p(j));%k
            b = -1.0;
            c = q(j) - a*p(j);%y1-kx1
        end
        cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
        cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
        if p(j)>p(j+1)
            maxx = p(j);
            minx = p(j+1);
        else
            maxx = p(j+1);
            minx = p(j);            
        end
        if q(j)>q(j+1)
            maxy = q(j);
            miny = q(j+1);
        else
            maxy = q(j+1);
            miny = q(j);            
        end
        if cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny %垂点在边�?
            dist = sqrt( (x_i-cdx)*(x_i-cdx) + (y_i - cdy)*(y_i - cdy) );
            if dist < R*nenda + epslon
                contact(i) =  contact(i)+1;
            end
        end
        %圆和容器�?
    end
    n = numP+3;
    while p(n) ~= 404
        numP = q(n);
        for k=n+1:n+numP
            dist = sqrt( (x_i-p(k))*(x_i-p(k)) + (y_i - q(k))*(y_i - q(k)) );
            if dist < R*nenda + epslon
                contact(i) =  contact(i)+1;
            end
            %圆和漏洞顶点
            if p(k+1) == p(k) %线段方程为Ax+By+C=0 �?x-x1=0
                a = 1.0;
                b = 0;
                c = -p(k);
            else                              %线段方程为Ax+By+C=0 �?kx-y+y1-kx1=0
                a = (q(k+1) - q(k)) / (p(k+1) - p(k));%k
                b = -1.0;
                c = q(k) - a*p(k);%y1-kx1
            end
            cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
            cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
            if p(k)>p(k+1)
                maxx = p(k);
                minx = p(k+1);
            else
                maxx = p(k+1);
                minx = p(k);
            end
            if q(k)>q(k+1)
                maxy = q(k);
                miny = q(k+1);
            else
                maxy = q(k+1);
                miny = q(k);
            end
            if cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny %垂点在边�?
                dist = sqrt( (x_i-cdx)*(x_i-cdx) + (y_i - cdy)*(y_i - cdy) );
                if dist < R*nenda + epslon
                    contact(i) =  contact(i)+1;
                end
            end
            %圆和漏洞�?
        end
        n = n+numP+2;
    end
end


for i = 1:Ncircles-1
    for j = i+1:Ncircles
        if sqrt((Circle(i,1)-Circle(j,1))*(Circle(i,1)-Circle(j,1)) + (Circle(i,2)-Circle(j,2))*(Circle(i,2)-Circle(j,2)))  < 2.0*R + epslon
          contact(i) =  contact(i)+1;
          contact(j) =  contact(j)+1;
        end 
    end
end

t= 0:0.001:2.0*pi;
for i =1:Ncircles
    c1= R*0.3*cos(t)+Circle(i,1);
    c2= R*0.3*sin(t)+Circle(i,2); 
    plot(c1,c2,'-k'); 
    switch contact(i)
    case 0
         fill(c1,c2,[1 0 0]); %[1 0 0] is RGB of color
    case 1
          fill(c1,c2,[255/255.0 105/255.0 180/255.0]);   
    case 2
          fill(c1,c2,[0 191/255.0 255/255.0]); 
    case  3
          fill(c1,c2,[0.667 0.667 1]);   
    case  4
          fill(c1,c2,'g');  
    case  5
          fill(c1,c2,[238/255.0 221/255.0 130/255.0]);   
    otherwise 
          fill(c1,c2,[255/255 165/255 0]);          
    end 
end 

for i = 1:Ncircles-1
    for j = i+1:Ncircles
        if sqrt((Circle(i,1)-Circle(j,1))*(Circle(i,1)-Circle(j,1)) + (Circle(i,2)-Circle(j,2))*(Circle(i,2)-Circle(j,2)))  < 2.0*R + epslon
          line([Circle(i,1),Circle(j,1)],[Circle(i,2),Circle(j,2)],'linestyle','-','color','r');
        end
    end
end

for i = 1:Ncircles
    x_i = Circle(i,1);
    y_i = Circle(i,2);
    numP = q(1);
    for j = 2:numP+1
        dist = sqrt( (x_i-p(j))*(x_i-p(j)) + (y_i - q(j))*(y_i - q(j)) );
        if dist < R*nenda + epslon
            line([Circle(i,1),p(j)],[Circle(i,2),q(j)],'linestyle','-','color','r'); 
        end
        %圆和容器顶点
        if p(j+1) == p(j) %线段方程为Ax+By+C=0 �?x-x1=0
            a = 1.0;
            b = 0;
            c = -p(j);
        else                              %线段方程为Ax+By+C=0 �?kx-y+y1-kx1=0
            a = (q(j+1) - q(j)) / (p(j+1) - p(j));%k
            b = -1.0;
            c = q(j) - a*p(j);%y1-kx1
        end
        cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
        cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
        if p(j)>p(j+1)
            maxx = p(j);
            minx = p(j+1);
        else
            maxx = p(j+1);
            minx = p(j);            
        end
        if q(j)>q(j+1)
            maxy = q(j);
            miny = q(j+1);
        else
            maxy = q(j+1);
            miny = q(j);            
        end
        if cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny %垂点在边�?
            dist = sqrt( (x_i-cdx)*(x_i-cdx) + (y_i - cdy)*(y_i - cdy) );
            if dist <= R*nenda + epslon
                line([Circle(i,1),cdx],[Circle(i,2),cdy],'linestyle','-','color','r'); 
            end
        end
        %圆和容器�?
    end
    n = numP+3;
    while p(n) ~= 404
        numP = q(n);
        for k = n+1:n+numP
            dist = sqrt( (x_i-p(k))*(x_i-p(k)) + (y_i - q(k))*(y_i - q(k)) );
            if dist < R*nenda + epslon
                line([Circle(i,1),p(k)],[Circle(i,2),q(k)],'linestyle','-','color','r'); 
            end
            %圆和漏洞顶点
            if p(k+1) == p(k) %线段方程为Ax+By+C=0 �?x-x1=0
                a = 1.0;
                b = 0;
                c = -p(k);
            else                              %线段方程为Ax+By+C=0 �?kx-y+y1-kx1=0
                a = (q(k+1) - q(k)) / (p(k+1) - p(k));%k
                b = -1.0;
                c = q(k) - a*p(k);%y1-kx1
            end
            cdx = (b*b*x_i - a*b*y_i - a*c)/(a*a+b*b);
            cdy = (a*a*y_i - a*b*x_i - b*c)/(a*a+b*b);
            if p(k)>p(k+1)
                maxx = p(k);
                minx = p(k+1);
            else
                maxx = p(k+1);
                minx = p(k);
            end
            if q(k)>q(k+1)
                maxy = q(k);
                miny = q(k+1);
            else
                maxy = q(k+1);
                miny = q(k);
            end
            if cdx <= maxx && cdx >= minx && cdy <= maxy && cdy >= miny %垂点在边�?
                dist = sqrt( (x_i-cdx)*(x_i-cdx) + (y_i - cdy)*(y_i - cdy) );
                if dist < R*nenda + epslon
                    line([Circle(i,1),cdx],[Circle(i,2),cdy],'linestyle','-','color','r');
                end
            end
            %圆和漏洞�?
        end
        n = n+numP+2;
    end
end

