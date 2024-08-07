function DCrr(x, y, r)
 alpha=0:pi/20:2*pi;%角度[0,2*pi]
 px=r*cos(alpha)+x;
 py=r*sin(alpha)+y;
 plot(px,py,'linewidth',3,'color','red');
 hold on
end