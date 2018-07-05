function [yy]=linear_interp(y0,y1,x0,x1,xx);

a=(y1-y0)/(x1-x0);
yy=y0+a*(xx-x0);