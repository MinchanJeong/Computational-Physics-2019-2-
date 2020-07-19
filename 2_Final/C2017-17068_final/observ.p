set term png
set term png size 640,480
set output 'obs.png'

f(x) = sin( 10 * (x**2 + 1)/(x**4 + x**2 + 1) * exp(-0.5*x**2))


g(x) = 10.0*x/(1.0+x*x);
h(x) = 10.0*x*x*exp(-1.0*x*x)/(1.0+x*x*x*x);
f(x) = sin(g(x)*sin(h(x)));
        
set title "Solution observation"
set title font "times new roman,20"

set nokey

set xrange[0.0:2.0]
set yrange[-1.5:1.5]
set xtics 0.0,0.1,2
set tics font "times new roman,10"

set zeroaxis
set grid

plot f(x) lt rgb "red"
