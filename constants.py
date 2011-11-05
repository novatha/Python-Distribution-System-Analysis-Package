from numpy import sin, cos, pi, array
a=cos(2*pi/3)+1j*sin(2*pi/3)
A=array([1.0,a,a*a])# creating three 120-shifted unity voltage.
P=array([[1,1,1],[1,a,a*a],[1,a*a,a]])/3
Pinv=array([[1,1,1],[1,a*a,a],[1,a,a*a]])
ps=cos(pi/6)+1j*sin(pi/6) # this is a 30 degree phase shift
PS=array([[ps,0,0],[0,ps,0],[0,0,ps]])


