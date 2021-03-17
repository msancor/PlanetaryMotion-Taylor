PROGRAM kepler
REAL::a,b,h
REAL::rp(0:20000),rp1(0:20000),rp2(0:20000),t(0:20000)
REAL::r(0:20000),r1(0:20000),r2(0:20000)
REAL::theta(0:20000),theta1(0:20000),theta2(0:20000)
REAL::x(0:20000),y(0:20000)
INTEGER::k,n

!We define our initial conditions and parameters

a=0.0
b=20.0
h=0.001
n=INT((b-a)/h)
t(0)=a
rp(0)=0.2
r(0)=1.0
theta(0)=0

!We use the second order Taylor method for R and theta

DO k=0,n
rp1(k)=1/((r(k))**3) - 1/((r(k))**2)
r1(k)=rp(k)
rp2(k)=((-3.0)/(r(k))**4)*rp(k)+((2.0)/(r(k))**3)*rp(k)
r2(k)=((1.0)/(r(k))**3)-((1.0)/(r(k))**2)
rp(k+1)=rp(k)+h*(rp1(k)+(h/2.0)*rp2(k))
r(k+1)=r(k)+h*(r1(k)+(h/2.0)*r2(k))
t(k+1)=t(k)+h
END DO


DO k=0,n
theta1(k)=((1.0)/((r(k))**2))
theta2(k)=((-2.0)/((r(k))**3))*rp(k)
theta(k+1)=theta(k)+h*(theta1(k)+(h/2.0)*theta2(k))
END DO

!We save the data for the trajectories in polar and cartesian coordinates

OPEN(1,FILE="elipse.txt",STATUS="new")
OPEN(2,FILE="r(t).txt",STATUS="new")
OPEN(3,FILE="theta(t).txt",STATUS="new")
OPEN(4,FILE="polares.txt",STATUS="new")
OPEN(5,FILE="rp(t).txt",STATUS="new")
DO k=0,n
x(k)=r(k)*COS(theta(k))
y(k)=r(k)*SIN(theta(k))
WRITE(1,*) x(k),y(k)
WRITE(2,*) t(k),r(k)
WRITE(3,*) t(k),theta(k)
WRITE(4,*) theta(k),r(k)
WRITE(5,*) t(k),rp(k)
END DO
CLOSE(1)
CLOSE(2)
CLOSE(3)
CLOSE(4)
CLOSE(5)

END PROGRAM
