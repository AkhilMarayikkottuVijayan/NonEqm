import numpy as np 
import matplotlib.pyplot as plt

# Define the functions we will be using

def m(t):
    if t <= tb1:
        return m0 - m_dot*t 
    else:
        return mpl

def cL(alpha):
    return -0.04 + 0.8 * alpha

def cD(alpha):
    return 0.012 - 0.01 * alpha + 0.6 * (alpha**2)

# Drag for launch vehicle 
def D(t,vrx,vry):
    
    if t <= tb1 + 75:
        return 0.5*rho*(vrx**2+vry**2)*S*cd
    else:
        return 0.5*rho*(vrx**2+vry**2)*0.015*cD(thetaf(t,vrx,vry))
    
# Lift for wave-glider
def L(t,vrx,vry):
    
    if t<= tb1 + 75:
        return 0
    else:
        return 0.5*rho*(vrx**2+vry**2)*0.015*cL(thetaf(t,vrx,vry))
    
def thetaf(t,vrx,vry):
    if t <= trail:
        return np.deg2rad(theta1)
    
    elif t <= 0.8406*tb1:
        return np.arctan(vry/vrx) - ((np.deg2rad(0.00325)*t)/0.01)
    
    elif tb1+75 < t <= tb1 + 75.72:
        return np.arctan(vry/vrx) + np.deg2rad(-4.1*(t-(tb1+75)))
    
    elif tb1+75.72 < t <= tgf:
        
        return np.deg2rad(-10)
    
    else:
        return np.arctan(vry/vrx)
    
        
def accel(dvrxdt,dvrydt,t,vrx,vry):
    
    if t <= trail:
        
        dvrxdt=(T/m(t)- D(t,vrx,vry)/m(t)-g*np.sin(thetaf(t,vrx,vry)))*np.cos(thetaf(t,vrx,vry))
        dvrydt=(T/m(t)- D(t,vrx,vry)/m(t)-g*np.sin(thetaf(t,vrx,vry)))*np.sin(thetaf(t,vrx,vry))    
    
    elif t <= tb1:
        
        dvrxdt = (T/m(t)-D(t,vrx,vry)/m(t))*np.cos(thetaf(t,vrx,vry))
        dvrydt = (T/m(t)-D(t,vrx,vry)/m(t))*np.sin(thetaf(t,vrx,vry)) - g
        
#Problem Section
        
    elif tb1 + 75 < t <= tgf:
        
        dvrxdt = (L(t,vrx,vry)*np.sin(thetaf(t,vrx,vry)) - (D(t,vrx,vry)*np.cos(thetaf(t,vrx,vry))))
        dvrydt = (L(t,vrx,vry)*np.cos(thetaf(t,vrx,vry)+np.deg2rad(15)) + D(t,vrx,vry)*np.sin(thetaf(t,vrx,vry)+np.deg2rad(15)))/(m(t)*g)
        
# Rest Works
        
    else:
        
        dvrxdt=(0-D(t,vrx,vry)/m(t))*np.cos(thetaf(t,vrx,vry))
        dvrydt=(0-D(t,vrx,vry)/m(t))*np.sin(thetaf(t,vrx,vry)) - g
    
    return dvrxdt,dvrydt


#Define the time-stepping method

def rk4(accel, t0, dvxdt0, dvydt0, vx0, vy0, x0, y0, theta0, tf, n):
    
# make the arrys we need
    dthrust = [0] * (n + 1)
    t2 = [0] * (n + 1)
    dm = [0] * (n + 1)
    vx = [0] * (n + 1)
    vy = [0] * (n + 1)
    v = [0] * (n + 1)
    Mach = [0] * (n + 1)
    dvxdt = [0] * (n + 1)
    dvydt = [0] * (n + 1)
    dx = [0] * (n + 1)
    dy = [0] * (n + 1)
    dtheta = [0] * (n + 1)
    Drag = [0] * (n + 1)
    Te = [0] * (n + 1)
    P = [0] * (n + 1)
    rho = [0] * (n + 1)
    sound = [0] * (n + 1)
    Q = [0] * (n + 1)
    Q1 = [0] *(n + 1)
    P1 = [0] * (n + 1)
    T1 = [0] * (n + 1)
    mu = [0] * (n + 1)
    dalpha = [0] * (n + 1)
    dcL = [0] * (n + 1)
    dcD = [0] * (n + 1)
    Lift = [0] * (n + 1)
    
# Define our time step
    h = (tf - t0) / float(n)
    
# Set initial values for Runge-Kutta 4th order
    
    dthrust[0] = T/m0
    dm[0] = m0
    t2[0] = t = t0
    vx[0] = vrx = vx0
    vy[0] = vry = vy0
    dvxdt[0] = dvrxdt = dvxdt0
    dvydt[0] = dvrydt = dvydt0
    dx[0] = x = x0
    dy[0] = y = y0
    dtheta[0] = np.rad2deg(theta0)
    Drag[0] = 0
    v[0] = 0
    Mach[0] = 0
    Te[0] = T0 = 288.15
    P[0] = P0 = 101.719
    rho[0] = 1.225
    sound[0] = 1830
    P1[0] = 0
    T1[0] = 0
    Q[0] = 0
    Q1[0] = 0
    mu[0] = 1.716e-5
    dalpha[0] = 0 
    dcL[0] = 0
    dcD[0] = 0
    Lift[0] = 0
    
# Now we use the method
    
    for i in range(1, n + 1):
        
        
        
# Define each k value for x and y 
        
        k1x = h * accel(dvrxdt, dvrydt, t, vrx, vry)[0]
        k1y = h * accel(dvrxdt, dvrydt, t, vrx, vry)[1]
        k2x = h * accel(dvrxdt,dvrydt, t + 0.5 * h, vrx + 0.5 * k1x,vry + 0.5 * k1y)[0]
        k2y = h * accel(dvrxdt,dvrydt, t + 0.5 * h, vrx + 0.5 * k1x,vry + 0.5 * k1y)[1]
        k3x = h * accel(dvrxdt,dvrydt, t + 0.5 * h, vrx + 0.5 * k2x,vry + 0.5 * k2y)[0]
        k3y = h * accel(dvrxdt,dvrydt, t + 0.5 * h, vrx + 0.5 * k2x,vry + 0.5 * k2y)[1]
        k4x = h * accel(dvrxdt,dvrydt, t + h, vrx + k3x, vry + k3y)[0]
        k4y = h * accel(dvrxdt,dvrydt, t + h, vrx + k3x, vry + k3y)[1]
         
#Return these variables to use and plot 
        
        t2[i] = t = t0 + i * h
        vx[i] = vrx = vrx + ((k1x + 2*k2x + 2*k3x + k4x) / 6)
        vy[i] = vry = vry + ((k1y + 2*k2y + 2*k3y + k4y) / 6)
        dx[i] = x = x + h*(vrx +((k1x + k2x + k3x) / 6)) 
        dy[i] = y = y + h*(vry +((k1y + k2y + k3y) / 6)) 
        v[i] = np.sqrt(vx[i]**2+vy[i]**2)
        dtheta[i] = np.rad2deg(thetaf(t,vrx,vry))
        dcL[i] = cL(dtheta[i])
        dcD[i] = cD(dtheta[i])
        dvxdt[i] = accel(dvrxdt, dvrydt, t, vrx, vry)[0]
        dvydt[i] = accel(dvrxdt, dvrydt, t, vrx, vry)[1]
        
#Fill the previously defined arrays with the values
        
## Determine Pressure/Temperature/Density as a function of altitude
        if dy[i] <= 11000:
                Te[i] = T0- 0.00649*dy[i]
                P[i] = P0*((Te[i])/288.08)**5.256
        elif 11000 < dy[i] <= 25000:
                Te[i] = -56.46 +273.15
                P[i] = 22.65*np.exp(1.73-0.000157*dy[i])
        else:
                Te[i] = (-131.21 + 0.00299*dy[i])+273.15
                P[i] = 2.488*((Te[i])/216.6)**-11.388
                
# Other variables to be plotted
        rho[i] = P[i]/(0.2869*Te[i]) #[kg/m^3]
        sound[i] = (np.sqrt((1.4*8314*Te[i])))
        Mach[i] = v[i]/sound[i]
        mu[i] = mu[0]*((Te[i]/Te[0])**1.5)*((Te[0]+110.4)/(Te[i]+110.4))
       
        
# Plot the Heat Fluxes
        
# Scott Model        
        Q[i] = 1.1813e-3*np.sqrt(rho[i])*v[i]**3.05
        
        
# Here we plot D(vrx,vry)/m(t)
        
        if t2[i] <= tb1:
            dm[i] = m0 - m_dot*t2[i]    
        else:
             dm[i] = mpl
             
        if t2[i] <= tb1:
            dthrust[i] = T/dm[i]
        else:
            dthrust[i] = 0
        
        Drag[i] = D(t2[i],vx[i],vy[i])     
        Lift[i] = L(t2[i],vx[i],vy[i])  
        
    return Lift, dcL,dcD,dalpha,Q1, Q, Te, P, rho, sound, dthrust, t2, dm, dvxdt, dvydt, v, vx, vy, Mach, dx, dy, dtheta, Drag
    
         
Lift, dcL, dcD, dalpha, Q1, Q, Te, P, rho, sound, dthrust, t2, dm, dvxdt, dvydt, v, vx, vy, Mach, dx, dy, dtheta, Drag = rk4(accel, 0, 0, 0, 0, 0, 0, 0, np.deg2rad(theta1), 1600, 16000)
