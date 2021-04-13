import numpy as np 
import matplotlib.pyplot as plt

# Define Initial variables that determine how the boost vehicle performs
g = 9.81
mp1 = 12280 # mass propellant 
mpl = 1200 # kg
m0 = mp1+mpl #kg
T = 300000 #N
Isp = 237 #s
m_dot = T/(Isp*g)
cd = 0.75
d = 0.15
S = (np.pi/4)*d**2
trail = 5 # s
rho = 0.8
tb1 = mp1/m_dot #s
mew = 1.79e-5 # [kg/m-s]
Cp = 1006
Tw= 298 # [K]
theta1 = 89.5

# Define the functions we will be using

def m(t):
    if t <= tb1:
        return m0 - m_dot*t 
    else:
        return mpl

def D(vrx,vry):
        return 0.5*rho*(vrx**2+vry**2)*S*cd


def thetaf(t,vrx,vry):
    if t <= trail:
        return np.deg2rad(theta1)
    elif t <= 0.8406*tb1:
        return np.arctan(vry/vrx) - ((np.deg2rad(0.00325)*t)/0.01)
    else:
        return np.arctan(vry/vrx)
    
# Def the acceleration function 
        
def accel(dvrxdt,dvrydt,t,vrx,vry):
    
    if t <= trail:
        
        dvrxdt=(T/m(t)- D(vrx,vry)/m(t)-g*np.sin(thetaf(t,vrx,vry)))*np.cos(thetaf(t,vrx,vry))
        dvrydt=(T/m(t)- D(vrx,vry)/m(t)-g*np.sin(thetaf(t,vrx,vry)))*np.sin(thetaf(t,vrx,vry))    
    
    elif t <= tb1:
        
        dvrxdt = (T/m(t)-D(vrx,vry)/m(t))*np.cos(thetaf(t,vrx,vry))
        dvrydt = (T/m(t)-D(vrx,vry)/m(t))*np.sin(thetaf(t,vrx,vry)) - g
          
    else:
        
        dvrxdt=(0-D(vrx,vry)/m(t))*np.cos(thetaf(t,vrx,vry))
        dvrydt=(0-D(vrx,vry)/m(t))*np.sin(thetaf(t,vrx,vry)) - g
    
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
    dtheta[0] = theta0
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
        dtheta[i] = np.arctan(vy[i]/vx[i])
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
        
        # Normal Shock Relations and sutherlands formula
#        if Mach[i] <=1:
#            P1[i] = P[i]
#            T1[i] = Te[i]
#        else:
#        P1[i] = P[i]*((((2.4*Mach[i]**2)/((0.4*Mach[i]**2)+2))**(1.4/0.4))*((2.4/((2*1.4*Mach[i]**2)-0.4))**(1/0.4)))
#        T1[i] = Te[i]*((((2*1.4*Mach[i]**2)-0.4)*((0.4*Mach[i]**2)+2))/((2.4**2)*Mach[i]**2))
        mu[i] = mu[0]*((Te[i]/Te[0])**1.5)*((Te[0]+110.4)/(Te[i]+110.4))
        
        # Plot the Heat Fluxes
#        Tw = 298 #[K]
        Q[i] = 1.1813e-3*np.sqrt(rho[i])*v[i]**3.05
        
#        if Mach[i] < 1:
#            Q1[i] = 0
#        else:
#            Q1[i] =  (0.763*(0.744)**(-0.6))*np.sqrt((rho[i]*mu[i]))*(Cp*(T1[i]-Tw))*np.sqrt((1/(6e-3/2))*np.sqrt((2*(P1[i]-P[i]))/rho[i]))
        
        # Here we plot D(vrx,vry)/m(t)
        if t2[i] <= tb1:
            dm[i] = m0 - m_dot*t2[i]    
        else:
             dm[i] = mpl
             
        if t2[i] <= tb1:
            dthrust[i] = T/dm[i]
        else:
            dthrust[i] = 0
        
        Drag[i] = D(vx[i],vy[i])/dm[i]     
        
        
    return Q1, Q, Te, P, rho, sound, dthrust, t2, dm, dvxdt, dvydt, v, vx, vy, Mach, dx, dy, dtheta, Drag
    
         
Q1, Q, Te, P, rho, sound, dthrust, t2, dm, dvxdt, dvydt, v, vx, vy, Mach, dx, dy, dtheta, Drag = rk4(accel, 0, 0, 0, 0, 0, 0, 0, np.deg2rad(theta1), 274, 27400)
#for t, vrx, vry in list(zip(t2, vx, vy))[::10]:
    #print("%4.1f %10.5f %10.5f" % (t, vrx, vry))

dx[:] = [x/1000 for x in dx]
dy[:] = [x/1000 for x in dy]
v[:] = [x/1000 for x in v]
vx[:] = [x/1000 for x in vx]
vy[:] = [x/1000 for x in vy]
Drag[:] =[x/1000 for x in Drag]
sound[:] =[x/1000 for x in sound]
Q[:] = [x/1e6 for x in Q]
Q1[:] = [x/1e6 for x in Q1]

#fig, ax = plt.subplots(figsize=(8,8))  
#plt.plot(dx,dy,linestyle='dashed', color="blue", linewidth = 2.0)
#plt.xlabel('Range [km]',fontsize=20)
#plt.ylabel('Altitude [km]',fontsize=20)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#ax.set_xlim([0,4000])
#ax.set_ylim([0,65])
#plt.tight_layout()

##           
fig, ax = plt.subplots(figsize=(8,8))  
plt.plot(t2,dy,linestyle='dashed', color="blue", linewidth = 2.0)
plt.xlabel('t [s]',fontsize=20)
plt.ylabel('Altitude [km]',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
ax.set_xlim([0,1750])
ax.set_ylim([0,65])
plt.tight_layout()   
   
#fig, ax = plt.subplots(figsize=(8,8))  
#plt.plot(t2,dx,linestyle='dashed', color="blue", linewidth = 2.0)
#plt.xlabel('t [s]',fontsize=20)
#plt.ylabel('Range [km]',fontsize=20)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#ax.set_xlim([0,1750])
#ax.set_ylim([0,4000])
#plt.tight_layout()        
##
#fig, ax = plt.subplots(figsize=(8,8))  
#plt.plot(t2,vy,linestyle='dashed', color="blue", linewidth = 2.0)
#plt.xlabel('t [s]',fontsize=20)
#plt.ylabel('y_velocity [km/s]',fontsize=20)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#ax.set_xlim([0,1750])
#ax.set_ylim([0,1])
#plt.tight_layout()       
#
#fig, ax = plt.subplots(figsize=(8,8))  
#plt.plot(t2,vx,linestyle='dashed', color="blue", linewidth = 2.0)
#plt.xlabel('t [s]',fontsize=20)
#plt.ylabel('x_velocity [km/s]',fontsize=20)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#ax.set_xlim([0,1750])
#ax.set_ylim([0,4.5])
#plt.tight_layout()  
#
#fig, ax = plt.subplots(figsize=(8,8))  
#plt.plot(t2,Drag,linestyle='dashed', color="blue", linewidth = 2.0)
#plt.xlabel('t [s]',fontsize=20)
#plt.ylabel('Drag [kN]',fontsize=20)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#ax.set_xlim([0,1750])
#plt.tight_layout()  
#
#    
#fig, ax = plt.subplots(figsize=(8,8))  
#plt.plot(t2,np.rad2deg(dtheta),linestyle='dashed', color="blue", linewidth = 2.0)
#plt.xlabel('t [s]',fontsize=20)
#plt.ylabel('Theta [deg]',fontsize=20)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#ax.set_xlim([0,1750])
#plt.tight_layout()     
#
#fig, ax = plt.subplots(figsize=(8,8))  
#plt.plot(t2,v,linestyle='dashed', color="blue", linewidth = 2.0)
#plt.xlabel('t [s]',fontsize=20)
#plt.ylabel('Velocity Magnitude [km/s]',fontsize=20)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#ax.set_xlim([0,1750])
#ax.set_ylim([0,4.5])
#plt.tight_layout()
#
#fig, ax = plt.subplots(figsize=(8,8))  
#plt.plot(t2,Mach,linestyle='dashed', color="blue", linewidth = 2.0)
#plt.xlabel('t [s]',fontsize=20)
#plt.ylabel('Mach Number',fontsize=20)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#ax.set_xlim([0,1750])
#ax.set_ylim([0,4.5])
#plt.tight_layout()
#
#fig, ax = plt.subplots(figsize=(8,8))  
#plt.plot(dy,Te,linestyle='dashed', color="blue", linewidth = 2.0)
#plt.ylabel('Temperature [K]',fontsize=20)
#plt.xlabel('Altitude [km]',fontsize=20)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#plt.tight_layout()
#
#fig, ax = plt.subplots(figsize=(8,8))  
#plt.plot(dy,P,linestyle='dashed', color="blue", linewidth = 2.0)
#plt.ylabel('Pressure [kPa]',fontsize=20)
#plt.xlabel('Altitude [km]',fontsize=20)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#plt.tight_layout()
#
#fig, ax = plt.subplots(figsize=(8,8))  
#plt.plot(dy,rho,linestyle='dashed', color="blue", linewidth = 2.0)
#plt.ylabel('Density [kg/$m^3$]',fontsize=20)
#plt.xlabel('Altitude [km]',fontsize=20)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#plt.tight_layout()
#
#fig, ax = plt.subplots(figsize=(8,8))  
#plt.plot(dy,sound,linestyle='dashed', color="blue", linewidth = 2.0)
#plt.ylabel('Speed of Sound  [km/s]',fontsize=20)
#plt.xlabel('Altitude [km]',fontsize=20)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#plt.tight_layout()
#
#fig, ax = plt.subplots(figsize=(8,8))  
#plt.plot(t2,Q,linestyle='dashed', color="blue", linewidth = 2.0)
##ax.plot(t2,Q1,linestyle='dashed', color="red", linewidth = 2.0)
#plt.ylabel('Heat Rate [MW/$m^2$]',fontsize=20)
#plt.xlabel('time [s]',fontsize=20)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#plt.tight_layout()

#fig, ax = plt.subplots(figsize=(8,8))  
#plt.plot(v,dy,linestyle='dashed', color="blue", linewidth = 2.0)
#plt.ylabel('Altitude  [km]',fontsize=20)
#plt.xlabel('Velocity [km/s]',fontsize=20)
#plt.xticks(fontsize=20)
#plt.yticks(fontsize=20)
#plt.tight_layout()