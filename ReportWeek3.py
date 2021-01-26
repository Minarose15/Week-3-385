# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 13:02:52 2021

@author: Sabrina
"""

from numpy import sin, array, arange, exp, pi
from pylab import plot, show, axes#, xkcd

#This function returns the derivative of each element in the input array. 
def derivs(vars,Fd,t):
    theta = vars[0] #separates out position components
    omega = vars[1] #separates out velocity
    g = 9.8 #gravity
    l = 9.8 #length of pendulum in meters. 
    q = 0.5
    #q = 0.5 #drag coefficient, takes mass into account

    thetaDeriv = omega #set derivative
    omegaDeriv = -g/l*sin(theta)-q*omega+Fd*sin((2/3)*t) #set derivative
    return array([thetaDeriv,omegaDeriv]) #returns values

#This gets the trajectory using second-order Runge-Kutta
def RK2Trajectory(theta0,omega,dt,Fd):
    theta = [theta0]
    time = [0]
    vars = array([theta0,omega]) #defines array with all needed values
    while time[-1] <= 150:
        k1 = dt * derivs(vars,Fd,time[-1]) #call derivative function to complete euler's steps
        k2 = dt * derivs(vars + 1/2 * k1,Fd,time[-1]) #find values at half-step
        vars += k2 #complete runge_kutta step
        theta.append(vars[0])
        time.append(time[-1] + dt)

    return theta, time

dt = 0.04
thetai = 0.2
omega0 = 0.0

rangeRKa, timea = array(RK2Trajectory(thetai,omega0,dt,1.2)) 
rangeRKb, timeb = array(RK2Trajectory(thetai+0.001,omega0,dt,1.2))

#this finds the fit with the lyupanov coefficient lambda, guesswork
lamb = 0.15
timec = arange(0,82,0.04)
dtheta = 0.0001*exp(lamb*timec)

deltaTheta = abs(rangeRKa-rangeRKb)
        
#xkcd()
plot(timea, deltaTheta, 'b')
plot(timec, dtheta, "lavenderblush", linestyle = '--')
axes().set_yscale("log")
axes().set_facecolor('k')
#axes().set_aspect(20)
show() 
