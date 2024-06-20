# Solar-System
This repository is meant to showcase a personal project I tried to pursue after being introduced to numerical methods for solving ordinary differential equations. My main goal was to create a first approximative simulation of the Solar System which I will try to explain here in its functioning. Many important details are missing from the simulation, which is not meant to be complete in every aspect, but rather it is meant to be a sort of personal "playground" with which I tried to experiment something more difficult than a straight up integration of a differential equation. I decided to put up a repository for that so that anyone interested is welcomed to use, copy, modify what I've done or maybe suggest corrections or point out errors in my approach. 

## Core functioning and initial-data
The core of the integrator is the following funcion, which I quote here:
```python
def rk4(t, timestep,vars,pars,funcs):
    k1=[]
    k2=[]
    k3=[]
    k4=[]
    temp=[]
    for i in range(0,len(vars)):
        k1.append(timestep*funcs[i](t,vars,pars))

    for i in range(0,len(vars)):
        temp.append(vars[i]+k1[i]/2.)

    for i in range(0,len(vars)):
        k2.append(timestep*funcs[i](t+timestep/2.,temp,pars))

    for i in range(0,len(vars)):
        temp[i]=vars[i]+k2[i]/2.

    for i in range(0,len(vars)):
        k3.append(timestep*funcs[i](t+timestep/2.,temp,pars))

    for i in range(0,len(vars)):
        temp[i]=vars[i]+k3[i]

    for i in range(0,len(vars)):
        k4.append(timestep*funcs[i](t+timestep,temp,pars))

    for i in range(0,len(vars)):
        vars[i]+=(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.

    return vars
```
Here "t" is the current value of the independent variable of integration (which is interpreted to be time). 
"timestep" is the integration step.
"vars" is meant to be a vector containing all the values of the variables which are integrated at the current step of integration. This vector si updated throughrout the function and returned thereafter.
"pars" is meant to be a vector containing all the parameters needed for the problem (in this case, masses)
"funcs" is a vector of functions which are used in the integration process

The starting point of this dynamical system is the following vector, which encapsulates the initial data
