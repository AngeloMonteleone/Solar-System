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
```python
init_data = list(list([0,0,0])+mercury_rotated+venus_rotated+earth_rotated+moon_rotated+mars_rotated+phobos_rotated+deimos_rotated+jupiter_rotated+saturn_rotated+uranus_rotated+neptune_rotated
+list([0,0,0])+mercury_speed_rotated+venus_speed_rotated+earth_speed_rotated+moon_speed_rotated+mars_speed_rotated+phobos_speed_rotated+deimos_speed_rotated+jupiter_speed_rotated+saturn_speed_rotated+uranus_speed_rotated+neptune_speed_rotated)
```
The vector contains first the positions of all the planets and then all the velocities. Each name (like "mercury_rotated", "venus_rotated", ...) is a three components vector which contains the x,y and z coordinates of each planet for the positions and the x,y and z components of the velocity of each planet. Here an example follows to show how they are calculated before intialinzing the big vector:

```python
mercury_rotation = R.from_euler("ZYX",[48.331,-3.38,29.124],degrees=True)
mercury_rotated = list(mercury_rotation.apply([mercury_dist,0,0]))
mercury_speed_rotated = list(mercury_rotation.apply([0,mercury_speed,0]))
mercury_z = mercury_rotated[2]
```
```python
earth_rotation = R.from_euler("ZYX",[-11.26064,-7.155,	114.20783],degrees=True)
earth_rotated = list(earth_rotation.apply([earth_dist,0,0]))
earth_speed_rotated = list(earth_rotation.apply([0,earth_speed,0]))
moon_rotated = list(earth_rotation.apply([earth_dist+moon_earth_dist,0,0]))
moon_speed_rotated = list(earth_rotation.apply([0,earth_speed+moon_earth_speed,0]))
earth_z=earth_rotated[2]
```
This is the example for mercury and for the Earth abd its moon: the data for every planet is initially arranged as a vector of coordinates where only the first one is different from zero, that is all the planets are originally allgigned on the x axis. For example mercury is [mercury_dist,0,0], venus is [venus_dist,0,0]. The moon is [earth_dist+moon_earth_dist,0,0] because its orbital parameters are defined relatively to the Earth. A similar reasoning is purused for the speed of each planet, which is originally directed along the y axis, that is, tangentially to the orbit's perihelion: for example mercury is [0,mercury_speed,0]. 
For each planet a rotation is defined to rotate the initial conditions according to the orbital parameters defining the orbit orientation: (Longitude of ascending node, inclination, argument of perihelion)

The proper constants, relating to the intial distance and velocities are initialized at the beginning of the programme, together whith other useful parameters.
```python
AU = 149597870700
mercury_dist = 0.307499*AU
venus_dist = 0.718440*AU
earth_dist = 147098450000
moon_earth_dist = 362600000
mars_dist = 1.3814*AU
phobos_mars_dist = 9234420
deimos_mars_dist = 23455500
jupiter_dist = 4.9506*AU
saturn_dist = 9.0412*AU
uranus_dist = 18.2861*AU
neptune_dist = 29.81*AU

mercury_speed = 58970
venus_speed = 35260
earth_speed = 30290
moon_earth_speed = 1082
mars_speed = 26500
phobos_mars_speed = 2138
deimos_mars_speed = 1351.3
jupiter_speed = 13720
saturn_speed = 10140
uranus_speed = 7130
neptune_speed = 5470
```

## Live integration implementation
To implement the integration with Runge-Kutta 4 alongside the animation process the matplotlib package was employed. The last two lines of the program in fact read:
```python
ani = animation.FuncAnimation(fig=fig, func=update, frames=500, interval=1)
plt.show()
```
the function update, which is called here implements the integration method in the following loop (defined inside the function):
```python
for j in range(0,repeat_integration):
        #here the data regarding positions and velocities is updated
        init_data=rk4(instant,step,init_data,pars,funcs)
        #here the data about the positions is stored in a different array, which keeps tracks of the trajectories, that is the
        #subsequent positions of each planet at each instant of time
        actual_data.append([[init_data[0+3*i],init_data[1+3*i],init_data[2+3*i]].copy() for i in range(0,number_of_bodies)])

        #if the number of points in the trajectory is above a certain limit the first one is removed
        if(len(actual_data)>draw_limit):
            actual_data.remove(actual_data[0])
        instant=instant+step
        study_orbits()
```
The first line inside the loop updates the big vector containing all the data and then other calculations are performed, such as the ones in the "study_orbits" function. It provides calculations for the "angular position" of the planet (calculated through a scalar product between the actual position and the initial one) and the displacements which are used to calculate when an orbit is "complete": when the displacement between the actual position and the initial one reacheas a minimum an orbit is thought as "completed". Although it is a bit rough method it provides good result, as shown in the following section.
The program can also "focus" on different things when different buttons are pressed:
+ 0: Sun
+ 1: Mercury
+ 2: Venus
+ 3: Earth
+ 4: Moon
+ 5: Mars
+ 6: Phobos
+ 7: Deimos
+ 8: Jupiter
+ 9: Saturn
+ A: Neptune
+ B: Uranus
+ a: center of mass of the Earth and the Moon
+ b: center of mass of Mars and its moons

At last, before the return statement a string is constructed in order to report on the animation some data, such as the total energy, total angular momentum, the planet position...

```python
string = "Day={:.3f}\n".format(instant/86400) +r'$E_{\text{tot}}$' + "={:.5E}\nComputation time={:.3f}s, Animation time={:.3f}s".format(
        total_energy(pars,init_data),computation_time,animation_time)
    
    string=center_on_planet(current_body,current_limits,string)
    string += "\n"+r'$L_{\text{tot}}$' + " = {}".format(total_angular_momentum())
    string += "\n"+r'$\vert L_{\text{tot}}\vert$' + " = {:.7E}".format(np.linalg.norm(total_angular_momentum()))
    data.set_text(string)
    planet_name.set_text(additional_text)
```
## Some interesting results
This is an example where the animation switches between all possible bodies:

![Example with the center of mass of the Earth and the Moon](images/all_planets.gif)

The code also computes when a planet has completed an orbit. This is implemented by looking at the angle of the planet between the initial position vector and the current position vector. The employed function to compute the angle returns a value between 0° and 180°, where the 0° means that the planet is at the beginning of the orbit and 180° that it is at the end. Obviously after the first orbit the planets will not occupy the exact point in space where they started, due to numerical errors which make it so that the orbits are not perfectly closed:

```python
#function which collects useful information about the orbits, to be later used in other functions
def study_orbits():
    for i in allowed_planets:
        #Save the current displacement for the i_th planet
        displacements_r[i].append(np.linalg.norm(np.array(initial_positions[i])-np.array(planet_vector(i))))
        #Save the current angle for the i_th planet
        delta_theta[i].append(angle(np.array(planet_vector(i))-np.array(planet_vector(0)),initial_positions[i]))
        #distances from the sun
        delta_r[i].append(np.array(planet_vector(i))-np.array(planet_vector(0)))
        if(len(delta_theta[i])==4):
            #only three subsequent values are stored for the angle
            delta_theta[i].remove(delta_theta[i][0])
        if(len(delta_r[i])==3):
            #only two subsequent values are stored for the angle
            delta_r[i].remove(delta_r[i][0])
        if(len(displacements_r[i])==4):
            #only three subsequent values are stored for the displacement
            displacements_r[i].remove(displacements_r[i][0])
        #These two condition seem to work in the same way, that is searching for a minimum of the angle or for a minimum of the displacement
        # if(len(displacements_r[i])==3 and displacements_r[i][1]<displacements_r[i][2] and displacements_r[i][0]>displacements_r[i][1]):
        if(len(delta_theta[i])==3 and delta_theta[i][1]<delta_theta[i][2] and delta_theta[i][0]>delta_theta[i][1]):
            years[i]+=1
            print("Minimum dispacement: {}".format(displacements_r[i][1]))
            print(names[i]+"Period #{}: {}".format(years[i],instant/86400))
```

Notwithstanding this, the code produced results which are coherent with the actual ones. For example in a test run with an integration step of 8640 seconds (1/10 of a day) the results for the duration of the first orbit were the following. All data is reported in days (in parenthesis the real value):
+ Mercury: 88.1 (87.96)
+ Venus: 225.0 (224.7)
+ Earth: 366.2 (365.25)
+ Mars: 687.7 (686.96)
Test with integration step of 86400 seconds (1 day):
+ Jupiter: 4342 (4332)
+ Saturn: 10625 (10756)
