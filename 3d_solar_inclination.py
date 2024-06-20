import matplotlib.pyplot as plt
import numpy as np
import time
import matplotlib.animation as animation
from scipy.spatial.transform import Rotation as R
import math

#https://nssdc.gsfc.nasa.gov/planetary/factsheet/
#https://ssd.jpl.nasa.gov/sats/elem/#legend

# started = time.time()

#NOTES:
#Some planets behave strangely if the integration step is to large. This is the case of the moons of mars "phobos" and "deimos"
#when, for example, the step is set to 86400 (that is, a single earth day). Due to the fact that the orbit of mars' moons last less
#than a day this large integration step may result in an inaccurate way to depict the orbit

#Gravitational constant
G=6.67*1e-11

#USEFUL PARAMETERS
number_of_bodies = 12
draw_limit = 365
repeat_integration = 10
instant=0
step=864.00
computation_time=0
animation_time=0

#initial data (coordinates and velocities) of the system of planets. Every planet begins at its perihelion and with its
#maximum velocity
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

distances = [0,mercury_dist,venus_dist,earth_dist,earth_dist+moon_earth_dist,mars_dist,mars_dist+phobos_mars_dist,mars_dist+deimos_mars_dist,
             jupiter_dist,saturn_dist,uranus_dist,neptune_dist]

mercury_speed = 58970
venus_speed = 35260
earth_speed = 30290
moon_earth_speed = 1082
mars_speed = 26500
phobos_mars_speed = 2138 #may be wrong
deimos_mars_speed = 1351.3 #may be wrong
jupiter_speed = 13720
saturn_speed = 10140
uranus_speed = 7130
neptune_speed = 5470

#The data for every planet is initially arranged as a vector of coordinates where only the first one is different from zero, 
#that is all the planets are originally allgigned on the x axis. For example mercury is [mercury_dist,0,0], venus is [venus_dist,0,0].
#the moon is [earth_dist+moon_earth_dist,0,0] because its orbital parameters are defined relatively to the Earth. A similar reasoning
#is purused for the speed of each planet, which is originally directed along the y axis, that is, tangentially to the orbit's
#perihelion: for example mercury is [0,mercury_speed,0].
#For each planet a rotation is defined to rotate the initial conditions according to the orbital parameters defining the orbit
#orientation: (Longitude of ascending node, inclination, argument of perihelion)
mercury_rotation = R.from_euler("ZYX",[48.331,-3.38,29.124],degrees=True)
mercury_rotated = list(mercury_rotation.apply([mercury_dist,0,0]))
mercury_speed_rotated = list(mercury_rotation.apply([0,mercury_speed,0]))
mercury_z = mercury_rotated[2]

venus_rotation = R.from_euler("ZYX",[76.680,-3.86,54.884],degrees=True)
venus_rotated = list(venus_rotation.apply([venus_dist,0,0]))
venus_speed_rotated = list(venus_rotation.apply([0,venus_speed,0]))
venus_z = venus_rotated[2]

earth_rotation = R.from_euler("ZYX",[-11.26064,-7.155,	114.20783],degrees=True)
earth_rotated = list(earth_rotation.apply([earth_dist,0,0]))
earth_speed_rotated = list(earth_rotation.apply([0,earth_speed,0]))
moon_rotated = list(earth_rotation.apply([earth_dist+moon_earth_dist,0,0]))
moon_speed_rotated = list(earth_rotation.apply([0,earth_speed+moon_earth_speed,0]))
earth_z=earth_rotated[2]

mars_rotation = R.from_euler("ZYX",[49.57854,-5.65,	286.5],degrees=True)
mars_rotated = list(mars_rotation.apply([mars_dist,0,0]))
mars_speed_rotated = list(mars_rotation.apply([0,mars_speed,0]))
phobos_rotated = list(mars_rotation.apply([mars_dist+phobos_mars_dist,0,0]))
phobos_speed_rotated = list(mars_rotation.apply([0,mars_speed+phobos_mars_speed,0]))
deimos_rotated = list(mars_rotation.apply([mars_dist+deimos_mars_dist,0,0]))
deimos_speed_rotated = list(mars_rotation.apply([0,mars_speed+deimos_mars_speed,0]))
mars_z=mars_rotated[2]

jupiter_rotation = R.from_euler("ZYX",[100.464,-6.09,273.867],degrees=True)
jupiter_rotated = list(jupiter_rotation.apply([jupiter_dist,0,0]))
jupiter_speed_rotated = list(jupiter_rotation.apply([0,jupiter_speed,0]))
jupiter_z=jupiter_rotated[2]

saturn_rotation = R.from_euler("ZYX",[113.665,-5.51,339.392],degrees=True)
saturn_rotated = list(saturn_rotation.apply([saturn_dist,0,0]))
saturn_speed_rotated = list(saturn_rotation.apply([0,saturn_speed,0]))
saturn_z=saturn_rotated[2]

uranus_rotation = R.from_euler("ZYX",[	74.006,-6.48,96.998857],degrees=True)
uranus_rotated = list(uranus_rotation.apply([uranus_dist,0,0]))
uranus_speed_rotated = list(uranus_rotation.apply([0,uranus_speed,0]))
uranus_z=uranus_rotated[2]

neptune_rotation = R.from_euler("ZYX",[131.783,-6.43,273.187],degrees=True)
neptune_rotated = list(neptune_rotation.apply([neptune_dist,0,0]))
neptune_speed_rotated = list(neptune_rotation.apply([0,neptune_speed,0]))
neptune_z=neptune_rotated[2]

zs = [0,mercury_z,venus_z,earth_z,earth_z,mars_z,mars_z,mars_z,jupiter_dist,saturn_dist,uranus_dist,neptune_dist]

allowed_planets = [1,2,3,5,8,9,10,11]#Index of "real planets" (excluding moons and the sun)
names = ["Sun","Mercury","Venus","Earth","Moon","Mars","Phobos","Deimos","Jupiter","Saturn","Uranus","Neptune"]
years = [0,0,0,0,0,0,0,0,0,0,0,0]#number of revolutions around the sun for each body
displacements_r = [[],[],[],[],[],[],[],[],[],[],[],[]]#displacements from the initial position for each body
delta_theta = [[],[],[],[],[],[],[],[],[],[],[],[]]#angle between the initial position and the actual position for each body
delta_r = [[],[],[],[],[],[],[],[],[],[],[],[]]
initial_positions = [[0,0,0],mercury_rotated,venus_rotated,earth_rotated,moon_rotated,mars_rotated,phobos_rotated,deimos_speed_rotated,
                     jupiter_rotated,saturn_rotated,uranus_rotated,neptune_rotated]

init_data = list(list([0,0,0])+mercury_rotated+venus_rotated+earth_rotated+moon_rotated+mars_rotated+phobos_rotated+deimos_rotated+jupiter_rotated+saturn_rotated+uranus_rotated+neptune_rotated
+list([0,0,0])+mercury_speed_rotated+venus_speed_rotated+earth_speed_rotated+moon_speed_rotated+mars_speed_rotated+phobos_speed_rotated+deimos_speed_rotated+jupiter_speed_rotated+saturn_speed_rotated+uranus_speed_rotated+neptune_speed_rotated)

#Masses of the bodies
sun_mass = 1.9891*1e30
earth_mass = 5.972168*1e24
moon_mass = 7.342*1e22
mercury_mass = 0.055*earth_mass
venus_mass = 0.815*earth_mass
mars_mass = 0.107*earth_mass
phobos_mass = 1.060*1e16
deimos_mass = 1.51*1e15
jupiter_mass = 317.8*earth_mass
saturn_mass = 95.159*earth_mass
uranus_mass = 14.536*earth_mass
neptune_mass = 17.147*earth_mass
pars = [sun_mass,mercury_mass,venus_mass,earth_mass,moon_mass,mars_mass,phobos_mass,deimos_mass,jupiter_mass,saturn_mass,uranus_mass,neptune_mass]

#FOR ANIMATION PURPOSE:
#current planet on which the plot is centered
current_body=0
#current values used to define axes limits as "position of the current planet"+-"current_limits[index]"
current_limits = [1.33*neptune_dist,1.33*neptune_dist,1.33*neptune_z]

#Slicing (coordinates only) of the initial data which is later useful 
actual_data=[[[init_data[0+3*i],init_data[1+3*i],init_data[2+3*i]] for i in range(0,number_of_bodies)]]

#Set up for the animation
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
ax.view_init(elev=10,azim=45)

#Scatter points at the bodies' positions
graph = ax.scatter([init_data[0+3*i] for i in range(0,number_of_bodies)],
                   [init_data[1+3*i] for i in range(0,number_of_bodies)],
                   [init_data[2+3*i] for i in range(0,number_of_bodies)],
                   c=["Yellow","Gray","Orange","Blue","Gray","Red","Gray","Gray","Orange","Yellow","Blue","Blue"])

#Definitiion of the event used to change view in the animation. If the right key is pressed the animation switches
#to the desired body, updating the axes limits consequently
def onclick(event):
    global current_body,current_limits
    allowed =['0','1','2','3','4','a','5','6','7','b','8','9','A','B']
    if(event.key in allowed):
        print(event.key)
        if(event.key=='0'):#Sun
            current_body = int(event.key)
            current_limits = [1.33*neptune_dist,1.33*neptune_dist,1.33*neptune_z]
        elif(event.key=='1'):#Mercury
            current_body = int(event.key)
            current_limits = [1.33*mercury_dist,1.33*mercury_dist,1.33*mercury_z]
        elif(event.key=='2'):#Venus
            current_body = int(event.key)
            current_limits = [1.33*venus_dist,1.33*venus_dist,1.33*venus_z]
        elif(event.key=='3'):#Earth
            current_body = int(event.key)
            current_limits = [1.33*earth_dist,1.33*earth_dist,1.33*earth_z]
        elif(event.key=='4'):#Moon
            current_body = int(event.key)
            current_limits = [1.33*moon_earth_dist,1.33*moon_earth_dist,1.33*earth_z]
        elif(event.key=='a'):#Earth-Moon center of mass
            current_body=-1
            current_limits = [1.33*moon_earth_dist,1.33*moon_earth_dist,1.33*earth_z]
        elif(event.key=='5'):#Mars
            current_body = int(event.key)
            current_limits = [1.33*mars_dist,1.33*mars_dist,1.33*mars_z]
        elif(event.key=='6'):#Phobos
            current_body = int(event.key)
            current_limits = [1.33*phobos_mars_dist,1.33*phobos_mars_dist,1.33*mars_z]
        elif(event.key=='7'):#Deimos
            current_body = int(event.key)
            current_limits = [1.33*deimos_mars_dist,1.33*deimos_mars_dist,1.33*mars_z]
        elif(event.key=='b'):#Phobos-Deimos_Mars center of mass
            current_body = -2
            current_limits = [1.33*deimos_mars_dist,1.33*deimos_mars_dist,1.33*mars_z]
        elif(event.key=='8'):#Jupiter
            current_body = int(event.key)
            current_limits = [1.33*jupiter_dist,1.33*jupiter_dist,1.33*jupiter_z]
        elif(event.key=='9'):#Saturn
            current_body = int(event.key)
            current_limits = [1.33*saturn_dist,1.33*saturn_dist,1.33*saturn_z]
        elif(event.key=='A'):#Uranus
            current_body = 10
            current_limits = [1.33*uranus_dist,1.33*uranus_dist,1.33*uranus_z]
        elif(event.key=='B'):#Neptune
            current_body = 11
            current_limits = [1.33*neptune_dist,1.33*neptune_dist,1.33*neptune_z]


cid = fig.canvas.mpl_connect('key_press_event', onclick)

def on_scroll(event):
    # print(event.button, event.step)
    increment = event.step
    scale = 0

    for i in range(0,2):
        if(i!=2):
            scale = increment * distances[current_body]/10
        else:
            scale = increment * zs[current_body]/10

        if(scale>current_limits[i] and event.button=="down"):
            print("Cannot scale")
        else:
            current_limits[i]+=scale

fig.canvas.mpl_connect('scroll_event', on_scroll)

#function calculating the angle between two vectors
def angle(v_1,v_2):
    v_1=v_1/np.linalg.norm(v_1)
    v_2=v_2/np.linalg.norm(v_2)
    return np.arccos(np.clip(np.dot(v_1, v_2), -1.0, 1.0))*180/math.pi

#function which returns a three dimensional vecto containing the current coordinates of a specific body
def planet_vector(index):
    return [init_data[0+3*index],init_data[1+3*index],init_data[2+3*index]]

def planet_velocity(index):
    return [init_data[0+3*(index+number_of_bodies)],init_data[1+3*(index+number_of_bodies)],init_data[2+3*(index+number_of_bodies)]]

def COM(planet_indexes):
    total_mass=0
    center_of_mass=0
    for i in planet_indexes:
        total_mass+=pars[i]
        center_of_mass+=pars[i]*np.array(planet_vector(i))
    return center_of_mass/total_mass

#Function to call in the animation later. It updates the axes when the user wants to change view. The array "init_data"
#is the one which is later updated in the animation.
def center_on_planet(planet_index,limits,string):
    if(planet_index==-1):
        # total_mass = pars[3]+pars[4]
        # center_of_mass = (pars[3]/total_mass)*np.array(planet_vector(3))+(pars[4]/total_mass)*np.array(planet_vector(4))
        center_of_mass = COM([3,4])
        ax.set(xlim3d=(center_of_mass[0]-limits[0], center_of_mass[0]+limits[0]))
        ax.set(ylim3d=(center_of_mass[1]-limits[1], center_of_mass[1]+limits[1]))
        ax.set(zlim3d=(center_of_mass[2]-limits[2], center_of_mass[2]+limits[2]))
        string+="\nEarth-Moon COM: ({:.3E},{:.3E},{:.3E})".format(center_of_mass[0]/AU,center_of_mass[1]/AU,center_of_mass[2]/AU)
    elif(planet_index==-2):
        # total_mass = pars[5]+pars[6]+pars[7]
        # center_of_mass = (pars[5]/total_mass)*np.array(planet_vector(5))+(pars[6]/total_mass)*np.array(planet_vector(6))+(pars[7]/total_mass)*np.array(planet_vector(7))
        center_of_mass=COM([5,6,7])
        ax.set(xlim3d=(center_of_mass[0]-limits[0], center_of_mass[0]+limits[0]))
        ax.set(ylim3d=(center_of_mass[1]-limits[1], center_of_mass[1]+limits[1]))
        ax.set(zlim3d=(center_of_mass[2]-limits[2], center_of_mass[2]+limits[2]))
        string+="\nMars with moons COM: ({:.3E},{:.3E},{:.3E})".format(center_of_mass[0]/AU,center_of_mass[1]/AU,center_of_mass[2]/AU)
    else:
        #update of the axes
        ax.set(xlim3d=(init_data[0+3*planet_index]-limits[0], init_data[0+3*planet_index]+limits[0]))
        ax.set(ylim3d=(init_data[1+3*planet_index]-limits[1], init_data[1+3*planet_index]+limits[1]))
        ax.set(zlim3d=(init_data[2+3*planet_index]-limits[2], init_data[2+3*planet_index]+limits[2]))

        string+="\n"+names[planet_index]
        string+="({:.3E},{:.3E},{:.3E})".format(init_data[0+3*planet_index]/AU,init_data[1+3*planet_index]/AU,init_data[2+3*planet_index]/AU)
        dist=np.linalg.norm(np.array(planet_vector(planet_index))
                            -np.array(planet_vector(0)))
        if(planet_index in allowed_planets):
            string+="\nDistance from the sun= {:.3E}".format(np.linalg.norm(delta_r[planet_index][1])/AU)
        
        if(planet_index in allowed_planets):
            string+="\nAngle= {:.3f}".format(delta_theta[planet_index][1])
            # if(len(delta_theta[planet_index])==3):
            #     #The angular speed is caclulated like Δtheta/Δt
            #     momentum = pars[planet_index]*(np.linalg.norm(planet_vector(planet_index))**2)*np.abs(delta_theta[planet_index][2]-delta_theta[planet_index][1])/float(step)
            #     string+="\nAngular momentum= {:.5E}".format(momentum)
    return string

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

plots=[]
for i in range(0,number_of_bodies):
    plot=ax.plot(init_data[0+3*i],init_data[1+3*i],init_data[2+3*i],linewidth=.5,color="gray")[0]
    plots.append(plot)

data = ax.text2D(0.05, 0.85, "t={},frame={}".format(instant,0), transform=ax.transAxes)

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

def f_vel(i,axis,masses,vars):
    planet_coords = [[vars[0+3*i],vars[1+3*i],vars[2+3*i]] for i in range(0,number_of_bodies)]

    ret = 0.0
    den = 0.0
    for j in range(0,number_of_bodies):
        if(j==i):
            continue
        den = (planet_coords[i][0] - planet_coords[j][0])**2 + (planet_coords[i][1] - planet_coords[j][1])**2 + (planet_coords[i][2] - planet_coords[j][2])**2
        den = den**1.5
        ret -= masses[j]*G*(planet_coords[i][axis] - planet_coords[j][axis])/den
    return ret

def f_coord(i,axis,masses,vars):
    planet_vels = [[vars[0+3*i],vars[1+3*i],vars[2+3*i]] for i in range(number_of_bodies,2*number_of_bodies)]
    return planet_vels[i][axis]

def total_energy(masses, vars):
    planet_coords = [np.array([vars[0+3*i],vars[1+3*i],vars[2+3*i]]) for i in range(0,number_of_bodies)]
    planet_vels = [np.array([vars[0+3*i],vars[1+3*i],vars[2+3*i]]) for i in range(number_of_bodies,2*number_of_bodies)]

    #computation of the total kinetic energy
    tot_kinetic = 0
    for i in range(0,number_of_bodies):
        tot_kinetic += .5*masses[i]*(np.linalg.norm(planet_vels[i])**2)

    #computation of the total potential energy
    tot_potential = 0
    for i in range(0,number_of_bodies):
        for j in range(0,number_of_bodies):
            if(j==i):
                continue
            tot_potential-=.5*(G*masses[i]*masses[j]/np.linalg.norm(planet_coords[i]-planet_coords[j]))
    
    return (tot_kinetic+tot_potential)

def total_angular_momentum():
    L=0
    for planet_index in range(number_of_bodies):
        if(len(delta_r[planet_index])==2):
            velocity = (delta_r[planet_index][1]-delta_r[planet_index][0])/step
            radius = delta_r[planet_index][1]
            L += np.cross(radius,pars[planet_index]*velocity)

    return L

funcs=[]

#The functions f_coord and f_vel are created as a general template for the functions used in numerical integration to "update"
#the value of the variables (which in this case are the components of the position or the componenents of the speed of the different)
#planets. These loops serve as a specialization of those functions so that the specifics axis and the specific planet are set.
#It is important to add the commands i=i and axis=axis, even though I still do not understand how those command work, otherwise
#all the stored lambdas will share the last used value in the loop ranges
for i in range(0,number_of_bodies):
    for axis in range(0,3):
        func = lambda t,vars,pars,i=i,axis=axis: f_coord(i,axis,pars,vars)
        funcs.append(func)

for i in range(0,number_of_bodies):
    for axis in range(0,3):
        func = lambda t,vars,pars,i=i,axis=axis: f_vel(i,axis,pars,vars)
        funcs.append(func)

plt.rcParams['text.usetex'] = True

def update(frame):
    global instant,init_data,computation_time,animation_time
    #the outer loop is simply a repetition of the integration more times. By integrating more steps before plotting the animation can
    #be made faster
    start = time.time()
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
    
    computation_time += time.time()-start

    start = time.time()
    graph._offsets3d = ([init_data[0+3*i] for i in range(0,number_of_bodies)],
                   [init_data[1+3*i] for i in range(0,number_of_bodies)],
                   [init_data[2+3*i] for i in range(0,number_of_bodies)])

    bodies = [[] for i in range(0,number_of_bodies)]
    #for simplicity, due to the fact that the 3 plots in the figure refer to the 3 planets, here the data of the trajectories
    #is separated into three different arrays, which are all stored into the array "bodies". 
    for i in range(0,len(actual_data)):
        for j in range(0,number_of_bodies):
            bodies[j].append(actual_data[i][j].copy())
    
    for i in range(0,number_of_bodies):
        np_arr=np.array(bodies[i])
        #Each bodies[something] is a trajectory whose data is later unpacked (I do not know if this is the right word)
        #to update the plots
        plots[i].set_data(np_arr[:len(np_arr), :2].T)
        plots[i].set_3d_properties(np_arr[:len(np_arr), 2])

    animation_time += time.time()-start

    string = "Day={:.3f}\n".format(instant/86400) +r'$E_{\text{tot}}$' + "={:.5E}\nComputation time={:.3f}s, Animation time={:.3f}s".format(
        total_energy(pars,init_data),computation_time,animation_time)
    
    string=center_on_planet(current_body,current_limits,string)
    # test = r'\nL_{\text{tot}}'
    string += "\n"+r'$\vec{L_{\text{tot}}}$' + " = {}".format(total_angular_momentum())
    string += "\n"+r'$\vert \vec{L_{\text{tot}}}\vert$' + " = {:.7E}".format(np.linalg.norm(total_angular_momentum()))
    data.set_text(string)
    # print(frame)
    return (plots,data)

        
ani = animation.FuncAnimation(fig=fig, func=update, frames=500, interval=1)
plt.show()
# ani.save(filename="test_solar_inclination.gif", writer="imagemagick")