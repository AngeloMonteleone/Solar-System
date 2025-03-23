import numpy as np
import matplotlib.pyplot as plt
import os

allowed_planets = [1,2,3,5,8,9,10,11]#Index of "real planets" (excluding moons and the sun)
names = ["Sun","Mercury","Venus","Earth","Moon","Mars","Phobos","Deimos","Jupiter","Saturn","Uranus","Neptune"]
true_periods = [0,87.9691,224.701,365.256363004,0,686.980,0,0,4332.59,10755.70,30688.5,60195]

xs = []
ys = []

for i in allowed_planets:
    if(os.path.isfile("output_" + names[i])):
        with open("output_" + names[i],"rb") as f:
            data = np.load(f)
            print(data)
            x = []
            y = []
            for point in data:
                x.append(point[0])
                y.append(abs(point[1]-true_periods[i]))

            xs.append(x)
            ys.append(y)
            plt.scatter(x,y,s=10)
            plt.plot(x,y,label=names[i])
            print(names[i])
            
plt.legend()
plt.xlabel("Elapsed days")
plt.ylabel("Absolute deviation from true period (days)")
plt.show()