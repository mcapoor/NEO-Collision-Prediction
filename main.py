import math
import requests
from bs4 import BeautifulSoup
import lxml
import re

#Creates list of NEOs from risk list 
names = []
with open("C:\\Users\\capoo\\Downloads\\esa_risk_list_20220116_1836.txt", 'r+') as file:
    for line in file.readlines(): 
        names.append(line.strip())

print(names[6])

def C(x, a):
    z = (x**2) / a

    if (z < 0):
        return (1 - math.cosh(math.sqrt(-z))) / z
    else:
        return (1 - math.cos(math.sqrt(math.radians(z))) / z)

def S(x, a):
    z = (x**2) / a

    if (z < 0):
        return (math.sinh(math.sqrt(-z)) - math.sqrt(-z)) / math.sqrt((-z)**3)
    else:
        return (math.sqrt(z) - math.sin(math.sqrt(math.radians(z)))) / (z ** 1.5)

def dt(x, a, vec_r_0, vec_v_0, r_0): 
    return (1 / math.sqrt(mu)) * (((x ** 2) * C(x, a)) + (((vec_r_0[0] * vec_v_0[0]) + (vec_r_0[1] * vec_v_0[1])) / math.sqrt(mu)) * x * (1 - (((x**2) / a) * S(x, a))) + (r_0 * (1 - (((x ** 2) / a) * C(x, a)))))
    
def nu_epoch(M, e):
    return M + ((2 * e - (0.25 * e**3)) * math.sin(math.radians(M))) + (1.25 * (e**2) * math.sin(math.radians(2 * M))) + ((13/12) * (e**3) * math.sin(math.radians(3 * M)))

def t_n(x, a, vec_r_0, vec_v_0, r_0):
    return (1 / math.sqrt(mu)) * (((x ** 2) * C(x, a)) * (((vec_r_0[0] * vec_v_0[0]) + (vec_r_0[1] * vec_v_0[1])) / math.sqrt(mu)))  + (1 - (r_0 / a) * ((x**3) * S(x, a))) + (r_0 * x) 

def propagate(a, vec_r_0, vec_v_0, r_0):
    x_0 = math.sqrt(mu) / a

    t = 0
    while (abs(1 - t) > 0.01):
        #t = t_n(x_0, a, vec_r_0, vec_v_0, r_0) 
        x_1 = x_0 + ((1 - t) / dt(x_0, a, vec_r_0, vec_v_0, r_0))
        if (x_1 > 1.4 or x_1 < 0):
            x_1 = x_0 
            break
        else:
            x_0 = x_1
        
    x = x_1 

    r =  math.sqrt(mu) * dt(x, a, vec_r_0, vec_v_0, r_0)

    f = 1 - ((x **2) / r_0) * C(x, a)
    g = t - ((x **3) / math.sqrt(mu)) * S(x, a)
    r_final = ((f * vec_r_0[0] + g * vec_v_0[0]), (f * vec_r_0[1] + g * vec_v_0[1]))

    f_dot = (math.sqrt(mu) / (r_0 * r)) * x * ((x**2 / a) * S(x, a) - 1)
    g_dot = 1 - ((x**2 / r) * C(x, a))
    v_final = ((f_dot * vec_r_0[0] + g_dot * vec_v_0[0]), (f_dot * vec_r_0[1] + g_dot * vec_v_0[1]))

    return r_final, v_final, x, r

mu = 1.327 #* (10 ** 20) 

#Creates array of Earth Positions
earth_epoch_0 = 51544 #J2000
epoch_final = 96105 #Jan 1, 2122

earth_e = 0.016709
earth_a = 1
earth_inclination = 0
earth_orbital_node = 0
earth_M = 356.0470

Earth_pos = []

for time in range(epoch_final - earth_epoch_0):
    nu = nu_epoch(earth_M, earth_e)
    r_0 = (earth_a * (1 - earth_e**2)) / (1 + earth_e * math.cos(math.radians(nu))) 

    vec_r_0 = (r_0 * math.cos(math.radians(nu)), r_0 * math.sin(math.radians(nu)))
    vec_v_0 = (math.sqrt(mu / (earth_a * (1 - (earth_e**2)))) * (-math.sin(math.radians(nu))), math.sqrt(mu / (earth_a * (1 - (earth_e**2)))) * (earth_e + math.cos(math.radians(nu))))

    vec_r, vec_v, x, r = propagate(earth_a, vec_r_0, vec_v_0, r_0)
    
    Earth_pos.append(vec_r)

    earth_M += math.sqrt(mu / (earth_a**3))
    print("Earth: ", time, "/", (epoch_final - earth_epoch_0 ))

pos_collection = []
with open('C:\\users\\capoo\\documents\\github\\internal assessments\\physics\\positions.txt', 'w+') as f:
    for index in range(len(names)):
        html_doc = requests.get("https://newton.spacedys.com/neodys/index.php?pc=1.1.1&n=" + names[index - 1])
        soup = BeautifulSoup(html_doc.content, 'lxml')

        epoch = float(re.search('[0-9]{2,}', str(soup.find('caption'))).group(0))
        a = float(soup.find(title="Semimajor axis").find_next().contents[0])
        e = float(soup.find(title="Eccentricity").find_next().contents[0])
        i = float(soup.find(title="Inclination").find_next().contents[0])
        ascending_node = float(soup.find(title="Ascending node").find_next().contents[0])
        M = float(soup.find(title="Anomaly").find_next().contents[0])
    
        NEO_dist = []
        for time in range(int(epoch_final - epoch)):
            nu = nu_epoch(M, e)
            r_0 = (a * (1 - e**2)) / (1 + e * math.cos(math.radians(nu))) 
            vec_r_0 = (r_0 * math.cos(math.radians(nu)), r_0 * math.sin(math.radians(nu)))
            vec_v_0 = (math.sqrt(mu / (a * (1 - e**2))) * (-math.sin(math.radians(nu))), math.sqrt(mu / (a * (1 - e**2))) * (e + math.cos(math.radians(nu))))
        
            vec_r, vec_v, x, r = propagate(a, vec_r_0, vec_v_0, r_0)
            r = math.sqrt(mu) * dt(x, a, vec_r, vec_v, r_0) 
            
            M += math.sqrt(mu / (a**3))
            print("NEO", index, "/", len(names), ": t =", time, "/", int(epoch_final - epoch))

            dist = math.sqrt((Earth_pos[time][0] - vec_r[0])**2 + (Earth_pos[time][1] - vec_r[1])**2)
            NEO_dist.append(dist)
        
        pos_collection.append(min(NEO_dist))
        f.write(f"{(names[index], min(NEO_dist))}\n")

near_collisions = []
for entry in pos_collection:
    if (entry < 0.05):
        near_collisions.append(entry)

print(len(near_collisions))
#print(f"{sorted(near_collisions, key = lambda x: x[1])[:10]}\n")
