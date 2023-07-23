import math
import numpy as np
import requests
from bs4 import BeautifulSoup
import lxml
import re

#CLASS: Helper Functions

def true_anomaly(M, e):
    M = math.radians(M)
    if (M < math.pi):
        E = M + e/2 
    else:
        E = M - e/2

    ratio = 1
    while (abs(ratio) > 0.0001):
        ratio = (E - e*(math.sin(E)) - M)/(1 - e*math.cos(E))
        E = E - ratio
        
    return math.acos((math.cos(E) - e) / (1 - e * math.cos(E)))

def stumpff_S(x, a):
    z = (x**2) * a
    if (z > 0):
        s = (math.sqrt(z) - math.sin(math.radians(math.sqrt(z)))) / ((math.sqrt(math.radians(z)))**3)
    elif (z < 0):
        s = (math.sinh(math.radians(math.sqrt(-z))) - math.sqrt(-z)) / ((math.sqrt(-z)) ** 3)
    else:
        s = 1/6

    return s
        
def stumpff_C(x, a):
    z = (x**2) * a
    if (z > 0):
        c = (1 - math.cos(math.radians(math.sqrt(z)))) / z
    elif (z < 0):
        c = (math.cosh(math.radians(math.sqrt(-z))) - 1) / (-z)
    else: 
        c = 1/2
    
    return c

def lagrange_f(x, a, r_0):
    return 1 - (x**2) / r_0 * stumpff_C(x, a)

def lagrange_g(x, a, t):
    return t - (1 / math.sqrt(mu)) * (x**3) * stumpff_S(x, a)

def lagrange_f_dot(x, r, r_0, a): 
    z = a * (x**2)
    return math.sqrt(mu) / r / r_0 * (z * stumpff_S(x, a) - 1) * x

def lagrange_g_dot(x, a, r): 
    return 1 - (x**2) / r * stumpff_C(x, a)

def universal_var(t, r_0, v_0, a):
    x = math.sqrt(mu) * a * t
    ratio = 1
    while (abs(ratio) > 0.001):
        C = stumpff_C(x, a)
        S = stumpff_S(x, a)

        F = (r_0 * v_0 / math.sqrt(mu)) * (x**2) * C + (1 - (a * r_0)) * (x**3) * S + (r_0 * x) - math.sqrt(mu) * t 
        
        dF = (r_0 * v_0) / math.sqrt(mu) * x *(1 - (a * (x**2) * C + r_0)) + (1 - (a * r_0)) * (x**2) * C + r_0

        ratio = F / dF
        x = x - ratio
    return x

#CLASS: Main functions
def propagate(vec_r_0, vec_v_0, t, a):
    norm_vec_r_0 = np.linalg.norm(vec_r_0)
    rad_v_0 = np.dot(vec_r_0, vec_v_0) / norm_vec_r_0

    x = universal_var(t, norm_vec_r_0, rad_v_0, a)

    f = lagrange_f(x, a, norm_vec_r_0)
    g = lagrange_g(x, a, t)

    vec_r_final = np.dot(f, vec_r_0) + np.dot(g, vec_v_0)

    norm_vec_r_final = np.linalg.norm(vec_r_final)

    f_dot = lagrange_f_dot(x, norm_vec_r_final, norm_vec_r_0, a)
    g_dot = lagrange_g_dot(x, a, norm_vec_r_final)
    vec_v_final = np.dot(f_dot, vec_r_0) + np.dot(g_dot, vec_v_0)

    return vec_r_final, vec_v_final

def state_vectors(h, e, omega, i, w, nu):
    k = ((h**2) / mu) * ((1 / (1 + e * math.cos(math.radians(nu)))))
    perifocal_r = (k * math.cos(math.radians(nu)) * np.array([[1], [0], [0]])) + (math.sin(math.radians(nu)) * np.array([[0], [1], [0]]))

    perifocal_v = (mu / h) * (-math.sin(math.radians(nu)) * np.array([[1], [0], [0]])) + ((e + math.cos(math.radians(nu))) * np.array([[0], [1], [0]]))

    z_rotation_omega = np.array([[math.cos(math.radians(omega)), math.sin(math.radians(omega)), 0], [-math.sin(math.radians(omega)), math.cos(math.radians(omega)), 0], [0, 0, 1]])

    x_rotation_i = np.array([[1, 0, 0], [0, math.cos(math.radians(i)), math.sin(math.radians(i))], [0, -math.sin(math.radians(i)), math.cos(math.radians(i))]])

    z_rotation_w = np.array([[math.cos(math.radians(w)), math.sin(math.radians(w)), 0], [-math.sin(math.radians(w)), math.cos(math.radians(w)), 0], [0, 0, 1]])

    transform_mat = z_rotation_omega.T * x_rotation_i.T * z_rotation_w

    r_final = np.dot(transform_mat, perifocal_r)
    v_final = np.dot(transform_mat, perifocal_v)
    print("State vectors calculated")
    return r_final.reshape((3,)), v_final.reshape((3,))

# ----------------------------------------------------------- #

mu = 1.327

earth_epoch_0 = 51544 #J2000
epoch_final = 96105 #Jan 1, 2122

earth_e = 0.01671022
earth_a = 1 / 1.00000011
earth_i = 0.00005
earth_omega = -11.26064
earth_M_0 = 356.0470
earth_w = 102.94719
earth_h = math.sqrt(mu * (earth_a * (1 - (earth_e **2))))

Earth_pos = []

nu = true_anomaly(earth_M_0, earth_e)
vec_r_0, vec_v_0 = state_vectors(earth_h, earth_e, earth_omega, earth_i, earth_w, nu)
for time in range(epoch_final - earth_epoch_0):
    vec_r, vec_v = propagate(vec_r_0, vec_v_0, 1, earth_a)
    vec_r = np.array(vec_r)
    vec_v = np.array(vec_v)

    Earth_pos.append(vec_r)

    vec_r_0, vec_v_0 = vec_r, vec_v
    print("Earth: ", time, "/", (epoch_final - earth_epoch_0 ))


names = []
with open("C:\\Users\\capoo\\Downloads\\esa_risk_list_20220116_1836.txt", 'r+') as file:
    for line in file.readlines(): 
        names.append(line.strip())

pos_collection = []
for name in range(len(names)):
    html_doc = requests.get("https://newton.spacedys.com/neodys/index.php?pc=1.1.1&n=" + names[name - 1])
    soup = BeautifulSoup(html_doc.content, 'lxml')

    epoch = float(re.search('[0-9]{2,}', str(soup.find('caption'))).group(0))
    a = 1 / float(soup.find(title="Semimajor axis").find_next().contents[0])
    e = float(soup.find(title="Eccentricity").find_next().contents[0])
    i = float(soup.find(title="Inclination").find_next().contents[0])
    omega = float(soup.find(title="Ascending node").find_next().contents[0])
    w = float(soup.find(title="Argument of perihelion").find_next().contents[0])
    M_0 = float(soup.find(title="Anomaly").find_next().contents[0])
    h = math.sqrt(mu * (a * (1 - (e **2))))
   
    NEO_pos = []

    nu = true_anomaly(M_0, e)
    vec_r_0, vec_v_0 = state_vectors(h, e, omega, i, w, nu)
    for time in range(int(epoch_final - epoch)):
        vec_r, vec_v = propagate(vec_r_0, vec_v_0, 1, a)
        NEO_pos.append(vec_r)
        vec_r_0, vec_v_0 = vec_r, vec_v
        print("NEO", name, "/", len(names), ": t =", time, "/", int(epoch_final - epoch))

    pos_collection.append(NEO_pos)
    
near_collisions = []
for obj in pos_collection:
    min_dist = 0 
    for pos in range(len(obj)):
        dist = abs(obj[pos] - Earth_pos[pos])
        if (dist < min_dist):
            min_dist = dist
            print(obj, min_dist, pos)
        if (dist < 0.3):
            near_collisions.append((obj, pos))
            
print(near_collisions)






