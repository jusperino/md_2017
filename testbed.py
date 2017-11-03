import math
import numpy as np

data = open("planetdata.txt")

num = 4
dim = 2
dt = 0.015
end = 485

F_old = [np.array([0 for d in range(dim)]) for n in range(num)]
F = [np.array([0 for d in range(dim)]) for n in range(num)]
x = [np.array([0 for d in range(dim)]) for n in range(num)]
v = [np.array([0 for d in range(dim)]) for n in range(num)]
names = ["" for n in range(num)]
m = [0 for n in range(num)]

def toint(ls):
    for i in range(len(ls)):
        ls[i] = float(ls[i])
    return ls
        
def read_data():
    i = 0
    for line in data:
        line = line.strip()
        line = line.split(" ")
        print(line)
        names[i] = line[0]
        m[i] = float(line[1])
        x[i] = np.array(toint(line[2:4]))
        v[i] = np.array(toint(line[4:6]))
        i+=1
    
def force(p,q):
    global F
    addforce = m[p]*m[q]*(x[q] - x[p])*(np.linalg.norm(x[q]-x[p]))**-3
    F[p] = F_old[p] + addforce
    F[q] = F_old[q] - addforce

def comp_forces():
    for i in range(num):
        for j in range(i+1, num):
            force(i,j)
                
def update_x():
    for i in range(num):
        x[i] = x[i] + dt*v[i] + F[i]*(dt**2)/(2*m[i])
        F_old[i] = F[i]

def update_v():
    for i in range(num):
        v[i] = v[i] + (F[i] + F_old[i])*dt/(2*m[i])
        
def verlet():
    t = 0
    comp_forces()
    while t<end:
        update_x()
        comp_forces()
        update_v()
        t += dt
        
read_data()
verlet()