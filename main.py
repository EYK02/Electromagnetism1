import numpy as np
import matplotlib.pyplot as plt

R = 0.155       #Radie (m)
I = 150         #Ström (A)
Ntheta = 200    #Antal segment
deltatheta = 2*np.pi/Ntheta
theta = deltatheta*np.arange(1,Ntheta+1)  #1-kolumnvektor, observera fnutten

# Avståndsvärden för z- och x-led i meter
z_dis = np.arange(-0.2, 0.21, 0.01)
x_dis = np.array([-0.12, -0.10, -0.08, -0.06, -0.04, -0.02, 0, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13])

# Bilda lista med teoretiskfältstyrka för varje position i z-led
B_lst = []
for step in z_dis:
  B = 0
  P = np.array([0,0,step])
  for k in np.arange(Ntheta):
    r = np.array([P[0]-R*np.cos(theta[k]), P[1]-R*np.sin(theta[k]), P[2]])
    ds = deltatheta*np.array([-R*np.sin(theta[k]), R*np.cos(theta[k]), 0])
    dB = I*np.cross(ds,r)/np.linalg.norm(r)**3
    B = B + dB
  B_lst.append(B[2]*1e-7*1e4) #med mu 0/4pi (G)


# Bilda lista med teoretiskfältstyrka för varje position i x-led
B_lst_x = []
for step in x_dis:
  B = 0
  P = np.array([step,0,0])
  for k in np.arange(Ntheta):
    r = np.array([P[0]-R*np.cos(theta[k]), P[1]-R*np.sin(theta[k]), P[2]])
    ds = deltatheta*np.array([-R*np.sin(theta[k]), R*np.cos(theta[k]), 0])
    dB = I*np.cross(ds,r)/np.linalg.norm(r)**3
    B = B + dB
  B_lst_x.append(B[2]*1e-7*1e4) #med mu 0/4pi (G)


# Mätdata av spänningsvärden
z_values =     np.array([5.7, 6.4, 7.0, 7.5, 8.7, 9.5, 10.6, 11.8, 13.0, 14.3, 15.8, 17.4, 19.2, 20.8, 22.5, 24.2, 25.8, 27.4, 28.6, 29.5, 30.1, 30.4, 30.3, 29.8, 29.1, 27.9, 26.6, 25.1, 23.3, 21.7, 20.0, 18.3, 16.8, 15.3, 14.0, 12.7, 11.4, 11.4, 10.4, 9.4, 8.5]) #mV
z_values =     np.divide(z_values, 5) # till G

x_values =     np.array([50.3, 41.6, 36, 32.6, 30.8, 29.8, 29.7, 29.6, 30.1, 31.8, 34.7, 39.5, 48.5, 52]) #mV
x_values =     np.divide(x_values, 5) # till G


# Figur 1
plt.figure(1)
plt.scatter(z_dis, z_values, marker="x", label="Mätdata")
plt.plot(z_dis, B_lst, color="red", label="Teoretiska värden")

plt.ylabel("Gauss (G)")
plt.xlabel("Avstånd från centrum (m)")
plt.legend()

# Figur 2
plt.figure(2)
plt.scatter(x_dis, x_values, marker="x", label="Mätdata")
plt.plot(x_dis, B_lst_x, color="red", label="Teoretiska värden")

plt.ylabel("Gauss (G)")
plt.xlabel("Avstånd från centrum (m)")
plt.legend()

plt.show()
