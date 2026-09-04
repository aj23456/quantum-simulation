#finte square wall
#importing library like matplotlib and njumpy
import numpy as np
import matplotlib.pyplot as plt
v_x=float(input("enter the potential of the well"))
n=int(input("number of state" ))
width=float(input("eneter width of finite square welleg(10e-9)"))
me=9.1e-31
#calculating dx for the aproximation or for matrix method and hcut for calculation
dx=width/(n+1)
hcut=1.0545718e-34
vmatrix=np.full(n,v_x)
well_end=int(n*0.6)
well_strat=int(n*0.4)
#creating a matrix for the  potential of the well and setting the value of the potential zero  in the middle of the well

vmatrix[well_strat:well_end]=0
well_width=2e-9
main_matrix= np.ones(n)
left_matrix=np.ones(n-1)
#creating eiganvector so we can get the eigan value for the energy
H=(2*np.diag(main_matrix)+np.diag(-1*left_matrix,1)+np.diag(-1*left_matrix,-1))*(hcut**2/(2*me*(dx**2)))+np.diag(vmatrix)
eiganvalue,eiganvector=np.linalg.eigh(H)
n_value=np.arange(n)*dx
num_state=min(5,n)
#original value of the eiganvalue for infinite square well is given by the formula e_n=(n^2)*(hcut^2)*(pi^2)/(2*me*(width^2))
for i in range(num_state):
    eiganvector[:,i] /= np.sqrt(np.sum(np.abs(eiganvector[:,i])**2)*dx)
    plt.plot(n_value,eiganvector[:,i],label=f"n_value={i+1}")
#comarision of the eiganvalue from the matrix method and the original value from the formula
plt.xlabel("position x(m)")
plt.ylabel("wave function")
plt.title("first 5 wave function")
plt.grid(True,alpha=0.3 )
plt.legend()
plt.show()
#ploting probablity density 
for i in range(num_state):
    plt.plot(n_value,np.abs(eiganvector[:,i])**2,label=f"n_value={i+1}")

plt.xlabel("position x(m)")
plt.ylabel("eiganvalue probablity density")
plt.title("eiganvalue density")
plt.grid(True,alpha=0.3 )
plt.legend()
plt.show()
# number of bound state
num=0
i=0
list=[]
while i<n:
    if eiganvalue[i]<v_x:
        list.append(eiganvalue[i])
        num+=i
    i+=1
print(f"{num} out of {n}")
print(f"eigannvalue of bound state={list}")