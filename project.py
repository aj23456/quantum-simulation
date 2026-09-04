#infinte square wall
#importing library like matplotlib and njumpy
import numpy as np
import matplotlib.pyplot as plt
n=int(input("number of state" ))
width=float(input("eneter width of infinite square well"))
me=9.1e-31
#calculating dx for the aproximation or for matrix method and hcut for calculation
dx=width/(n+1)
hcut=6.626e-34/(2*np.pi)
main_matrix= np.ones(n)
left_matrix=np.ones(n-1)
#creating eiganvector so we can get the eigan value for the energy
H=(2*np.diag(main_matrix)+np.diag(-1*left_matrix,1)+np.diag(-1*left_matrix,-1))*(hcut**2/(2*me*(dx**2)))
eiganvalue,eiganvector=np.linalg.eigh(H)

n_value=np.arange(1,n+1)
#original value of the eiganvalue for infinite square well is given by the formula e_n=(n^2)*(hcut^2)*(pi^2)/(2*me*(width^2))
eiganlvalue_real=(n_value**2)*(hcut**2)*(np.pi**2)/(2*me*(width**2))
#comarision of the eiganvalue from the matrix method and the original value from the formula
plt.plot(n_value,eiganvalue,"r--",label="esiteted")
plt.plot(n_value,eiganlvalue_real,"b",label="calculated")
plt.xlabel("n")
plt.ylabel("eiganvalue")
plt.title("comparision of the eiganvalue")
plt.grid(True,alpha=0.3 )
plt.show()
