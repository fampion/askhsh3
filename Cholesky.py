from math import sqrt
import numpy as np
def cholesky(A):
    n=len(A)#Παιρνω την διασταση του πινακα A, nxn
    L=np.zeros((n,n)) #Δημιουργώ πίνακα nxn με μηδενικά
    for i in range(n):
        for j in range(n):
            sum=0
            if i==j: #Αθροίσματα για διαγώνιο
                for k in range(j):
                    sum+=pow(L[j,k],2)
                L[j][j]=int(sqrt(A[j][j]-sum))
            else: #Αθροίσματα για υπόλοιπα στοιχεία, δηλαδή για i>j
                for k in range(j):
                    sum+=(L[i][k]*L[j][k])
                if(L[j][j]>0):
                    L[i][j]=int(( A[i][j]-sum) /L[j][j])

    return L

A=np.array([[4,12,-16],[12,37,-43],[-16,-43,98]])
print("A=\n"+str(A)+"\n")
L=cholesky(A)
print("L=\n"+str(L))