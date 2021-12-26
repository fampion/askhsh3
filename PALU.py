import numpy as np

#pa_lu δέχεται εναν πινακα A και επιστρέφει τους πίνακες P,L και U
def pa_lu(A):
    n = len(A)  #n= οι  διαστάσεις του πίνακα
    P = np.identity(n)  #Αρχικα ο P είναι ο μοναδιαίος nxn
    L = np.identity(n)  #Αρχικα και ο L είναι ο μοναδιαίος nxn
    for y in range(n):  #Οσες επαναλήψεις όσες και οι στηλες που έχει ο πίνακας,δηλαδη n,κάθε επανάληψη μηδενίζει τα στοιχεία
        #κάτω απο τον οδηγό ολης της στήλης
        grammh_odhgou = y  #Εστω οτι ο μεγαλύτερος οδηγός είναι ο πρώτος,δηλαδή το στοιχείο [y][y]
        for x in range(y + 1,n): #Ευρεση μεγαλύτερου οδηγού της y στηλης
            if abs(A[x][y]) > abs(A[grammh_odhgou][y]):
                grammh_odhgou = x
        if A[grammh_odhgou][y] == 0: #Αν ο μεγαλύτερος οδηγός είναι 0,τελειώσαμε με την στηλη,παμε στην επόμενη
            continue
        if grammh_odhgou != y : #Αν ο οδηγός δεν είναι στην διαγώνιο,αλλάζουμε τις 2 γραμμές,μετά τον οδήγο
            #πριν ειναι 0 ολα ετσι και αλλιως.
            for i in range(y,n,1):
                (A[y][i],A[grammh_odhgou][i])=(A[grammh_odhgou][i],A[y][i])
            #Πρέπει να αποθηκέυσω και στον πίνακα P οτι έκανα μια αλλαγή γραμμών
            P[[y,grammh_odhgou]]=P[[grammh_odhgou,y]]
            #Πρεπει τις αλλαγές μου να τις αποθηκέυσω και στον πίνακα L
            for i in range(y):
                (L[y][i],L[grammh_odhgou][i]) = (L[grammh_odhgou][i],L[y][i])
        #Τωρα ξερουμε την γραμμη του οδηγου, παμε να μηδενισουμε ολα τα απο κατω
        for x in range(y+1,n,1): #Για καθε στοιχείο κάτω απο τον οδηγο:
            timh=A[x][y] #Παιρνουμε την τιμη του

            if timh==0: #Αν ειναι ηδη μηδενισμενο, παμε στο απο κάτω
                continue
            odhgos=A[y][y]
            c=timh/odhgos
            #Μηδενίζουμε τον οδηγό, και αφαιρούμε απτα στοιχεια που εναι κάτω του το c*την τιμη στοιχειου  γραμμης οδηγου
            A[x][y]=0
            for i in range(y+1,n,1):
                A[x][i] -= float(c)*A[y][i]
            #Γραφουμε και το c στον πίνακα L
            L[x][y]=c
    return (P,L,A)

#Δέχεται τον L και τον b και επιστρέφει τον y
def forward_subtitution(L,b):
    n=len(b) #παίρνω μέγεθος του συστήματος
    x=np.empty(n) #παίρνω μέγεθος του συστήματος
    for i in range(n):
        if L[i][i] == 0:  #Αν η διαγώνιος ειναι 0 τότε
            x[i]=0
            continue
        #Υπολογισμος της τιμής της i-οστης μεταβλητής:
        timh=b[i]
        for j in range(i):
            timh -= L[i][j] * x[j]
        #Διαιρώ με τον συντελεστή απο τον πίνακα L
        timh /=L[i][i]
        #αποθηκέυω την τιμή στο διάνυσμα x
        x[i]=timh
    return x

#Δέχεται τον U και τον y και επιστρτρέφει το διάνυσμα x, σχεδόν ίδια με forward, μονο που αντι να πάει απο πάνω -> κάτω
#πάει απο κάτω -> πανω
def backward_subtitution(U,b):
    n=len(b)
    x=np.empty(n)
    for i in range(n-1,-1,-1):
        if U[i][i] == 0:
            x[i]=0
            continue
        timh=b[i]
        for j in range(i+1,n,1):
            timh-=U[i][j] * x[j]
        timh/=U[i][i]
        x[i]=timh
    return x

def system_solver_with_pa_lu(A,b):
    print("Πίνακας Α:\n"+str(A))
    print("Διάνυσμα b:\n"+ str(b))
    (P,L,U) = pa_lu(A)
    b = np.matmul(P,b)#Πολλαπλασιάζω και τον b με τον P που βρήκα απτην αποσύνθεση PA=LU
    #print("Λύνω το Ly=Pb ως προς y, μέσω της forward_subtitution")
    print("P=\n"+str(P))
    print("L=\n"+str(L))
    print("U=\n"+str(U))
    y=forward_subtitution(L,b)
    #print("μετα λύνω την y=ux ως προς x μέσω της backward_subtitution")
    x=backward_subtitution(U,y)
    return x


#Σαν input A θα βάλω τον παρακάτω,για τον οποίον είχαμε βρεί τους P,L και U σαν παράδειγμα στο μάθημα
A = np.array([[2,1,5],[4,4,-4],[1,3,1]])
#A = np.array([[-1,1,6],[-4,-8,6],[2,16,23]])
#A = np.array([[1,2,-1],[1,-1,2],[2,2,2]])
#Σαν input για τον b θα βάλω το διάνυσμα αυτο
b = [3,4,5]
result=system_solver_with_pa_lu(A,b)
print("Το δίανυσμα x είναι το :"+str(result))


