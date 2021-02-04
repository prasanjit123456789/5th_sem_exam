from . import UTILITY as ut
def gauss_jordan(ar,b):
    n=len(ar)
    if(n!=len(b)):
        print("Invalid input")
        return b
    for i in range(n):
        a=ut.partial_pivot(ar,b)
        if not a:
            print("Determinant=0")
            return b
        pivot=ar[i][i]
        if type(b[0])!=list:b[i]/=pivot
        else: b[i] = [b[i][t] / pivot for t in range(n)]
        ar[i] = [ar[i][t] / pivot for t in range(n)]
        for j in range(n):
            if i==j or ar[j][i]==0:
                continue
            else:
                factor=ar[j][i]
                if type(b[0])!=list:b[j]-=b[i]*factor
                else: b[j] = [b[j][t] - b[i][t] * factor for t in range(n)]
                ar[j] = [ar[j][t] - ar[i][t] * factor for t in range(n)]
    return b
#------------------****************------------------#
def LU(A,b=None):
    row_rotations=0 #it counts no of row swaps
    n=len(A)
    i=0;s=0;t=False
    while i<n :
        if not t:s=i
        new_row=A[i][:]
        j=0
        while j<n:
            if i>j:
                sum=0
                for k in range(j):sum+=A[i][k]*A[k][j]
                A[i][j]=(A[i][j]-sum)/A[j][j]
                j+=1
            elif i==j:
                sum = 0
                for k in range(i):sum += A[i][k] * A[k][j]
                A[i][j]-= sum
                if A[i][j]!=0:
                    j+=1
                    t=False
                else: #if u_ii goes to zero
                    j=n
                    A[i]=new_row #undo the previous changes for the row
                    t=True
            else:
                sum = 0
                for k in range(i):sum += A[i][k] * A[k][j]
                A[i][j] -= sum
                j+=1
        if not t: i+=1
        else:#as the u_ii was zero so we need to swap this row with below rows as long as u_ii!=0 or row<n
            s+=1
            if s==n: return -1 #hence the determinant is zero
            ut.partial_pivot_swap(A,i,s,b)
            row_rotations+=1
    return row_rotations
#-----------------------*************************----------------------#
def lu_and_det(A,b=None): #b is not required for determinant. but b wiil be involved in partial pivoting
    p=LU(A,b) #so calling b is to do the LU decomposition once
    mult=1
    if p==-1: return 0
    for i in range(len(A)):
        mult*=A[i][i]
    return mult*(-1)**p
#---------------------------********************--------------------------#
def forward_substitution(A,b):
    n=len(A)
    for i in range(n):#the forward part
        sum=0
        for j in range(i):
            sum+=A[i][j]*b[j]
        b[i]-=sum
    return b
#------------------*******************-----------------------#
def backward_substitution(A,b):
    n=len(A)
    i = n - 1
    while (i > -1):#the backward part
        sum = 0
        for j in range(i + 1, n):
            sum += A[i][j] * b[j]
        b[i] -= sum
        b[i] /= A[i][i]
        i -= 1
    return b
#------------------------***********************-----------------------------#
def inv(matrix):
    identity=[[1 if i==j else 0 for i in range(len(matrix))] for j in range(len(matrix))]
    return gauss_jordan(matrix,identity)
#---------------------------************************------------------------#
def solve(A,b):
    value = lu_and_det(A, b)  # lu_and_det function decomposes Ar to LU and calculates the determinant
    if value == 0:
        print("The Determinant of operator is zero. Hence No unique solution is possible.")
    else:
        b = forward_substitution(A, b)
        b = backward_substitution(A, b)
    return b
#--------------------********************--------------#
def inv_lu(A):
    n = len(A)
    inv_mat = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    value = lu_and_det(A, inv_mat)
    if value == 0:
        print("The Determinant of operator is zero. Hence Inverse does not exist.")
    else:
        for i in range(n):
            column = [inv_mat[j][i] for j in range(n)]
            column = forward_substitution(A, column)
            column = backward_substitution(A, column)
            for j in range(n): inv_mat[j][i] = column[j]
    return inv_mat
