#accepting different number of variables and constraints as input
import numpy as np
import pandas as pd
import os
import time

#INPUT THE CONSTRAINTS AND OBJECTIVE FUNCTION
n = int(input("Enter number of constraints : "))
m = int(input("Enter numbers of variables again : "))
table1 = np.zeros((n+1,n+m+2))
constraint=[0 for x in range(n)]
row_lim = n+1
col_lim=n+m+2
for i in range(0,n) :
    obj_func = input("Enter constraint coefficient : ").split()
    obj_func = [int(x) for x in obj_func]
    constraint[i] = int(input("Enter 1 for greater than, 2 for less than and 3 for equals : "))
    for j in range(0,m) :
        table1[i][j] = obj_func[j]
    table1[i][(col_lim)-1] = obj_func[len(obj_func)-1]

for i in range(0,n):
    if constraint[i] != 3 :
        if constraint[i] == 1:
            table1[i][m + i] = -1
        else :
            table1[i][m + i] = 1


obj_func = input("Enter objective function coefficient : ").split()
obj_func = [int(x) for x in obj_func]
for j in range(0, m):
    table1[n][j] = -1*obj_func[j]
table1[n][col_lim-2] = 1
table1[n][col_lim-1] = 0

row1=[]
for i in range(0,n):
    row1.append(f"Constraint-{i+1}")
row1.extend(["Function"])

column1=[]
for i in range(0,m):
    column1.append(f"x{i+1}")
for i in range(0,n):
    column1.append(f"s{i+1}")
column1.extend(["z", "c"])

df = pd.DataFrame(table1, index=row1, columns=column1)
print(df)
print("\n")

#CHECKING FOR NEGATIVE VALUE IN LAST ROW
def check(arr2):
    n1 = 0
    pos = -1
    for i in range(0, col_lim-1) :
        if arr2[i]<n1:
            n1 = arr2[i]
            pos = i
    return pos

#CHECKING FOR SMALLEST RATIO
def check_ratio(arr3,pk):
    min_ratio = 1000.0
    pt = 0
    for i in range(0,n):
        if arr3[i][pk] > 0:
            ratio = arr3[i][-1] / arr3[i][pk]
            if ratio < min_ratio and ratio>0:
                min_ratio = ratio
                pt = i
    return pt

#WORKING
def simplex(table1, row_lim, n):
    while check(table1[n]) != -1 :
        pc = check(table1[n])
        pr = check_ratio(table1,pc)
        pivot = table1[pr][pc]
        table1[pr, :] /= pivot
        for i in range(0, row_lim):
            if i != pr :
                table1[i, :] -= (table1[i][pc]*table1[pr, :])
        print(df)
        print("\n")
        time.sleep(3)

simplex(table1, row_lim, n)
