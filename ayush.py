print("Start")
z = []

n = int(input("Enter number of variables : "))

#input objective function
print("Enter coeffs : ")
for i in range(0,n):
    z.append(int(input(f"Enter {i+1} coeff : ")))

#input constraints
const = []
for i in range(0,n):
    coeff_const = []
    for j in range(0,n+1):
        coeff_const.append(int(input(f"Enter coeff of {j+1} variable of {i+1} constraint : ")))
    const.append(coeff_const)

#burst variables
for i in range(0,n) :
    z.append(0)
    for j in range(0,n):
        if i == j :
            const[i].append(1)
        else :
            const[i].append(0)


#finding the key 
# z_j = []
# for i in range(0,n):
#     sum = 0
#     for j in range(0,n):
#         sum += const[i][j]*z[]