import numpy as np

def main():
    num_vars = int(input("Enter number of variables: "))
    num_constraints = int(input("Enter number of constraints: "))
    
    print("Enter the coefficients of the objective function separated by spaces:")
    obj_coeffs = list(map(float, input().split()))
    if len(obj_coeffs) != num_vars:
        print("Error: number of coefficients does not match number of variables.")
        return
    constraints = []  
    print("Enter each constraint on one line in the format:")
    print("a1 a2 ... an [<=,>=,=] b")
    for i in range(num_constraints):
        line = input(f"Constraint {i+1}: ")
        tokens = line.split()
        if len(tokens) != num_vars + 2:
            print("Error: Incorrect number of tokens in constraint.")
            return
        
        coeffs = list(map(float, tokens[:num_vars]))
        op = tokens[num_vars]
        rhs = float(tokens[num_vars+1])
        if op == ">=":
            coeffs = [-c for c in coeffs]
            rhs = -rhs
            add_slack = True 
        elif op == "<=":
            add_slack = True 
        elif op == "=":
            add_slack = False  
        else:
            print("Error: Operator must be <=, >=, or =.")
            return

        constraints.append((coeffs, rhs, add_slack))

    slack_count = sum(1 for (_, _, add_slack) in constraints if add_slack)

    total_cols = num_vars + slack_count + 1

    total_rows = num_constraints + 1

    tableau = np.zeros((total_rows, total_cols), dtype=float)
    
    tableau[0, :num_vars] = -np.array(obj_coeffs)

    slack_index = num_vars  
    for i, (coeffs, rhs, add_slack) in enumerate(constraints):
        row_index = i + 1  
        tableau[row_index, :num_vars] = coeffs
        if add_slack:
            tableau[row_index, slack_index] = 1.0
            slack_index += 1 
        tableau[row_index, -1] = rhs

    header = [f"x{j+1}" for j in range(num_vars)]
    header += [f"s{k+1}" for k in range(slack_count)]
    header.append("RHS")
    
    print("\nSimplex Tableau:")
    for name in header:
        print(f"{name:>8}", end=" ")
    print()
    print("-" * (9 * len(header)))

    for row in tableau:
        for value in row:
            print(f"{value:8.2f}", end=" ")
        print()

if __name__ == "__main__":
    main()
