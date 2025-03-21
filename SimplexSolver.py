class SimplexSolver:
    def __init__(self, objective_func, constraint, ineq, tolerance=1e-6):

        self.objective_func = objective_func
        self.constraints = constraint
        self.ineq = ineq
        self.num_variables = len(objective_func)
        self.num_constraints = len(constraint)
        self.tolerance = tolerance
        self.log=""
        self.col_vectors = []
        self.variable_types = []  
        self.variable_names = []  

        for i in range(self.num_variables):
            self.variable_types.append('original')
            self.variable_names.append(f'x{i+1}')

        for i in range(self.num_constraints):
            if self.ineq[i] == 1:
                col_vector = [0] * self.num_constraints
                col_vector[i] = 1
                self.col_vectors.append(col_vector)
                self.variable_types.append('slack')
                self.variable_names.append(f's{i+1}')
            
            elif self.ineq[i] == 3:
                col_vector = [0] * self.num_constraints
                col_vector[i] = -1
                self.col_vectors.append(col_vector)
                self.variable_types.append('surplus')
                self.variable_names.append(f'sp{i+1}')
                col_vector = [0] * self.num_constraints
                col_vector[i] = 1
                self.col_vectors.append(col_vector)
                self.variable_types.append('artificial')
                self.variable_names.append(f'a{i+1}')
            
            elif self.ineq[i] == 2: 
                col_vector = [0] * self.num_constraints
                col_vector[i] = 1
                self.col_vectors.append(col_vector)
                self.variable_types.append('artificial')
                self.variable_names.append(f'a{i+1}')

        self.create_simplex_table()

    def is_zero(self, value):
        return abs(value) < self.tolerance

    def is_positive(self, value):
        return value > -self.tolerance
    
    def is_negative(self, value):
        return value < +self.tolerance

    def create_simplex_table(self):

        self.total_variables = self.num_variables + len(self.col_vectors)

        self.simplex_table = []
        for i in range(self.num_constraints + 1):
            row = [0.0 for _ in range(self.total_variables + 1)]
            self.simplex_table.append(row)

        for i in range(self.num_constraints):
            for j in range(self.num_variables):
                self.simplex_table[i][j] = float(self.constraints[i][j])
            self.simplex_table[i][-1] = float(self.constraints[i][-1])
        col_idx = self.num_variables
        for col_vector in self.col_vectors:
            for i in range(self.num_constraints):
                self.simplex_table[i][col_idx] = float(col_vector[i])
            col_idx += 1

        for j in range(self.num_variables):
            self.simplex_table[-1][j] = -float(self.objective_func[j])

        col_idx = self.num_variables
        for var_type in self.variable_types[self.num_variables:]:
            if var_type == 'artificial':
                self.simplex_table[-1][col_idx] = 99999.999 
            col_idx += 1

        self.basis_variables = []
        self.non_basis_variables = list(range(self.total_variables))

        col_idx = self.num_variables
        for i in range(self.num_constraints):
            if self.ineq[i] == 1:  
                self.basis_variables.append(col_idx)
                self.non_basis_variables.remove(col_idx)
                col_idx += 1
            elif self.ineq[i] == 3:  
                col_idx += 1
                self.basis_variables.append(col_idx)
                self.non_basis_variables.remove(col_idx)
                col_idx += 1
            elif self.ineq[i] == 2: 
                self.basis_variables.append(col_idx)
                self.non_basis_variables.remove(col_idx)
                col_idx += 1
    
    def get_tableau_string(self):
        """Return a string representation of the current tableau."""
        s = ""
        header = []
        for i in range(self.total_variables):
            header.append(self.variable_names[i] if i < len(self.variable_names) else f"x{i+1}")
        header.append("RHS")
        s += "   " + " ".join([f"{h:8s}" for h in header]) + "\n"
        
        for i in range(self.num_constraints):
            basis_var = self.basis_variables[i]
            basis_name = self.variable_names[basis_var] if basis_var < len(self.variable_names) else f"x{basis_var+1}"
            row_str = f"{basis_name}: "
            for j in range(self.total_variables + 1):
                row_str += f"{self.simplex_table[i][j]:7.3f} "
            s += row_str + "\n"
        
        s += "z: " + " ".join([f"{self.simplex_table[-1][j]:7.3f}" for j in range(self.total_variables + 1)]) + "\n"
        return s
    
    def scale_row(self, row_idx, scale_factor):
        if self.is_zero(scale_factor):
            raise ValueError("Scale factor cannot be zero")
        
        for j in range(self.total_variables + 1):
            self.simplex_table[row_idx][j] *= scale_factor
            if self.is_zero(self.simplex_table[row_idx][j]):
                self.simplex_table[row_idx][j] = 0.0
    
    def row_operation(self, target_row, pivot_row, factor):
        for j in range(self.total_variables + 1):
            self.simplex_table[target_row][j] -= factor * self.simplex_table[pivot_row][j]
            # Make almost zero number zero
            if self.is_zero(self.simplex_table[target_row][j]):
                self.simplex_table[target_row][j] = 0.0
    
    def pivot(self, pivot_row, pivot_col):
        pivot_element = self.simplex_table[pivot_row][pivot_col]
        self.scale_row(pivot_row, 1.0 / pivot_element)
        for i in range(self.num_constraints + 1):
            if i != pivot_row:
                factor = self.simplex_table[i][pivot_col]
                self.row_operation(i, pivot_row, factor)
        if pivot_col in self.non_basis_variables:
            self.non_basis_variables.remove(pivot_col)
            entering_var = pivot_col
            leaving_var = self.basis_variables[pivot_row]
            self.basis_variables[pivot_row] = entering_var
            self.non_basis_variables.append(leaving_var)
            self.non_basis_variables.sort()
    
    def find_pivot_column(self):
        min_val = 0.0
        pivot_col = -1
        
        for j in self.non_basis_variables:
            if self.is_negative(self.simplex_table[-1][j]) and self.simplex_table[-1][j] < min_val:
                min_val = self.simplex_table[-1][j]
                pivot_col = j
        
        return pivot_col
    
    def find_pivot_row(self, pivot_col):
        min_ratio = float('inf')
        pivot_row = -1
        
        for i in range(self.num_constraints):
            if self.is_positive(self.simplex_table[i][pivot_col]):
                ratio = self.simplex_table[i][-1] / self.simplex_table[i][pivot_col]
                if ratio < min_ratio:
                    min_ratio = ratio
                    pivot_row = i
        
        return pivot_row
    
    def solve(self):
        max_iterations = 100 
        iteration = 0
        # Log the initial tableau.
        self.log += "Initial Tableau:\n" + self.get_tableau_string() + "\n"
        
        while iteration < max_iterations:
            pivot_col = self.find_pivot_column()

            if pivot_col == -1:
                break

            pivot_row = self.find_pivot_row(pivot_col)
            if pivot_row == -1:
                self.log += "Unbounded solution detected during iteration.\n"
                return "Unbounded solution"

            self.pivot(pivot_row, pivot_col)
            iteration += 1
            self.log += f"After iteration {iteration}:\n" + self.get_tableau_string() + "\n"
        
        # Log the final tableau.
        self.log += "Final Tableau:\n" + self.get_tableau_string() + "\n"
        
        for i, basis_var in enumerate(self.basis_variables):
            var_index = basis_var
            if var_index < len(self.variable_types) and self.variable_types[var_index] == 'artificial':
                if self.is_positive(self.simplex_table[i][-1]):
                    self.log += "Infeasible solution detected: artificial variable remains in basis.\n"
                    return "Infeasible solution"
                
        if iteration < max_iterations:
            print(self.log)
        
        solution = [0.0] * self.num_variables
        for i, var in enumerate(self.basis_variables):
            if var < self.num_variables:
                solution[var] = self.simplex_table[i][-1]
        objective_value = self.simplex_table[-1][-1]
        
        return {
            "solution": solution,
            "objective_value": objective_value,
            "tableau": self.simplex_table,
            "basis_variables": self.basis_variables,
            "non_basis_variables": self.non_basis_variables,
            "iterations": iteration
        }
    
    def print_tableau(self):
        print(self.get_tableau_string())
        
objective_function = [3, 5]


constraints = [
    [2, 1, 6], 
    [1, 2, 8]  
]

inequality_types = [1, 1]

solver = SimplexSolver(objective_function, constraints, inequality_types)

result = solver.solve()

print("Solution:")
print(f"x1 = {float(result['solution'][0])}")
print(f"x2 = {float(result['solution'][1])}")
print(f"Objective value = {float(result['objective_value'])}")


objective_function = [2, 3, 4]

constraints = [
    [1, 1, 1, 40],   
    [2, 1, -1, 10],  
    [1, 3, 1, 30]    
]

inequality_types = [1, 3, 2] 

solver = SimplexSolver(objective_function, constraints, inequality_types)

result = solver.solve()

print("Solution:")
try:
    for i, val in enumerate(result['solution']):
        if i < len(result['solution']):
            print(f"x{i+1} = {float(val)}")
    print(f"Objective value = {float(result['objective_value'])}")
except:
    print("Infeasibile Solution")
objective_function = [2, 3, 4]

constraints = [
    [1, 1, 1, 10],   
    [1, 1, 1, 15],  
    [1, 3, 1, 30]    
]

inequality_types = [2, 2, 1] 

solver = SimplexSolver(objective_function, constraints, inequality_types)

result = solver.solve()

print("Solution:")
try:
    for i, val in enumerate(result['solution']):
        if i < len(result['solution']):
            print(f"x{i+1} = {float(val)}")
    print(f"Objective value = {float(result['objective_value'])}")
except:
    print("Infeasibile Solution")
objective_function = [2, 3, 4]