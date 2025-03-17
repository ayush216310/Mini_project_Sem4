class SimplexTable:
    def __init__(self, matrix, column_vecs):
        self.matrix = matrix
        self.column_vecs = column_vecs
        
        if not (type(self.matrix) == list and all(type(row) == list for row in self.matrix)):
            raise ValueError("Needs to be a 2D list")

        if not all(all((type(value) == int or type(value) == float) for value in row) for row in self.matrix):
            raise ValueError("All elements in the matrix must be int or float")

        self.n=len(self.matrix[0])

        if not (type(self.column_vecs) == list and all(type(col) == list for col in self.column_vecs)):
            raise ValueError("Column vectors must be a list of lists")

        if not all(len(col) == len(self.matrix) for col in self.column_vecs):
            raise ValueError("Each column vector must have the same number of rows as the matrix")
        
        #after checking reformat to better fit tableau
        self.matrix = [row[:-1] for row in matrix]
        self.column_vecs = [col[:] for col in column_vecs]

        rhs = [row[-1] for row in matrix]
        self.column_vecs.append(rhs)

    @property
    def tableau(self):

        result = []
        for i in range(len(self.matrix)):

            result.append(self.matrix[i] + [col[i] for col in self.column_vecs])
        return result
    
    def display_tableau(self):
        for row in self.tableau:
            print("  ".join(f"{value:0.3f}" for value in row))

    def scale_row(self, row_index, scalar):

        self.matrix[row_index] = [value * scalar for value in self.matrix[row_index]]

        for col in self.column_vecs:
            col[row_index] = col[row_index] * scalar

    def subtract_row_multiple(self, target_row, source_row, factor):

        self.matrix[target_row] = [
            target - factor * source 
            for target, source in zip(self.matrix[target_row], self.matrix[source_row])
        ]

        for col in self.column_vecs:
            col[target_row] = col[target_row] - factor * col[source_row]
            
    
    def find_pivot_col(self):
        last_row = self.tableau[-1][:-1]
        min_val=min(last_row)
        if min_val<0:
            return last_row.index(min_val)
        else:
            return None
        
    def find_pivot_row(self, pivot_col):
        min_ratio = float('inf')
        pivot_row = None
        
        for i in range(len(self.tableau) - 1):
            denominator = self.tableau[i][pivot_col]
            rhs = self.tableau[i][-1]
            
            if denominator!= 0: 
                ratio = rhs / denominator
                if ratio >= 0 and ratio < min_ratio:
                    min_ratio = ratio
                    pivot_row = i

        return pivot_row
    
    def get_solution(self):
        solution = {}
        tolerance = 1e-4
        num_decision = len(self.matrix[0])
        num_constraints = len(self.tableau)
        
        for j in range(num_decision):
            count_one = 0
            count_zero = 0
            unit_row = None
            for i in range(num_constraints):
                val = self.tableau[i][j]
                if abs(val - 1) < tolerance:
                    count_one += 1
                    unit_row = i
                elif abs(val) < tolerance:
                    count_zero += 1
            is_unit_vec = (count_one == 1 and (count_one + count_zero == num_constraints))
            if is_unit_vec:
                solution[f"x{j+1}"] = self.tableau[unit_row][-1]
            else:
                solution[f"x{j+1}"] = 0

        solution["max_z"] = self.tableau[-1][-1]
        return solution
    
    def solve(self):
        while True:
            self.display_tableau()
            print("-" * 30)

            pivot_col = self.find_pivot_col() 
            if pivot_col is None:
                break

            pivot_row = self.find_pivot_row(pivot_col)
            if pivot_row is None:
                print("Unbounded solution.")
                return

            pivot_value = self.tableau[pivot_row][pivot_col]

            self.scale_row(pivot_row, 1 / pivot_value)

            for i in range(len(self.tableau)):
                if i != pivot_row:
                    factor = self.tableau[i][pivot_col]
                    self.subtract_row_multiple(i, pivot_row, factor)

        print("Final Tableau:")
        self.display_tableau()
        print("Solution:")
        sol=self.get_solution()
        for var,val in sol.items():
            print(f"{var}={val:0.3f}")
            
    

        
#Testing Code
matrix = [
    [4, 2, 7, 30],    # Constraint 1: 4x1 + 2x2 + 7x3 <= 30
    [3, 6, 2, 25],    # Constraint 2: 3x1 + 6x2 + 2x3 <= 25
    [5, 4, 3, 28],    # Constraint 3: 5x1 + 4x2 + 3x3 <= 28
    [2, 7, 4, 26],    # Constraint 4: 2x1 + 7x2 + 4x3 <= 26
    [-3, -5, -2, 0]   # Objective row: maximize z = 3x1 + 5x2 + 2x3 
]

column_vecs = [
    [1, 0, 0, 0, 0],
    [0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0] 
]
table = SimplexTable(matrix, column_vecs)

table.display_tableau()
extended_tableau = table.tableau
print("Tableau as a 2D list:")
print(extended_tableau)
pc=table.find_pivot_col()
print(pc)
pr=table.find_pivot_row(pc)
print(pr)
print(table.tableau[pr][pc])
table.solve()