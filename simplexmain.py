import tkinter as tk
from tkinter import ttk, messagebox

class SimplexSolver:
    def __init__(self, objective_func, constraint, ineq, is_maximization=True, tolerance=1e-6): 
        self.objective_func = objective_func
        self.constraints = constraint
        self.ineq = ineq
        self.is_maximization = is_maximization  
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
        return value > self.tolerance

    def is_negative(self, value):
        return value < -self.tolerance

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

        if not self.is_maximization:
            self.objective_func = [-x for x in self.objective_func]

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

        M = 99999.999  
        for i in range(self.num_constraints):
            var_idx = self.basis_variables[i]
            if var_idx < len(self.variable_types) and self.variable_types[var_idx] == 'artificial':
                for j in range(self.total_variables + 1):
                    self.simplex_table[-1][j] -= M * self.simplex_table[i][j]

    def get_tableau_string(self):
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

        self.log += "Final Tableau:\n" + self.get_tableau_string() + "\n"

        for i, basis_var in enumerate(self.basis_variables):
            var_index = basis_var
            if var_index < len(self.variable_types) and self.variable_types[var_index] == 'artificial':
                if self.is_positive(self.simplex_table[i][-1]):
                    self.log += "Infeasible solution detected: artificial variable remains in basis.\n"
                    return "Infeasible solution"

        solution = [0.0] * self.num_variables
        for i, var in enumerate(self.basis_variables):
            if var < self.num_variables:
                solution[var] = self.simplex_table[i][-1]
        objective_value = self.simplex_table[-1][-1]

        if not self.is_maximization:
            objective_value = -objective_value

        return {
            "solution": solution,
            "objective_value": objective_value,
            "tableau": self.simplex_table,
            "basis_variables": self.basis_variables,
            "non_basis_variables": self.non_basis_variables,
            "iterations": iteration
        }

class SimplexSolverGUI:
    def __init__(self, master):
        self.master = master
        master.title("Simplex Solver - Simplex Method")

        self.num_variables = 3
        self.num_constraints = 3

        self.objective_entries = []
        self.constraint_entries = []
        self.constraint_ineq_vars = []
        self.problem_type = tk.StringVar(value="Maximize") 
        self.create_problem_type_selector() 

        tk.Label(master, text="Number of variables:").pack()
        self.var_entry = tk.Entry(master)
        self.var_entry.pack()
        self.var_entry.insert(0, str(self.num_variables))
        tk.Button(master, text="Set Variables", command=self.modify_variables).pack()

        tk.Label(master, text="Number of constraints:").pack()
        self.con_entry = tk.Entry(master)
        self.con_entry.pack()
        self.con_entry.insert(0, str(self.num_constraints))
        tk.Button(master, text="set constraints", command=self.modify_constraints).pack()

        tk.Label(master, text="Objective Function:").pack()
        self.objective_frame = tk.Frame(master)
        self.objective_frame.pack()
        self.create_objective_entries()

        tk.Label(master, text="constraints:").pack()
        self.constraints_frame = tk.Frame(master)
        self.constraints_frame.pack()
        self.create_constraint_entries()

        self.Solve_button = tk.Button(master, text="Solve", command=self.run_simplex)
        self.Solve_button.pack()

        self.result_text = tk.Text(master, height=25, width=80)
        self.result_text.pack()

    def create_problem_type_selector(self):
        frame = tk.Frame(self.master)
        frame.pack(pady=5)

        maximize_radio = tk.Radiobutton(frame, text="Maximize", variable=self.problem_type, value="Maximize")
        minimize_radio = tk.Radiobutton(frame, text="Minimize", variable=self.problem_type, value="Minimize")

        maximize_radio.pack(side=tk.LEFT, padx=5)
        minimize_radio.pack(side=tk.LEFT, padx=5)

    def create_objective_entries(self):
        for widget in self.objective_frame.winfo_children():
            widget.destroy()
        self.objective_entries = []
        for i in range(self.num_variables):
            entry = tk.Entry(self.objective_frame, width=5)
            entry.pack(side=tk.LEFT)
            entry.insert(0, "1")
            self.objective_entries.append(entry)
            if i < self.num_variables - 1:
                tk.Label(self.objective_frame, text=f"x{i+1} +").pack(side=tk.LEFT)
            else:
                tk.Label(self.objective_frame, text=f"x{i+1}").pack(side=tk.LEFT)

    def create_constraint_entries(self):
        for widget in self.constraints_frame.winfo_children():
            widget.destroy()
        self.constraint_entries = []
        self.constraint_ineq_vars = []
        for i in range(self.num_constraints):
            frame = tk.Frame(self.constraints_frame)
            frame.pack(pady=2)
            row_entries = []
            for j in range(self.num_variables):
                entry = tk.Entry(frame, width=5)
                entry.pack(side=tk.LEFT)
                entry.insert(0, "1")
                row_entries.append(entry)
                if j < self.num_variables - 1:
                    tk.Label(frame, text=f"x{i+1} +").pack(side=tk.LEFT)
                else:
                    tk.Label(frame, text=f"x{i+1}").pack(side=tk.LEFT)
            ineq = ttk.Combobox(frame, values=["<=", ">=", "="], width=3)
            ineq.pack(side=tk.LEFT, padx=5)
            ineq.set("<=")
            self.constraint_ineq_vars.append(ineq)
            rhs = tk.Entry(frame, width=5)
            rhs.pack(side=tk.LEFT)
            rhs.insert(0, "0")
            row_entries.append(rhs)
            self.constraint_entries.append(row_entries)

    def modify_variables(self):
        try:
            self.num_variables = int(self.var_entry.get())
            self.create_objective_entries()
            self.create_constraint_entries()
        except ValueError:
            messagebox.showerror("Error", "Enter a valid number of variables.")

    def modify_constraints(self):
        try:
            self.num_constraints = int(self.con_entry.get())
            self.create_constraint_entries()
        except ValueError:
            messagebox.showerror("Error", "Enter a valid number of constraints.")

    def run_simplex(self):
        try:
            m = self.num_variables
            n = self.num_constraints

            objective = []
            for entry in self.objective_entries:
                objective.append(float(entry.get()))

            constraints_coeffs = []
            constraints_rhs = [] 

            constraints_type = []
            for i, row in enumerate(self.constraint_entries):
                coeffs = []
                for j in range(m):
                    coeffs.append(float(row[j].get()))  

                rhs_val = float(row[-1].get())   
                coeffs.append(rhs_val) 

                constraints_coeffs.append(coeffs) 
                sign = self.constraint_ineq_vars[i].get()

                if sign == "<=":
                    constraints_type.append(1)
                elif sign == "=":
                    constraints_type.append(2)
                elif sign == ">=":
                    constraints_type.append(3)
                else:
                    constraints_type.append(1)  

            is_maximization = self.problem_type.get() == "Maximize"
            st = SimplexSolver(objective, constraints_coeffs, constraints_type, is_maximization, tolerance=1e-6)

            solution = st.solve()

            if isinstance(solution, str) and "Unbounded solution" in solution:
                solution_text = "Solution is unbounded"
                iteration_text = st.log
            elif isinstance(solution, str) and "Infeasible solution" in solution:
                solution_text = "Solution is infeasible"
                iteration_text = st.log
            else:
                solution_text = "Solution:\n"
                if solution and "solution" in solution:
                   for i, val in enumerate(solution["solution"]):
                        solution_text += f"x{i+1} = {val:.2f}\n"
                else:
                    solution_text = "No solution found.\n"

                if solution and "objective_value" in solution:
                    solution_text += f"Objective value: {solution['objective_value']:.2f}\n"
                else:
                    solution_text += "Objective value: N/A\n"

                iteration_text = "Iterations:\n"
                iteration_text += st.log
            self.result_text.delete("1.0", tk.END)
            self.result_text.insert(tk.END, solution_text + "\n" + iteration_text)

        except Exception as e:
            messagebox.showerror("Error", str(e))

root = tk.Tk()
gui = SimplexSolverGUI(root)
root.mainloop()
