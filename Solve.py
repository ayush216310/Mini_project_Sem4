import tkinter as tk  
from tkinter import ttk, messagebox

class SimplexTable:
    def __init__(self, num_constraints, num_variables, objective_coefficients, constraint_coefficients, rhs_values, constraint_types):
        self.num_constraints = num_constraints
        self.num_variables = num_variables
        self.objective_coefficients = objective_coefficients
        self.constraint_coefficients = constraint_coefficients
        self.rhs_values = rhs_values
        self.constraint_types = constraint_types

        self.num_slack_surplus = self.num_constraints
        self.num_artificial = 0

        self.num_total_variables = self.num_variables + self.num_slack_surplus + self.num_artificial
        self.num_rows = self.num_constraints + 1

        self.tableau = [[0.0 for _ in range(self.num_total_variables + 1)]
                        for _ in range(self.num_rows)]

        self.variable_names = []
        for i in range(self.num_variables):
            self.variable_names.append(f"x{i+1}")

        for i in range(self.num_slack_surplus):
            self.variable_names.append(f"s{i+1}")

        for i in range(self.num_artificial):
            self.variable_names.append(f"a{i+1}")

        for i in range(self.num_variables):
            self.tableau[self.num_rows - 1][i] = - self.objective_coefficients[i]

        for i in range(self.num_constraints):
            for j in range(self.num_variables):
                self.tableau[i][j] = self.constraint_coefficients[i][j]
            self.tableau[i][self.num_total_variables] = self.rhs_values[i]

            if constraint_types[i] == 3: 
                self.tableau[i][self.num_variables + i] = 1
            elif constraint_types[i] == 1: 
                self.tableau[i][self.num_variables + i] = -1

        print("Tableau after data entry:", self.tableau)

    def find_pivot_column(self):
        objective_row = self.tableau[-1]
        most_negative = 0
        pivot_column = -1

        for i in range(len(objective_row) - 1):
            if objective_row[i] < most_negative:
                most_negative = objective_row[i]
                pivot_column = i

        return pivot_column

    def find_pivot_row(self, pivot_column):
        min_ratio = 1e9
        pivot_row = -1
        num_rows = len(self.tableau) - 1

        for i in range(num_rows):
            if self.tableau[i][pivot_column] > 0:
                ratio = self.tableau[i][-1] / self.tableau[i][pivot_column]
                if ratio < min_ratio:
                    min_ratio = ratio
                    pivot_row = i

        return pivot_row

    def perform_row_operations(self, pivot_row, pivot_column):
        pivot_value = self.tableau[pivot_row][pivot_column]

        for j in range(len(self.tableau[pivot_row])):
            self.tableau[pivot_row][j] /= pivot_value

        for i in range(len(self.tableau)):
            if i != pivot_row:
                factor = self.tableau[i][pivot_column]
                for j in range(len(self.tableau[i])):
                    self.tableau[i][j] -= factor * self.tableau[pivot_row][j]

    def solve(self):
        iterations = []
        while True:
            iterations.append([row[:] for row in self.tableau])  
            pivot_col = self.find_pivot_column()
            if pivot_col == -1:
                break  

            pivot_row = self.find_pivot_row(pivot_col)

            if pivot_row == -1:
                print("Unbounded")
                return {"status": "unbounded", "iterations": iterations}  

            self.perform_row_operations(pivot_row, pivot_col)
            print("Tableau after step:", self.tableau)

        
        solution = {}
        
        tol = 1e-8
        for j in range(self.num_variables):
            column = [self.tableau[i][j] for i in range(self.num_constraints)]
            ones = [abs(x - 1) < tol for x in column] 
            zeros = [abs(x) < tol for x in column] 
            if sum(ones) == 1 and sum(zeros) == self.num_constraints - 1:
                row_index = ones.index(True)  
                value = self.tableau[row_index][self.num_total_variables]
            else:
                value = 0
            solution[f"x{j+1}"] = value

        objective_value = self.tableau[-1][-1]  
        solution["objective"] = objective_value
        solution["iterations"] = iterations

        return solution    
    
class SimplexSolverGUI:
    def __init__(self, master):
        self.master = master
        master.title("Simplex Solver - Simplex Method")
        
        self.num_variables = 3
        self.num_constraints = 3
        
        self.objective_entries = []
        self.costraint_entries = []
        self.constraint_ineq_entries = []
        
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
                    tk.Label(frame, text=f"x{j+1} +").pack(side=tk.LEFT)
                else:
                    tk.Label(frame, text=f"x{j+1}").pack(side=tk.LEFT)
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
                constraints_coeffs.append(coeffs)
                constraints_rhs.append(rhs_val)
                sign =self.constraint_ineq_vars[i].get()
                if sign == "<=":
                    constraints_type.append(2)
                elif sign == ">=":
                    constraints_type.append(1)
                elif sign == "=":
                    constraints_type.append(3)
                else:
                    constraints_type.append(2)
                    
            st = SimplexTable(n, m, objective, constraints_coeffs, constraints_rhs, constraints_type)
            
            solution = st.solve()
            
            if "status" in solution and solution["status"] == "unbounded":
                solution_text = "solution is unbounded"
                iteration_text = ""
            else:
                solution_text = "solution:\n"
                for var, val in solution.items():
                    if var != "objective" and var != "iterations":
                        solution_text += f"{var} = {val:.2f}\n"
                solution_text += f"objective value: {solution['objective']:.2f}\n"
                
                iteration_text = "Iterations:\n"
                for idx, tableau in enumerate(solution["iterations"]):
                    iteration_text += f"Iteration {idx + 1}:\n"
                    for row in tableau:
                        iteration_text += "\t".join(f"{val:.2f}" for val in row) + "\n"
                    iteration_text += "\n"
            
            self.result_text.delete("1.0", tk.END)
            self.result_text.insert(tk.END, solution_text + "\n" + iteration_text)
        
        except Exception as e:
            messagebox.showerror("Error", str(e))                    
        
root = tk.Tk()
gui = SimplexSolverGUI(root)
root.mainloop()
