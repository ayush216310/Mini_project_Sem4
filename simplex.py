import tkinter as tk
from tkinter import ttk, scrolledtext

class SimplexSolverGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Simplex Solver with Big M Method")
        self.root.geometry("800x600")

        # Variables
        self.num_vars = tk.StringVar(value="2")
        self.num_constraints = tk.StringVar(value="2")
        self.objective_entries = []
        self.constraint_entries = []
        self.constraint_signs = []
        self.constraint_rhs = []

        # Main layout
        main_frame = ttk.Frame(root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        opt_frame = ttk.Frame(main_frame)
        opt_frame.pack(anchor=tk.W, pady=5)
        self.opt_type = tk.StringVar(value="maximize")
        ttk.Radiobutton(opt_frame, text="Maximize", variable=self.opt_type, value="maximize").pack(side=tk.LEFT, padx=10)
        ttk.Radiobutton(opt_frame, text="Minimize", variable=self.opt_type, value="minimize").pack(side=tk.LEFT, padx=10)

        var_frame = ttk.Frame(main_frame)
        var_frame.pack(fill=tk.X, pady=5)
        ttk.Label(var_frame, text="Number of variables:").pack(side=tk.LEFT, padx=5)
        ttk.Entry(var_frame, textvariable=self.num_vars, width=5).pack(side=tk.LEFT)
        ttk.Label(var_frame, text="Number of constraints:").pack(side=tk.LEFT, padx=10)
        ttk.Entry(var_frame, textvariable=self.num_constraints, width=5).pack(side=tk.LEFT)
        ttk.Button(var_frame, text="Set", command=self.setup_inputs).pack(side=tk.LEFT, padx=10)

        self.objective_frame = ttk.LabelFrame(main_frame, text="Objective Function")
        self.objective_frame.pack(fill=tk.X, pady=5)

        self.constraints_frame = ttk.LabelFrame(main_frame, text="Constraints")
        self.constraints_frame.pack(fill=tk.X, pady=5)

        ttk.Button(main_frame, text="Solve", command=self.solve_simplex).pack(pady=10)

        self.output_frame = ttk.LabelFrame(main_frame, text="Output")
        self.output_frame.pack(fill=tk.BOTH, expand=True, pady=10)
        self.output_text = scrolledtext.ScrolledText(self.output_frame, wrap=tk.WORD)
        self.output_text.pack(fill=tk.BOTH, expand=True)

        self.setup_inputs()

    def setup_inputs(self):
        for widget in self.objective_frame.winfo_children():
            widget.destroy()
        for widget in self.constraints_frame.winfo_children():
            widget.destroy()

        try:
            num_vars = int(self.num_vars.get())
            num_constraints = int(self.num_constraints.get())

            self.objective_entries = []
            self.constraint_entries = []
            self.constraint_signs = []
            self.constraint_rhs = []

            for i in range(num_vars):
                entry = ttk.Entry(self.objective_frame, width=5)
                entry.insert(0, "1")
                entry.pack(side=tk.LEFT)
                ttk.Label(self.objective_frame, text=f"x{i+1} +").pack(side=tk.LEFT)
                self.objective_entries.append(entry)
            ttk.Label(self.objective_frame, text="Z").pack(side=tk.LEFT)

            for i in range(num_constraints):
                row_frame = ttk.Frame(self.constraints_frame)
                row_frame.pack(fill=tk.X, pady=2)
                row_entries = []
                for j in range(num_vars):
                    entry = ttk.Entry(row_frame, width=5)
                    entry.insert(0, "1")
                    entry.pack(side=tk.LEFT)
                    ttk.Label(row_frame, text=f"x{j+1} +").pack(side=tk.LEFT)
                    row_entries.append(entry)
                sign_var = tk.StringVar(value="<=")
                sign_menu = ttk.Combobox(row_frame, textvariable=sign_var, values=["<=", ">=", "="], width=5, state="readonly")
                sign_menu.pack(side=tk.LEFT, padx=5)
                rhs = ttk.Entry(row_frame, width=5)
                rhs.insert(0, "0")
                rhs.pack(side=tk.LEFT, padx=5)

                self.constraint_entries.append(row_entries)
                self.constraint_signs.append(sign_var)
                self.constraint_rhs.append(rhs)

        except ValueError:
            pass

    def solve_simplex(self):
        try:
            self.output_text.delete(1.0, tk.END)

            num_vars = int(self.num_vars.get())
            objective = [float(entry.get()) for entry in self.objective_entries]
            maximize = self.opt_type.get() == "maximize"
            
            # Store original objective for final output calculation
            original_objective = objective[:]
            
            if not maximize:
                objective = [-c for c in objective]

            # Prepare the tableau structure
            num_constraints = len(self.constraint_entries)
            
            # Count slack and artificial variables
            num_slack = 0
            num_artificial = 0
            constraint_types = []  # Store constraint types for each row
            
            for i in range(num_constraints):
                sign = self.constraint_signs[i].get()
                constraint_types.append(sign)
                
                if sign == "<=":
                    num_slack += 1
                elif sign == ">=":
                    num_slack += 1
                    num_artificial += 1
                elif sign == "=":
                    num_artificial += 1
            
            # Initialize variable lists
            slack_vars = []  # Will hold the column indices of slack variables
            artificial_vars = []  # Will hold the column indices of artificial variables
            
            # Create the tableau
            tableau = []
            for i in range(num_constraints):
                coef = [float(entry.get()) for entry in self.constraint_entries[i]]
                rhs = float(self.constraint_rhs[i].get())
                sign = constraint_types[i]
                
                # Ensure RHS is non-negative
                if rhs < 0:
                    coef = [-c for c in coef]
                    rhs = -rhs
                    if sign == "<=":
                        sign = ">="
                        constraint_types[i] = sign
                    elif sign == ">=":
                        sign = "<="
                        constraint_types[i] = sign
                    
                # Create row with original variables
                row = coef[:]
                
                # Add slack variable columns
                for j in range(num_constraints):
                    if j == i:
                        if sign == "<=":
                            row.append(1.0)  # Slack for <=
                        elif sign == ">=":
                            row.append(-1.0)  # Surplus for >=
                        else:  # sign == "="
                            row.append(0.0)   # No slack for =
                    else:
                        row.append(0.0)
                
                # Add artificial variable columns
                for j in range(num_constraints):
                    if j == i and (sign == ">=" or sign == "="):
                        row.append(1.0)  # Artificial variable for >= or =
                    else:
                        row.append(0.0)
                
                # Add RHS
                row.append(rhs)
                tableau.append(row)
            
            # Create slack and artificial variable names
            slack_names = [f"s{i+1}" for i in range(num_slack)]
            artificial_names = [f"a{i+1}" for i in range(num_artificial)]
            
            # Track the column indices of slack and artificial variables
            slack_vars = list(range(num_vars, num_vars + num_slack))
            artificial_vars = list(range(num_vars + num_slack, num_vars + num_slack + num_artificial))
            
            # Create cost row (objective function)
            cost = objective + [0] * num_slack + [0] * num_artificial
            
            # Initialize basic variables and their costs
            basic_vars = []
            cb = []
            
            M = 1000000  # Big M value
            
            # Initialize tableau with basic variables
            art_index = 0
            slack_index = 0
            
            for i, sign in enumerate(constraint_types):
                if sign == "<=":
                    # Slack variable is basic
                    basic_vars.append(f"s{slack_index+1}")
                    cb.append(0)  # Cost for slack is 0
                    slack_index += 1
                elif sign == ">=":
                    # Artificial variable is basic
                    basic_vars.append(f"a{art_index+1}")
                    cb.append(M)  # Cost for artificial is M
                    # Make the corresponding artificial variable cost M in the objective
                    if art_index < len(artificial_vars):
                        cost[artificial_vars[art_index]] = M
                    art_index += 1
                    slack_index += 1  # Skip the surplus variable
                elif sign == "=":
                    # Artificial variable is basic
                    basic_vars.append(f"a{art_index+1}")
                    cb.append(M)  # Cost for artificial is M
                    # Make the corresponding artificial variable cost M in the objective
                    if art_index < len(artificial_vars):
                        cost[artificial_vars[art_index]] = M
                    art_index += 1
            
            # Create variable names for all columns
            var_names = [f"x{i+1}" for i in range(num_vars)] + slack_names + artificial_names
            
            output = self.simplex_solver(tableau, cost, basic_vars, cb, var_names, artificial_vars, original_objective, maximize)
            self.output_text.insert(tk.END, output)

        except Exception as e:
            self.output_text.insert(tk.END, f"Error: {str(e)}\n")
            import traceback
            self.output_text.insert(tk.END, traceback.format_exc())

    def simplex_solver(self, tableau, cost, basic_vars, cb, var_names, artificial_vars=None, original_objective=None, maximize=True, M=1000000):
        output = ""
        
        def print_table(step, zj_cj, ratios=None):
            nonlocal output
            output += f"\nSimplex Tableau - Step {step}\n"
            
            # Create header with proper variable names
            header = "Cb\tBasic\t" + "\t".join(var_names) + "\tRHS"
            output += header + "\n"
            
            for i, row in enumerate(tableau):
                line = f"{cb[i]:.2f}\t{basic_vars[i]}\t" + "\t".join([f"{x:.2f}" for x in row])
                if ratios and i < len(ratios):
                    line += f"\t{ratios[i]:.2f}" if ratios[i] != float('inf') else "\t-"
                output += line + "\n"
                
            if zj_cj:
                output += "Zj - Cj:\t\t" + "\t".join([f"{x:.2f}" for x in zj_cj]) + "\n"

        def calculate_zj_cj():
            zj = [0] * (len(tableau[0]) - 1) 
            for j in range(len(zj)):
                for i in range(len(tableau)):
                    zj[j] += cb[i] * tableau[i][j]
            
            zj_cj = [zj[j] - cost[j] for j in range(len(zj))]
            return zj_cj

        def calculate_ratios(pivot_col):
            ratios = []
            for i in range(len(tableau)):
                if tableau[i][pivot_col] <= 1e-10:  # Use epsilon to handle floating point errors
                    ratios.append(float('inf'))
                else:
                    ratios.append(tableau[i][-1] / tableau[i][pivot_col])
            return ratios

        def pivot(pivot_row, pivot_col):
            pivot_element = tableau[pivot_row][pivot_col]
            # Avoid division by zero
            if abs(pivot_element) < 1e-10:
                pivot_element = 1e-10
                
            # Scale the pivot row
            tableau[pivot_row] = [x / pivot_element for x in tableau[pivot_row]]
            
            # Update other rows
            for i in range(len(tableau)):
                if i != pivot_row:
                    factor = tableau[i][pivot_col]
                    tableau[i] = [a - factor * b for a, b in zip(tableau[i], tableau[pivot_row])]

        step = 1
        max_iterations = 100  # Safety limit to prevent infinite loops
        
        while step <= max_iterations:
            zj_cj = calculate_zj_cj()
            
            # Print current tableau
            print_table(step, zj_cj)
            
            # For maximization, we look for the most negative Zj-Cj value
            # (we're already handling minimization by inverting coefficients)
            if any(z < -1e-10 for z in zj_cj):  # Use a small epsilon to handle floating point errors
                pivot_col = min(range(len(zj_cj)), key=lambda i: zj_cj[i])
            else:
                # No negative values, optimal solution found
                break
                
            # Calculate ratios for minimum ratio test
            ratios = calculate_ratios(pivot_col)
            
            # Find minimum positive ratio - this is the leaving variable
            valid_ratios = [(i, r) for i, r in enumerate(ratios) if r >= 0 and r < float('inf')]
            if not valid_ratios:
                output += "\nUnbounded solution.\n"
                return output
                
            pivot_row = min(valid_ratios, key=lambda x: x[1])[0]
            
            # Update basic variable and its cost
            basic_vars[pivot_row] = var_names[pivot_col]
            cb[pivot_row] = cost[pivot_col]
            
            # Perform the pivot operation
            pivot(pivot_row, pivot_col)
            step += 1
            
        if step > max_iterations:
            output += "\nReached maximum iterations. Solution may not be optimal.\n"
            
        # Check if any artificial variables remain in the basis with non-zero values
        infeasible = False
        for i, var in enumerate(basic_vars):
            if var.startswith("a") and abs(tableau[i][-1]) > 1e-6:
                infeasible = True
                break
                
        if infeasible:
            output += "\nNo feasible solution found (artificial variables remain in the optimal solution).\n"
            return output
            
        # Extract solution values
        num_vars = len(original_objective)
        solution = [0] * num_vars
        
        for i, var in enumerate(basic_vars):
            if var.startswith("x") and var[1:].isdigit():
                index = int(var[1:]) - 1
                if 0 <= index < len(solution):
                    solution[index] = tableau[i][-1]
        
        # Calculate final objective value based on the original objective function
        z_value = sum(original_objective[i] * solution[i] for i in range(min(len(original_objective), len(solution))))
        
        # If minimizing, we need to negate the result back
        if not maximize:
            z_value = -z_value
        
        # Generate output
        output += "\nFinal Solution:\n"
        for i, val in enumerate(solution):
            output += f"x{i+1} = {val:.2f}\n"
        
        # Show values of slack and artificial variables
        slack_and_art_vars = {}
        for i, var in enumerate(basic_vars):
            if var.startswith("s") or var.startswith("a"):
                slack_and_art_vars[var] = tableau[i][-1]
        
        # Print slack variables first, then artificial
        for i, name in enumerate(var_names[num_vars:]):
            if name in slack_and_art_vars:
                output += f"{name} = {slack_and_art_vars[name]:.2f}\n"
            else:
                # Find the variable in the tableau if it's not basic
                found = False
                for row_idx, row in enumerate(tableau):
                    col_idx = num_vars + i
                    if col_idx < len(row) - 1 and abs(row[col_idx] - 1.0) < 1e-10:
                        if all(abs(tableau[r][col_idx]) < 1e-10 for r in range(len(tableau)) if r != row_idx):
                            output += f"{name} = 0.00\n"
                            found = True
                            break
                if not found:
                    output += f"{name} = 0.00\n"
        
        output += f"Z = {z_value:.2f} ({'maximized' if maximize else 'minimized'})\n"
        return output

root = tk.Tk()
app = SimplexSolverGUI(root)
root.mainloop()