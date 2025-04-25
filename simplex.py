import tkinter as tk
from tkinter import ttk, scrolledtext

class SimplexSolverGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Simplex Solver with Big M Method")
        self.root.geometry("1920x1080")

        # Variables for user input
        self.num_vars = tk.StringVar(value="2")
        self.num_constraints = tk.StringVar(value="2")
        self.objective_entries = []
        self.constraint_entries = []
        self.constraint_signs = []
        self.constraint_rhs = []
        self.maximize_var = tk.BooleanVar(value=True)

        # Main layout
        main_frame = ttk.Frame(root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Input for number of variables and constraints
        var_frame = ttk.Frame(main_frame)
        var_frame.pack(fill=tk.X, pady=5)
        ttk.Label(var_frame, text="Number of variables:").pack(side=tk.LEFT, padx=5)
        ttk.Entry(var_frame, textvariable=self.num_vars, width=5).pack(side=tk.LEFT)
        ttk.Label(var_frame, text="Number of constraints:").pack(side=tk.LEFT, padx=10)
        ttk.Entry(var_frame, textvariable=self.num_constraints, width=5).pack(side=tk.LEFT)
        ttk.Button(var_frame, text="Set", command=self.setup_inputs).pack(side=tk.LEFT, padx=10)

        # Maximize/Minimize option
        obj_type_frame = ttk.Frame(main_frame)
        obj_type_frame.pack(fill=tk.X, pady=5)
        ttk.Radiobutton(obj_type_frame, text="Maximize", variable=self.maximize_var, value=True).pack(side=tk.LEFT, padx=5)
        ttk.Radiobutton(obj_type_frame, text="Minimize", variable=self.maximize_var, value=False).pack(side=tk.LEFT)

        # Objective function input
        self.objective_frame = ttk.LabelFrame(main_frame, text="Objective Function")
        self.objective_frame.pack(fill=tk.X, pady=5)

        # Constraints input
        self.constraints_frame = ttk.LabelFrame(main_frame, text="Constraints")
        self.constraints_frame.pack(fill=tk.X, pady=5)

        ttk.Button(main_frame, text="Solve", command=self.solve_simplex).pack(pady=10)

        # Output text box for displaying results
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

            obj_label = "Max Z = " if self.maximize_var.get() else "Min Z = "
            ttk.Label(self.objective_frame, text=obj_label).pack(side=tk.LEFT)
            
            for i in range(num_vars):
                entry = ttk.Entry(self.objective_frame, width=5)
                entry.insert(0, "1")
                entry.pack(side=tk.LEFT)
                
                if i < num_vars - 1:
                    ttk.Label(self.objective_frame, text=f"x{i+1} +").pack(side=tk.LEFT)
                else:
                    ttk.Label(self.objective_frame, text=f"x{i+1}").pack(side=tk.LEFT)
                    
                self.objective_entries.append(entry)

            # Setup constraint inputs
            for i in range(num_constraints):
                row_frame = ttk.Frame(self.constraints_frame)
                row_frame.pack(fill=tk.X, pady=2)
                row_entries = []
                
                for j in range(num_vars):
                    entry = ttk.Entry(row_frame, width=5)
                    entry.insert(0, "1")
                    entry.pack(side=tk.LEFT)
                    
                    if j < num_vars - 1:
                        ttk.Label(row_frame, text=f"x{j+1} +").pack(side=tk.LEFT)
                    else:
                        ttk.Label(row_frame, text=f"x{j+1}").pack(side=tk.LEFT)
                        
                    row_entries.append(entry)
                    
                sign_var = tk.StringVar(value="<=")
                sign_menu = ttk.Combobox(row_frame, textvariable=sign_var, values=["<=", ">=", "="], width=5, state="readonly")
                sign_menu.pack(side=tk.LEFT, padx=5)
                
                rhs = ttk.Entry(row_frame, width=5)
                rhs.insert(0, "10")
                rhs.pack(side=tk.LEFT, padx=5)

                self.constraint_entries.append(row_entries)
                self.constraint_signs.append(sign_var)
                self.constraint_rhs.append(rhs)

        except ValueError:
            self.output_text.delete(1.0, tk.END)
            self.output_text.insert(tk.END, "Please enter valid numbers for variables and constraints.")

    def solve_simplex(self):
        """Prepares and solves the LP problem using the Simplex method."""
        try:
            self.output_text.delete(1.0, tk.END)
            num_vars = int(self.num_vars.get())
            num_constraints = int(self.num_constraints.get())
            
            # Get objective coefficients
            objective = [float(entry.get()) for entry in self.objective_entries]
            maximize = self.maximize_var.get()
            
            # Store original objective coefficients for later use
            original_objective = objective.copy()
            
            # If minimizing, negate coefficients to convert to maximization problem
            if not maximize:
                objective = [-c for c in objective]
            
            # Convert constraints to standard form
            constraints_matrix = []
            constraint_types = []
            rhs_values = []
            
            for i in range(num_constraints):
                row = [float(entry.get()) for entry in self.constraint_entries[i]]
                sign = self.constraint_signs[i].get()
                rhs = float(self.constraint_rhs[i].get())
                
                # Handle negative RHS by multiplying constraint by -1
                if rhs < 0:
                    row = [-coef for coef in row]
                    if sign == "<=":
                        sign = ">="
                    elif sign == ">=":
                        sign = "<="
                    rhs = -rhs
                
                constraints_matrix.append(row)
                constraint_types.append(sign)
                rhs_values.append(rhs)
            
            # Build the full tableau and solve
            result = self.build_and_solve_tableau(
                num_vars, num_constraints, objective, original_objective, constraints_matrix, 
                constraint_types, rhs_values, maximize
            )
            
            self.output_text.insert(tk.END, result)
            
        except Exception as e:
            self.output_text.insert(tk.END, f"Error: {str(e)}\n")
            import traceback
            self.output_text.insert(tk.END, traceback.format_exc())

    def build_and_solve_tableau(self, num_vars, num_constraints, objective, original_objective, 
                              constraints_matrix, constraint_types, rhs_values, maximize):
        """Builds the initial tableau and solves using simplex method."""
        output = "Building Initial Tableau:\n"
        
        # Count artificial and slack variables needed
        num_slack = sum(1 for ctype in constraint_types if ctype == "<=")
        num_surplus = sum(1 for ctype in constraint_types if ctype == ">=")
        num_artificial = sum(1 for ctype in constraint_types if ctype in [">=", "="])
        
        # Total number of variables in the tableau
        total_vars = num_vars + num_slack + num_surplus + num_artificial
        
        # Prepare the tableau
        tableau = []
        basic_var_indices = []  # Indices of basic variables for each constraint
        artificial_indices = []  # Indices of artificial variables
        
        # First row is for the objective function (will be filled later)
        tableau.append([0] * (total_vars + 1))  # +1 for RHS
        
        # Big M value
        M = 1000  # A large value for the Big M method
        
        # Initialize variable positions
        slack_start = num_vars
        surplus_start = slack_start + num_slack
        artificial_start = surplus_start + num_surplus
        
        slack_count = 0
        surplus_count = 0
        artificial_count = 0
        
        # Build constraint rows in the tableau
        for i in range(num_constraints):
            row = [0] * (total_vars + 1)
            
            # Original variables
            for j in range(num_vars):
                row[j] = constraints_matrix[i][j]
                
            # RHS
            row[-1] = rhs_values[i]
            
            # Add slack, surplus, and artificial variables as needed
            if constraint_types[i] == "<=":
                row[slack_start + slack_count] = 1
                basic_var_indices.append(slack_start + slack_count)
                slack_count += 1
            elif constraint_types[i] == ">=":
                row[surplus_start + surplus_count] = -1
                row[artificial_start + artificial_count] = 1
                basic_var_indices.append(artificial_start + artificial_count)
                artificial_indices.append(artificial_start + artificial_count)
                surplus_count += 1
                artificial_count += 1
            elif constraint_types[i] == "=":
                row[artificial_start + artificial_count] = 1
                basic_var_indices.append(artificial_start + artificial_count)
                artificial_indices.append(artificial_start + artificial_count)
                artificial_count += 1
            
            tableau.append(row)
        
        # Fill objective row (Z row) - FIX: Correct negative signs for maximization problem
        for j in range(num_vars):
            tableau[0][j] = -objective[j]  # Always negative for maximize format
        
        # Apply Big M method for artificial variables - FIX: correct sign for Big M penalty
        for j in artificial_indices:
            tableau[0][j] = M  # Large positive penalty for artificial variables
        
        # FIX: Apply Big M method properly by subtracting artificial variable rows from objective
        for i in range(1, len(tableau)):
            for j in artificial_indices:
                if tableau[i][j] != 0:  # If this constraint has an artificial variable
                    # Subtract M times this row from objective row to ensure artificial vars = 0 at optimality
                    for k in range(total_vars + 1):
                        tableau[0][k] -= M * tableau[i][k]
                    break  # Only subtract each row once
        
        # Print initial tableau
        output += "\n--- INITIAL TABLEAU ---\n"
        output += self.tableau_to_string(tableau, basic_var_indices, total_vars, num_vars)
        output += "---------------------\n"
        
        # Solve using simplex method
        iterations, tableaux_history, final_basic_vars, final_tableau = self.simplex_algorithm(
            tableau, basic_var_indices, total_vars, num_vars
        )
        
        # Add tableaux history to output
        for i, (iter_tableau, iter_basic_vars, pivot_info) in enumerate(tableaux_history):
            output += f"\n--- ITERATION {i+1}: {pivot_info} ---\n"
            output += self.tableau_to_string(iter_tableau, iter_basic_vars, total_vars, num_vars)
            output += "---------------------\n"
        
        # Get the solution - FIX: Correctly extract values from the tableau
        solution = [0] * num_vars
        for i, basic_idx in enumerate(final_basic_vars):
            if basic_idx < num_vars:  # Only for original variables
                solution[basic_idx] = final_tableau[i+1][-1]  # Extract from RHS of constraint row
        
        # Compute optimal objective value - FIX: Use the final tableau objective value
        obj_value = final_tableau[0][-1]  # Z value is already correctly calculated in the tableau
        
        if not maximize:
            obj_value = -obj_value  # Convert back for minimization problems
        
        # Check for artificial variables in final solution
        has_artificial_in_basis = False
        for i, idx in enumerate(final_basic_vars):
            if idx in artificial_indices and final_tableau[i+1][-1] > 0.00001:
                has_artificial_in_basis = True
                break
        
        # Output the solution
        output += f"\n--- FINAL SOLUTION ---\n"
        if has_artificial_in_basis:
            output += "WARNING: Artificial variables remain in the basis with positive values.\n"
            output += "This indicates that the problem has no feasible solution.\n"
        else:
            for i in range(num_vars):
                output += f"x{i+1} = {solution[i]:.4f}\n"
            
            output += f"\nOptimal {'Maximum' if maximize else 'Minimum'} Value: {obj_value:.4f}\n"
            
            # Also show how the objective was calculated using original coefficients
            output += "\nObjective Calculation:\n"
            terms = []
            for i in range(num_vars):
                if abs(solution[i]) > 0.00001:  # Only show non-zero terms
                    coef = original_objective[i]
                    terms.append(f"{coef:.4f} × {solution[i]:.4f}")
            
            if terms:
                if maximize:
                    output += " + ".join(terms) + f" = {obj_value:.4f}\n"
                else:
                    output += " + ".join(terms) + f" = {obj_value:.4f}\n"
            else:
                output += "0 = 0\n"
        
        return output

    def simplex_algorithm(self, tableau, basic_var_indices, total_vars, num_vars):
        """Implements the simplex algorithm with iteration history."""
        iterations = 0
        max_iterations = 100  # Prevent infinite loops
        tableaux_history = []  # Store tableau and basic variables at each iteration
        
        while iterations < max_iterations:
            # Find pivot column (most negative entry in objective row)
            pivot_col = -1
            min_value = -0.00001  # Small negative threshold to account for floating point errors
            
            for j in range(len(tableau[0]) - 1):  # Skip RHS
                if tableau[0][j] < min_value:  # For maximization problems
                    min_value = tableau[0][j]
                    pivot_col = j
            
            # If no negative values in objective row, we have reached optimality
            if pivot_col == -1:
                return iterations, tableaux_history, basic_var_indices, tableau
            
            # Find pivot row using minimum ratio test
            pivot_row = -1
            min_ratio = float('inf')
            
            for i in range(1, len(tableau)):  # Skip objective row
                if tableau[i][pivot_col] > 0.00001:  # Small positive threshold
                    ratio = tableau[i][-1] / tableau[i][pivot_col]
                    if ratio < min_ratio:
                        min_ratio = ratio
                        pivot_row = i
            
            # If no positive coefficient found, problem is unbounded
            if pivot_row == -1:
                pivot_info = "Unbounded Problem"
                tableaux_history.append((tableau.copy(), basic_var_indices.copy(), pivot_info))
                return iterations, tableaux_history, basic_var_indices, tableau
            
            # Record pivot information
            entering_var = self.get_variable_name(pivot_col, num_vars)
            leaving_var = self.get_variable_name(basic_var_indices[pivot_row - 1], num_vars)
            pivot_info = f"Pivot at ({pivot_row}, {pivot_col}) - {entering_var} enters, {leaving_var} leaves"
            
            # Store current tableau before pivoting
            current_tableau = [row[:] for row in tableau]  # Deep copy of tableau
            current_basic_vars = basic_var_indices.copy()
            
            # Update basic variable
            basic_var_indices[pivot_row - 1] = pivot_col  # -1 because row 0 is objective
            
            # Perform pivot operation
            self.pivot_operation(tableau, pivot_row, pivot_col)
            
            # Store the tableau after this iteration
            tableaux_history.append((current_tableau, current_basic_vars, pivot_info))
            
            iterations += 1
        
        pivot_info = "Max iterations reached"
        tableaux_history.append((tableau.copy(), basic_var_indices.copy(), pivot_info))
        return iterations, tableaux_history, basic_var_indices, tableau

    def get_variable_name(self, idx, num_vars):
        """Returns a meaningful name for variables based on their index."""
        if idx < num_vars:
            return f"x{idx+1}"
        else:
            # For slack, surplus, and artificial variables
            var_type = "s"  # Default to slack
            
            # Count how many slack variables there would be
            num_slack = 0
            num_surplus = 0
            for sign_var in self.constraint_signs:
                sign = sign_var.get()
                if sign == "<=":
                    num_slack += 1
                elif sign == ">=":
                    num_surplus += 1
                    
            # Determine variable type based on index
            if idx < num_vars + num_slack:
                var_type = "s"  # Slack
            elif idx < num_vars + num_slack + num_surplus:
                var_type = "e"  # Surplus (excess)
            else:
                var_type = "a"  # Artificial
            
            # Calculate the proper variable number
            if var_type == "s":
                var_num = idx - num_vars + 1
            elif var_type == "e":
                var_num = idx - (num_vars + num_slack) + 1
            else:
                var_num = idx - (num_vars + num_slack + num_surplus) + 1
                
            return f"{var_type}{var_num}"

    def pivot_operation(self, tableau, pivot_row, pivot_col):
        """Performs the pivot operation on the tableau."""
        pivot_value = tableau[pivot_row][pivot_col]
        
        # Normalize pivot row
        for j in range(len(tableau[pivot_row])):
            tableau[pivot_row][j] /= pivot_value
        
        # Update other rows
        for i in range(len(tableau)):
            if i != pivot_row:
                factor = tableau[i][pivot_col]
                for j in range(len(tableau[i])):
                    tableau[i][j] -= factor * tableau[pivot_row][j]

    def tableau_to_string(self, tableau, basic_vars, total_vars, num_vars):
        """Converts the tableau to a formatted string for display with implicit alignment."""
        output = ""

        # Header
        header = f"{'Basis':<10}"  # Adjust initial spacing as needed
        for j in range(total_vars):
            header += f"{self.get_variable_name(j, num_vars):<10}" # Adjust initial spacing as needed
        header += f"{'RHS':<10}\n" # Adjust initial spacing as needed
        output += header

        # Row 0 (objective function)
        output += f"{'Z':<10}"
        for j in range(len(tableau[0])):
            output += f"{tableau[0][j]:<10.3f}"
        output += "\n"

        # Constraint rows
        for i in range(1, len(tableau)):
            basic_idx = basic_vars[i-1]
            var_name = self.get_variable_name(basic_idx, num_vars)

            output += f"{var_name:<10}"
            for j in range(len(tableau[i])):
                output += f"{tableau[i][j]:<10.3f}"
            output += "\n"

        return output

root = tk.Tk()
app = SimplexSolverGUI(root)
root.mainloop()