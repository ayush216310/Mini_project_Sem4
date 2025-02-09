import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import pandas as pd


class LinearOptimizerGUI:
    def __init__(self, master):
        self.master = master
        master.title("Linear Optimization - Simplex Method")

        # default values
        self.num_variables = 2  # Number of decision variables.
        self.num_constraints = 2  # Number of constraints.

        self.objective_entries = []
        self.constraint_entries = []
        self.constraint_ineq_vars = []

        tk.Label(master, text="Number of Variables:").pack()
        self.var_entry = tk.Entry(master)
        self.var_entry.pack()
        self.var_entry.insert(0, str(self.num_variables))
        tk.Button(master, text="Set Variables", command=self.modify_variables).pack()

        tk.Label(master, text="Number of Constraints:").pack()
        self.con_entry = tk.Entry(master)
        self.con_entry.pack()
        self.con_entry.insert(0, str(self.num_constraints))
        tk.Button(
            master, text="Set Constraints", command=self.modify_constraints
        ).pack()

        tk.Label(master, text="Objective Function (Maximize):").pack()
        self.objective_frame = tk.Frame(master)
        self.objective_frame.pack()
        self.create_objective_entries()

        tk.Label(master, text="Constraints:").pack()
        self.constraints_frame = tk.Frame(master)
        self.constraints_frame.pack()
        self.create_constraint_entries()

        self.solve_button = tk.Button(master, text="Solve", command=self.run_simplex)
        self.solve_button.pack()

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
            m = self.num_variables  # Number of decision variables.
            n = self.num_constraints  # Number of constraints.

            # --- Gather Objective Function Coefficients ---
            objective = []
            for entry in self.objective_entries:
                objective.append(float(entry.get()))

            # --- Gather Constraint Data ---
            constraints_coeffs = []
            constraints_rhs = []
            constraints_type = []  # 1 for ">=", 2 for "<=", 3 for "="
            for i, row in enumerate(self.constraint_entries):
                coeffs = []
                for j in range(m):
                    coeffs.append(float(row[j].get()))
                rhs_val = float(row[-1].get())
                constraints_coeffs.append(coeffs)
                constraints_rhs.append(rhs_val)
                sign = self.constraint_ineq_vars[i].get()
                if sign == "<=":
                    constraints_type.append(2)
                elif sign == ">=":
                    constraints_type.append(1)
                elif sign == "=":
                    constraints_type.append(3)
                else:
                    constraints_type.append(2)

            # --- Build the Initial Simplex Tableau ---
            row_lim = n + 1
            col_lim = m + n + 2
            table = np.zeros((row_lim, col_lim))
            for i in range(n):
                for j in range(m):
                    table[i][j] = constraints_coeffs[i][j]
                table[i][col_lim - 1] = constraints_rhs[i]
            for i in range(n):
                if constraints_type[i] != 3:
                    if constraints_type[i] == 1:  # for ">=" add surplus (–1)
                        table[i][m + i] = -1
                    else:  # for "<=" add slack (+1)
                        table[i][m + i] = 1

            for j in range(m):
                table[n][j] = -1 * objective[j]
            table[n][col_lim - 2] = 1
            table[n][col_lim - 1] = 0

            row_labels = [f"Constraint-{i+1}" for i in range(n)] + ["Function"]
            col_labels = (
                [f"x{j+1}" for j in range(m)]
                + [f"s{i+1}" for i in range(n)]
                + ["z", "c"]
            )

            # --- Helper Functions (with tolerance) ---
            tol = 1e-8

            def check(last_row):
                min_val = 0
                pos = -1
                for i in range(col_lim - 1):
                    if last_row[i] < -tol and last_row[i] < min_val:
                        min_val = last_row[i]
                        pos = i
                return pos

            def check_ratio(table, pivot_col):
                min_ratio = float("inf")
                pivot_row = -1
                for i in range(n):
                    if table[i][pivot_col] > tol:
                        ratio = table[i][col_lim - 1] / table[i][pivot_col]
                        # Here we add the extra condition that ratio must be > 0.
                        if ratio < min_ratio and ratio > 0:
                            min_ratio = ratio
                            pivot_row = i
                return pivot_row

            iteration_details = ""
            iteration = 0

            # Display initial tableau as ITERATION 0.
            df = pd.DataFrame(table, index=row_labels, columns=col_labels)
            iteration_details += f"ITERATION {iteration}:\n" + df.to_string() + "\n\n"

            while True:
                pivot_col = check(table[n])
                if pivot_col == -1:
                    # No negative coefficient found → optimality reached.
                    break
                pivot_row = check_ratio(table, pivot_col)
                if pivot_row == -1:
                    iteration_details += "The solution is unbounded.\n"
                    break
                pivot = table[pivot_row][pivot_col]
                table[pivot_row, :] = table[pivot_row, :] / pivot
                for i in range(row_lim):
                    if i != pivot_row:
                        table[i, :] = (
                            table[i, :] - table[i][pivot_col] * table[pivot_row, :]
                        )
                iteration += 1
                iteration_details += f"ITERATION {iteration}:\n"
                df = pd.DataFrame(table, index=row_labels, columns=col_labels)
                iteration_details += df.to_string() + "\n\n"

            solution_text = "Solution:\n"
            for j in range(m):
                column = table[:n, j]
                ones = np.isclose(column, 1, atol=tol)
                zeros = np.isclose(column, 0, atol=tol)
                if np.sum(ones) == 1 and np.sum(zeros) == n - 1:
                    row_index = np.where(ones)[0][0]
                    value = table[row_index, col_lim - 1]
                else:
                    value = 0
                solution_text += f"x{j+1} = {value:.2f}\n"
            max_z = table[n][col_lim - 1]
            solution_text += f"\nMaximized value of z: {max_z:.2f}\n"

            self.result_text.delete("1.0", tk.END)
            self.result_text.insert(tk.END, iteration_details + "\n" + solution_text)

        except Exception as e:
            messagebox.showerror("Error", str(e))


root = tk.Tk()
gui = LinearOptimizerGUI(root)
root.mainloop()
