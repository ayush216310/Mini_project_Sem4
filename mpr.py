import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np

class LinearOptimizerGUI:
    def __init__(self, master):
        self.master = master
        master.title("Linear Optimization - Simplex Method")

        self.num_variables = 2
        self.objective_entries = []
        self.constraint_entries = []
        self.constraint_ineq_vars = []
        self.constraint_frames = []
        self.objective_function_type = "Maximum"

        tk.Label(master, text="Number of Variables:").pack()
        self.var_entry = tk.Entry(master)
        self.var_entry.pack()
        self.var_entry.insert(0, "2")
        tk.Button(master, text="Set Variables", command=self.modify_variables).pack()

        tk.Label(master, text="Objective Function (Maximize):").pack()
        self.objective_frame = tk.Frame(master)
        self.objective_frame.pack()
        self.create_objective_entries()

        tk.Label(master, text="Constraints:").pack()
        self.constraints_frame = tk.Frame(master)
        self.constraints_frame.pack()
        self.create_constraint_entries()

        self.solve_button = tk.Button(master, text="Solve", command=self.solve)
        self.solve_button.pack()
        self.steps_button = tk.Button(master, text="Show Steps", command=self.show_steps)
        self.steps_button.pack()

        self.result_text = tk.Text(master, height=15, width=60)
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
            tk.Label(self.objective_frame, text=f"x{i+1} +").pack(side=tk.LEFT)

    def create_constraint_entries(self):
        for widget in self.constraints_frame.winfo_children():
            widget.destroy()
        self.constraint_entries = []
        self.constraint_ineq_vars = []
        for i in range(self.num_variables):
            frame = tk.Frame(self.constraints_frame)
            frame.pack()
            row = []
            for j in range(self.num_variables):
                entry = tk.Entry(frame, width=5)
                entry.pack(side=tk.LEFT)
                entry.insert(0, "1")
                row.append(entry)
                tk.Label(frame, text=f"x{j+1} +").pack(side=tk.LEFT)
            ineq = ttk.Combobox(frame, values=["<=", ">=", "="], width=3)
            ineq.pack(side=tk.LEFT)
            ineq.set("<=")
            rhs = tk.Entry(frame, width=5)
            rhs.pack(side=tk.LEFT)
            rhs.insert(0, "0")
            row.append(rhs)
            self.constraint_entries.append(row)
            self.constraint_ineq_vars.append(ineq)

    def modify_variables(self):
        try:
            self.num_variables = int(self.var_entry.get())
            self.create_objective_entries()
            self.create_constraint_entries()
        except ValueError:
            messagebox.showerror("Error", "Enter a valid number of variables.")

    def solve(self):
        self.result_text.delete("1.0", tk.END)
        self.result_text.insert(tk.END, "Solution process will be shown in steps. Click 'Show Steps'.")

    def show_steps(self):
        self.result_text.delete("1.0", tk.END)
        try:
            obj = [-float(e.get()) for e in self.objective_entries]
            constraints = []
            rhs = []
            for i, row in enumerate(self.constraint_entries):
                constraints.append([float(e.get()) for e in row[:-1]])
                rhs.append(float(row[-1].get()))

            A = np.array(constraints, dtype=float)
            b = np.array(rhs, dtype=float)
            c = np.array(obj, dtype=float)
            
            num_slack = len(b)
            tableau = np.hstack([A, np.eye(num_slack), b.reshape(-1, 1)])
            cost_row = np.hstack([c, np.zeros(num_slack + 1)])
            tableau = np.vstack([tableau, cost_row])
            
            self.result_text.insert(tk.END, "Initial Simplex Tableau:\n")
            self.display_tableau(tableau)
        except ValueError:
            messagebox.showerror("Error", "Invalid input values.")

    def display_tableau(self, tableau):
        for row in tableau:
            self.result_text.insert(tk.END, "\t".join(map(lambda x: f"{x:.2f}", row)) + "\n")

root = tk.Tk()
gui = LinearOptimizerGUI(root)
root.mainloop()
