import tkinter as tk
from tkinter import ttk, messagebox
import copy
import math

class SimplexSolver:
    def __init__(self, objective_func, constraint, ineq, is_maximization=True, tolerance=1e-6):
        self.is_min = not is_maximization
        self.sign = -1 if self.is_min else 1
        self.objective = [self.sign * c for c in objective_func]
        self.orig_objective = list(objective_func)
        self.constraints = [list(row) for row in constraint]
        self.ineq_codes = list(ineq)
        self.tol = tolerance
        self.log = ""
        self._build_tableau()

    def _fmt(self, x):
        return f"{x:9.3f}"

    def _build_tableau(self):
        m = len(self.constraints)
        n = len(self.objective)
        extra = []  
        for i, code in enumerate(self.ineq_codes):
            if code == 1:
                extra.append(('slack', i))
            elif code == 3:
                extra.append(('surplus', i))
                extra.append(('art', i))
            elif code == 2:
                extra.append(('art', i))
            else:
                raise ValueError(f"Unknown inequality code {code}")
        p = len(extra)
        # Build var names
        self.var_names = [f'x{j+1}' for j in range(n)]
        cnt = {'slack': 0, 'surplus': 0, 'art': 0}
        for typ, _ in extra:
            cnt[typ] += 1
            if typ == 'slack':     self.var_names.append(f's{cnt[typ]}')
            elif typ == 'surplus': self.var_names.append(f'sp{cnt[typ]}')
            elif typ == 'art':     self.var_names.append(f'a{cnt[typ]}')
        # Build tableau rows
        self.tableau = []
        for i, row in enumerate(self.constraints):
            coefs = row[:-1]
            b = row[-1]
            new = coefs + [0]*p + [b]
            for k, (typ, ci) in enumerate(extra):
                if ci == i:
                    col = n + k
                    if typ == 'slack':     new[col] = 1
                    elif typ == 'surplus': new[col] = -1
                    elif typ == 'art':     new[col] = 1
            self.tableau.append(new)
        # Objective row (last)
        cols = n + p + 1
        obj = [0]*cols
        for j, c in enumerate(self.objective):
            obj[j] = -c  
        M = 1e6
        for k, (typ, _) in enumerate(extra):
            if typ == 'art': obj[n + k] = -M
        obj[-1] = 0
        self.tableau.append(obj)
        # Basic vars
        self.basic_vars = []
        for i in range(m):
            for k, (typ, ci) in enumerate(extra):
                if ci == i and typ in ('slack', 'art'):
                    self.basic_vars.append(self.var_names[n + k])
                    break
        self._log_tableau()

    def _compute_zj_cj(self):
        cols = len(self.tableau[0])
        # build Cj array (internal maximize)
        cj = []
        # originals
        for c in self.objective:
            cj.append(c)
        # extras
        for name in self.var_names[len(self.objective):]:
            cj.append(-1e6 if name.startswith('a') else 0)
        cj.append(0)  # RHS
        zj = [0]*cols
        for i, bv in enumerate(self.basic_vars):
            bi = self.var_names.index(bv)
            for j in range(cols): zj[j] += cj[bi] * self.tableau[i][j]
        return [zj[j] - cj[j] for j in range(cols)]

    def _find_pivot(self):
        zjcj = self._compute_zj_cj()
        # entering: most negative
        candidates = [(v,i) for i,v in enumerate(zjcj[:-1]) if v < -self.tol]
        if not candidates:
            return None
        val, col = min(candidates)
        # leaving by min ratio
        ratios = []
        for i, row in enumerate(self.tableau[:-1]):
            a = row[col]
            if a > self.tol:
                ratios.append((row[-1] / a, i))
        if not ratios:
            raise ValueError('Unbounded solution')
        _, row = min(ratios)
        return row, col, self.tableau[row][col]

    def _pivot(self, r, c):
        pv = self.tableau[r][c]
        self.tableau[r] = [v / pv for v in self.tableau[r]]
        for i in range(len(self.tableau)):
            if i == r:
                continue
            fac = self.tableau[i][c]
            self.tableau[i] = [self.tableau[i][j] - fac * self.tableau[r][j]
                                for j in range(len(self.tableau[0]))]
        self.basic_vars[r] = self.var_names[c]

    def _log_tableau(self, pivot_info=None):
        col_width = 12
        header = ['BV'] + self.var_names + ['RHS']
        # Header line
        lines = [''.join(f"{h:>{col_width}}" for h in header)]
        # Constraint rows
        for i, row in enumerate(self.tableau[:-1]):
            values = row + []  # include RHS at end
            formatted = [self._fmt(v) for v in values]
            line = f"{self.basic_vars[i]:>{col_width}}" + ''.join(f"{val:>{col_width}}" for val in formatted)
            lines.append(line)
        # Zj-Cj row
        zjcj = self._compute_zj_cj()
        formatted_z = [self._fmt(v) for v in zjcj]
        zj_line = f"{'ZJ-CJ':>{col_width}}" + ''.join(f"{val:>{col_width}}" for val in formatted_z)
        lines.append(zj_line)
        # Append to log
        self.log += "\n".join(lines) + "\n"
        if pivot_info:
            pivot_str = f"Pivot:{self._fmt(pivot_info['value'])}, enter={pivot_info['enter']}, leave={pivot_info['leave']}"
            self.log += f"{pivot_str:^{col_width * len(header)}}\n"

    def solve(self):
        it = 0
        try:
            while True:
                p = self._find_pivot()
                if p is None:
                    break
                r, c, val = p
                leave = self.basic_vars[r]
                enter = self.var_names[c]
                self._pivot(r, c)
                it += 1
                self._log_tableau({'value': val, 'enter': enter, 'leave': leave})
            # infeasible check
            for i, bv in enumerate(self.basic_vars):
                if bv.startswith('a') and abs(self.tableau[i][-1]) > self.tol:
                    self.log += "Infeasible solution.\n"
                    return "Infeasible solution"
            # extract
            sol = [0]*len(self.orig_objective)
            for i, bv in enumerate(self.basic_vars):
                if bv.startswith('x'):
                    idx = int(bv[1:]) - 1
                    sol[idx] = self.tableau[i][-1]
            z_internal = self.tableau[-1][-1]
            # convert back for minimize
            z = self.sign * z_internal
            # log final
            self.log += "Final Solution:\n"
            for i, v in enumerate(sol):
                self.log += f"x{i+1} = {self._fmt(v)}\n"
            self.log += f"Objective = {self._fmt(z)}\n"
            basis = [self.var_names.index(b) for b in self.basic_vars]
            non_basis = [i for i in range(len(self.var_names)) if i not in basis]
            return {
                'solution': sol,
                'objective_value': z,
                'tableau': self.tableau,
                'basis_variables': basis,
                'non_basis_variables': non_basis,
                'iterations': it,
                'log': self.log
            }
        except ValueError as e:
            if 'Unbounded' in str(e):
                self.log += 'Unbounded solution.\n'
                return 'Unbounded solution'
            raise


class SimplexSolverGUI:
    def __init__(self, master):
        self.master = master
        master.title("Simplex Solver - Simplex Method")
        master.geometry("1920x1080")
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

        self.result_text = tk.Text(master, height=25, )
        self.result_text.pack(fill=tk.BOTH, expand=True)

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
