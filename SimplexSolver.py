class SimplexTable:
    def __init__(self, matrix, column_vecs):
        self.matrix = matrix
        self.column_vecs = column_vecs

        if not (type(self.matrix) == list and all(type(row) == list for row in self.matrix)):
            raise ValueError("Needs to be a 2D list")

        if not all(all((type(value) == int or type(value) == float) for value in row) for row in self.matrix):
            raise ValueError("All elements in the matrix must be int or float")

        num_rows = len(self.matrix)
        if not all(len(row) == num_rows for row in self.matrix):
            raise ValueError("Matrix must be square (n x n)")

        if not (type(self.column_vecs) == list and all(type(col) == list for col in self.column_vecs)):
            raise ValueError("Column vectors must be a list of lists")

        if not all(len(col) == num_rows for col in self.column_vecs):
            raise ValueError("Each column vector must have the same number of rows as the matrix")

    def display_tableau(self):

        for i in range(len(self.matrix)):
            print(self.matrix[i] + [col[i] for col in self.column_vecs])

    def tableau(self):

        result = []
        for i in range(len(self.matrix)):

            result.append(self.matrix[i] + [col[i] for col in self.column_vecs])
        return result

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
            
            
#Testing Code
matrix = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
column_vecs = [[10, 11, 12], [13, 14, 15]]
table = SimplexTable(matrix, column_vecs)

table.display_tableau()
extended_tableau = table.tableau()
print("Tableau as a 2D list:")
print(extended_tableau)
table.scale_row(0, 2)
print("After Scaling:")
table.display_tableau()
table.subtract_row_multiple(1, 0, 2)
print("After row operations:")
table.display_tableau()
