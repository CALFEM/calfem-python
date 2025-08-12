import numpy as np
from numpy.linalg import inv
import sys

class MatrixCompat:
    """
    Compatibility layer to replace np.matrix in calfem-python for NumPy 2.0 support.
    
    This class mimics the behavior of np.matrix, particularly:
    - Matrix multiplication with * operator
    - .I property for matrix inverse
    - .T property for matrix transpose
    - Array-like indexing and shape behavior
    """
    
    def __init__(self, input_array):
        """
        Initialize a MatrixCompat object
        
        Parameters:
            input_array: Can be a MatrixCompat instance, numpy array, list, or other compatible type
        """
        if isinstance(input_array, MatrixCompat):
            self.array = input_array.array.copy()
        elif isinstance(input_array, np.ndarray):
            self.array = input_array.copy()
        else:
            self.array = np.array(input_array, dtype=float)
        
        # Ensure that 1D inputs mimic np.matrix behavior: create a 1 x n row matrix
        # Original np.matrix(list) yields shape (1, n). Code elsewhere (e.g., np.matrix(vec).T)
        # depends on this to obtain an (n,1) column after transpose. Previous implementation
        # reshaped to (n,1), breaking expressions like np.matrix(q).T used to build column vectors
        # and causing dimension mismatches (see apply_traction_linear_element). Restore semantics.
        if self.array.ndim == 1:
            self.array = self.array.reshape(1, -1)
    
    # Negation
    def __neg__(self):
        """
        Implement unary negation -A
        """
        return MatrixCompat(-self.array)
    
    # Matrix multiplication
    def __mul__(self, other):
        """
        Implement matrix multiplication A * B
        """
        if isinstance(other, (int, float)):
            # Scalar multiplication
            return MatrixCompat(self.array * other)
        elif isinstance(other, MatrixCompat):
            # Matrix multiplication
            return MatrixCompat(np.matmul(self.array, other.array))
        elif isinstance(other, np.ndarray):
            # Matrix multiplication with numpy array
            if other.ndim == 1:
                # Handle 1D arrays specially
                return MatrixCompat(np.matmul(self.array, other.reshape(-1, 1)))
            return MatrixCompat(np.matmul(self.array, other))
        else:
            # Try to convert to array and multiply
            try:
                other_array = np.array(other, dtype=float)
                if other_array.ndim == 1:
                    other_array = other_array.reshape(-1, 1)
                return MatrixCompat(np.matmul(self.array, other_array))
            except:
                return NotImplemented
    
    def __rmul__(self, other):
        """
        Implement right multiplication: B * A
        """
        if isinstance(other, (int, float)):
            # Scalar multiplication
            return MatrixCompat(self.array * other)
        else:
            # Try to convert to array and multiply
            try:
                other_array = np.array(other, dtype=float)
                if other_array.ndim == 1:
                    other_array = other_array.reshape(1, -1)
                return MatrixCompat(np.matmul(other_array, self.array))
            except:
                return NotImplemented
    
    # Division operations
    def __truediv__(self, other):
        """
        Implement division A / B
        """
        if isinstance(other, (int, float)):
            # Scalar division
            return MatrixCompat(self.array / other)
        elif isinstance(other, MatrixCompat):
            # Matrix division (A / B = A * B^-1)
            return MatrixCompat(np.matmul(self.array, inv(other.array)))
        else:
            try:
                # Try to treat as scalar
                return MatrixCompat(self.array / other)
            except:
                return NotImplemented
    
    def __rtruediv__(self, other):
        """
        Implement right division: B / A
        """
        if isinstance(other, (int, float)):
            # Scalar / Matrix
            return MatrixCompat(other / self.array)
        else:
            try:
                other_array = np.array(other, dtype=float)
                # B / A = B * A^-1
                return MatrixCompat(np.matmul(other_array, inv(self.array)))
            except:
                return NotImplemented
    
    # Handle inverse property
    @property
    def I(self):
        """
        Matrix inverse property, similar to np.matrix.I
        """
        return MatrixCompat(inv(self.array))
    
    # Handle transpose property
    @property
    def T(self):
        """
        Matrix transpose property, similar to np.matrix.T
        """
        return MatrixCompat(self.array.T)
    
    # Support array-like indexing
    def __getitem__(self, idx):
        """
        Support array indexing: A[i, j] or A[i]
        """
        result = self.array[idx]
        if isinstance(result, np.ndarray) and result.ndim > 0:
            return MatrixCompat(result)
        return result
    
    def __setitem__(self, idx, value):
        """
        Support array assignment: A[i, j] = value
        """
        self.array[idx] = value
    
    # Additional utility methods
    @property
    def shape(self):
        """
        Return the shape of the array
        """
        return self.array.shape
    
    def toarray(self):
        """
        Convert to a regular numpy array
        """
        return self.array.copy()
    
    def __repr__(self):
        """
        String representation
        """
        return f"MatrixCompat({self.array.__repr__()})"
    
    def __str__(self):
        """
        String conversion
        """
        return str(self.array)
        
    # Support matrix operations with @ operator
    def __matmul__(self, other):
        """
        Support for @ operator (Python 3.5+)
        """
        if isinstance(other, MatrixCompat):
            return MatrixCompat(np.matmul(self.array, other.array))
        else:
            return MatrixCompat(np.matmul(self.array, other))
            
    def __rmatmul__(self, other):
        """
        Support for @ operator (Python 3.5+) on right side
        """
        return MatrixCompat(np.matmul(other, self.array))
        
    # Addition operation
    def __add__(self, other):
        """
        Addition: A + B
        """
        if isinstance(other, MatrixCompat):
            return MatrixCompat(self.array + other.array)
        else:
            return MatrixCompat(self.array + other)
            
    def __radd__(self, other):
        """
        Right addition: B + A
        """
        return MatrixCompat(other + self.array)
    
    # Subtraction operation
    def __sub__(self, other):
        """
        Subtraction: A - B
        """
        if isinstance(other, MatrixCompat):
            return MatrixCompat(self.array - other.array)
        else:
            return MatrixCompat(self.array - other)
            
    def __rsub__(self, other):
        """
        Right subtraction: B - A
        """
        return MatrixCompat(other - self.array)
    
    # Support for numpy functions
    def __array__(self, dtype=None):
        """
        Convert to numpy array when used in numpy functions
        """
        if dtype is None:
            return self.array
        else:
            return self.array.astype(dtype)
    
    # Additional matrix operations
    def dot(self, other):
        """
        Matrix dot product
        """
        if isinstance(other, MatrixCompat):
            return MatrixCompat(np.dot(self.array, other.array))
        else:
            return MatrixCompat(np.dot(self.array, other))
    
    def copy(self):
        """
        Return a copy of the matrix
        """
        return MatrixCompat(self.array.copy())
    
    # Implement common numpy methods
    def reshape(self, *args, **kwargs):
        """
        Reshape the array
        """
        return MatrixCompat(self.array.reshape(*args, **kwargs))
    
    def sum(self, *args, **kwargs):
        """
        Sum the array elements
        """
        result = self.array.sum(*args, **kwargs)
        if isinstance(result, np.ndarray):
            return MatrixCompat(result)
        return result
    
    # Support conversion to and from list
    def tolist(self):
        """
        Convert to list
        """
        return self.array.tolist()
    
    # Support hstack/vstack compatibility
    def __len__(self):
        """
        Return the length of the array
        """
        return self.array.shape[0]
        
    # Add flatten method
    def flatten(self):
        """
        Return a flattened copy of the array
        """
        return self.array.flatten()

# Now we patch numpy to use our compatibility layer, but only for NumPy >= 2.0
# where np.matrix has been removed

def _get_numpy_version():
    """Get numpy version as tuple of integers"""
    version_str = np.__version__.split('.')
    return tuple(int(x) for x in version_str[:2])

# Only patch if NumPy version is 2.0 or higher
if _get_numpy_version() >= (2, 0):
    # Store original if it exists (for potential restoration)
    np_matrix_original = getattr(np, 'matrix', None)
    np.matrix = MatrixCompat
    np.mat = MatrixCompat  # Alias for np.matrix
else:
    # For NumPy < 2.0, keep the original matrix
    np_matrix_original = np.matrix

# Optional cleanup function to restore original np.matrix if needed
def restore_numpy_matrix():
    """
    Restore the original np.matrix function
    """
    np.matrix = np_matrix_original
    np.mat = np_matrix_original  # Restore alias as well