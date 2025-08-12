# Type Hints Guide for CALFEM Python

## Overview
Type hints have been added to CALFEM Python to improve code documentation, IDE support, and catch potential type-related bugs early. This is done in a backward-compatible way that doesn't affect runtime behavior.

## Benefits
1. **Better IDE Support**: IDEs can provide better autocomplete, error detection, and refactoring
2. **Documentation**: Function signatures are self-documenting 
3. **Error Prevention**: Static type checkers can catch type-related bugs before runtime
4. **Backward Compatibility**: All existing code continues to work unchanged

## Type Hint Patterns

### Basic Types
```python
def func(param: int) -> str:
    return str(param)
```

### Array Types
```python
from numpy.typing import ArrayLike, NDArray
import numpy as np

# Input arrays (can be lists, tuples, numpy arrays)
def func(arr: ArrayLike) -> NDArray[np.floating]:
    return np.array(arr, dtype=float)
```

### Optional Parameters
```python
from typing import Optional

def func(required: ArrayLike, optional: Optional[ArrayLike] = None) -> NDArray[np.floating]:
    if optional is None:
        return np.array(required)
    return np.array(required) + np.array(optional)
```

### Functions with Conditional Returns
```python
from typing import Union, Tuple

# Returns either Ke alone or (Ke, fe) tuple
def element_func(ep: ArrayLike, eq: Optional[ArrayLike] = None) -> Union[NDArray[np.floating], Tuple[NDArray[np.floating], NDArray[np.floating]]]:
    Ke = compute_stiffness(ep)
    if eq is None:
        return Ke
    else:
        fe = compute_load(eq)
        return Ke, fe
```

### Sparse Matrix Support
```python
from scipy.sparse import csr_matrix, csc_matrix, lil_matrix

def assem_func(K: Union[NDArray[np.floating], csr_matrix, csc_matrix, lil_matrix], 
               Ke: ArrayLike) -> Union[NDArray[np.floating], csr_matrix, csc_matrix, lil_matrix]:
    # Function body
    pass
```

## Common CALFEM Type Patterns

### Element Stiffness Functions
- Input: `ex: ArrayLike, ey: ArrayLike, ep: ArrayLike, eq: Optional[ArrayLike] = None`
- Return: `Union[NDArray[np.floating], Tuple[NDArray[np.floating], NDArray[np.floating]]]`

### Element Stress/Force Functions  
- Input: `ex: ArrayLike, ey: ArrayLike, ep: ArrayLike, ed: ArrayLike`
- Return: `NDArray[np.floating]` or specific types like `float`

### Utility Functions
- Coordinate functions: `(int, int) -> NDArray[np.integer]`
- Warning/Error functions: `str -> None`

## Type Checking
Install mypy for static type checking:
```bash
pip install mypy
```

Run type checking:
```bash
mypy src/calfem/core.py
```

The mypy.ini configuration allows gradual typing adoption without breaking existing code.

## Implementation Strategy
1. Start with the most commonly used functions
2. Add type hints to new functions as they're written
3. Gradually add hints to existing functions during maintenance
4. Focus on public API functions first
5. Use `# type: ignore` comments for complex cases that mypy can't handle

## Examples in Core Module
The following functions already have type hints as examples:
- `spring1e`: Simple function with array input/output
- `spring1s`: Function returning scalar
- `bar1e`: Function with optional parameters and conditional return
- `beam2e`: Complex element function pattern
- `eigen`: Functions with multiple array inputs/outputs
- `assem`: Functions supporting both dense and sparse matrices
- `create_dofs`: Utility function with integer types

## Future Enhancements
- Add type hints to other modules (mesh, vis, etc.)
- Use Protocol types for more complex interfaces
- Add runtime type checking with libraries like pydantic
- Integration with documentation generation tools

## Compatibility
- Minimum Python version: 3.8 (supports all required typing features)
- No runtime impact: Type hints are ignored at runtime
- All existing code continues to work without modification
- Optional adoption: Teams can choose their level of type hint usage
