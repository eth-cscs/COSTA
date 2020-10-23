<p align="center"><img src="./docs/costa-logo.svg" width="50%"></p>

COSTA is a communication-optimal, highly-optimised algorithm for data redistribution accross multiple processors, using `MPI` and `OpenMP` and offering the possibility to transpose and scale some or all data. It implements scalapack routines for matrix scale & transpose operations (`sub(C) = alpha * sub(A)^T + beta * C`, provided by `pxtran(u)`) and data redistribution (`sub(C) = sub(A)`, provided by `pxgemr2d`) and outperforms other scalapack implementations by orders of magnitude in some cases. Unlike previous redistribution algorithms, COSTA will also propose the relabelling of MPI ranks that minimizes the data reshuffling cost, leaving to users to decide if they want to use it. This way, if the initial and the target data distributions differ up to a rank permutation, COSTA will perform no communication, whereas other algorithms will reshuffle all the data. Thanks to its optimizations, significant speedups will be achieved even if the proposed rank relabelling is not used.

What makes COSTA more general than scalapack routines is that it is not limited only to block-cyclic data distributions, but can deal with completely arbitrary and irregular matrix distributions and can be easily generalized for n-dimensional tensors. 

Thanks to its scalapack wrappers, scalapack users do not need to change their code in order to use COSTA: it is enough to link your library to COSTA before linking to scalapack and all `pxtran, pxtranu` and `pxgemr2d` routines will automatically be using the COSTA algorithm.

## Examples

### Block-cyclic (Scalapack) Layout

To represent an arbitrary block-cyclic (scalapack) layout, we can use the following function defined in `costa/layout.hpp` header:
```cpp
#include <costa/layout.hpp>
// ...
template <typename T>
grid_layout<T> block_cyclic_layout<double>(
                   const int m, const int n,             // global matrix dimensions
                   const int block_m, const int block_n, // block dimensions
                   const int i, const int j,             // submatrix start
                                                         // (1-based, scalapack-compatible)
                   const int sub_m, const int sub_n,     // submatrix size
                   const int p_m, const int p_n,         // processor grid dimension
                   const char order,                     // rank grid ordering ('R' or 'C')
                   const int rsrc, const int csrc,       // coordinates of ranks oweing 
                                                         // the first row (0-based)
                   T* ptr,                               // local data of matrix A 
                                                         // (not the submatrix)
                   const int lld,                        // local leading dimension
                   const int rank                        // processor rank
               );
```
The arguments meaning can be nicely visualized with the following figure, where the red submatrix is represented:
<p align="center"><img src="./docs/block-cyclic.svg" width="100%"></p>

In case we want to represent the full matrix (instead of a submatrix), it suffices to put:
```cpp
// start of the submatrix is the start of the full matrix
int i = 1; int j = 1 // 1-based due to scalapack-compatibility
// size of the submatrix is that size of the full matrix
int sub_m = m; int sub_n = n
```

For a complete example of transforming between two block-cyclic matrix layouts, please refer to `examples/example0.cpp`.


## Miniapps (for testing and benchmarking)

### Data-redistribution with pxgemr2d

COSTA implements ScaLAPACK `pxgemr2d` routines that transforms the matrix between two block-cyclic data layouts (`sub(C) = sub(A)`) where the two matrices do not necessarily have to belong to the same MPI communicators. In addition, COSTA will propose the MPI rank relabelling that minimizes the data reshuffling cost and that user is free to choose whether to use it. 

The miniapp consists of an executable `./build/examples/pxgemr2d_miniapp` which can be run as follows (assuming we are in the root folder of the project):

```bash
# set the number of threads to be used by each MPI rank
export OMP_NUM_THREADS=18
# if using CPU version with MKL backend, set MKL_NUM_THREADS as well
export MKL_NUM_THREADS=18 
# run the miniapp
mpirun -np 4 ./build/examples/pxgemr2d_miniapp -m 1000 -n 1000 \
                                            --block_a=128,128 \ 
                                            --block_c=128,128 \
                                            --p_grid=2,2 \
                                            --type=double \
                                            --algorithm=costa
```

The overview of all supported options is given below:
- `-m (--m_dim)` (default: `1000`): number of rows of matrices `A` and `C`.
- `-n (--n_dim)` (default: `1000`): number of columns of matrices `A` and `C`. 
- `--block_a` (optional, default: `128,128`): 2D-block size for matrix A. 
- `--block_c` (optional, default `128,128`): 2D-block size for matrix C.
- `-p (--p_grid)` (optional, default: `1,P`): 2D-processor grid. By default `1xP` where `P` is the total number of MPI ranks.
- `-r (--n_rep)` (optional, default: 2): number of repetitions.
- `-t (--type)` (optional, default: `double`): data type of matrix entries. Can be one of: `float`, `double`, `zfloat` and `zdouble`. The last two correspond to complex numbers.
- `--test` (optional): if present, the result of COSTA will be verified with the result of the available SCALAPACK.
- `--algorithm` (optional, default: `both`): defines which algorithm (`costa`, `scalapack` or `both`) to run.
- `-h (--help) (optional)`: print available options.

### Scale and Transpose with pxtran and pxtranu

COSTA implements ScaLAPACK `pxtran` and `pxtranu` routines that performs the scale and transpose operation, given by:
```sub(C) = alpha * sub(A)^T + beta * sub(C)```
In addition, COSTA will propose the MPI rank relabelling that minimizes the data reshuffling cost and that user is free to choose whether to use it. 

The miniapp consists of an executable `./build/examples/pxtran_miniapp` which can be run as follows (assuming we are in the root folder of the project):

```bash
# set the number of threads to be used by each MPI rank
export OMP_NUM_THREADS=18
# if using CPU version with MKL backend, set MKL_NUM_THREADS as well
export MKL_NUM_THREADS=18 
# run the miniapp
mpirun -np 4 ./build/examples/pxtran_miniapp -m 1000 -n 1000 -k 1000 \
                                            --block_a=128,128 \ 
                                            --block_c=128,128 \
                                            --p_grid=2,2 \
                                            --alpha=1 \
                                            --beta=1 \
                                            --type=double \
                                            --algorithm=costa
```

The overview of all supported options is given below:
- `-m (--m_dim)` (default: `1000`): number of rows of matrices `A` and `C`.
- `-n (--n_dim)` (default: `1000`): number of columns of matrices `A` and `C`. 
- `--block_a` (optional, default: `128,128`): 2D-block size for matrix A. 
- `--block_c` (optional, default `128,128`): 2D-block size for matrix C.
- `-p (--p_grid)` (optional, default: `1,P`): 2D-processor grid. By default `1xP` where `P` is the total number of MPI ranks.
- `--alpha` (optional, default: 1): alpha parameter in `sub(C) = alpha*sub(A)^T + beta*sub(C)`.
- `--beta` (optional, default: 0): beta parameter in `sub(C) = alpha*sub(A)^T + beta*sub(C)`.
- `-r (--n_rep)` (optional, default: 2): number of repetitions.
- `-t (--type)` (optional, default: `double`): data type of matrix entries. Can be one of: `float`, `double`, `zfloat` and `zdouble`. The last two correspond to complex numbers.
- `--test` (optional): if present, the result of COSTA will be verified with the result of the available SCALAPACK.
- `--algorithm` (optional, default: `both`): defines which algorithm (`costa`, `scalapack` or `both`) to run.
- `-h (--help) (optional)`: print available options.
