using CUDA

function gpu_add!(y, x)
    id = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    
    Nx, Ny = size(y)
    cind = CartesianIndices((Nx, Ny))

    # println(blockIdx().x, " ", threadIdx().x)
    for k in id:stride:Nx*Ny
        i = cind[k][1]
        j = cind[k][2]
        y[i,j] += x[i,j]
    end
    return nothing
end

N = 2^14
M = 2^14
x_d = CUDA.fill(1.0f0, (N, M))
y_d = CUDA.fill(2.0f0, (N, M))

kernel = @cuda launch=false gpu_add!(y_d, x_d)
config = launch_configuration(kernel.fun)
threads = min(N*M, config.threads)
blocks = cld(N*M, threads)

@time kernel(y_d, x_d; threads, blocks)

@time y_d .+= x_d