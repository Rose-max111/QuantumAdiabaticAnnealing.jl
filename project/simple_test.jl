using CairoMakie
function rule(cell, left, right)
    # 示例规则：模仿Rule 110
    return (left | right) ⊻ cell
end

# 初始化元胞自动机
function initialize_automaton(size)
    return rand(0:1, size)  # 随机初始化元胞状态
end

# 更新元胞自动机
function update_automaton(cells)
    new_cells = zeros(Int,length(cells))
    size = length(cells)
    
    # 更新奇数位置
    for i in 1:2:size
        left = cells[mod1(i-1, size)]
        right = cells[mod1(i+1, size)]
        new_cells[i] = rule(cells[i], left, right)
    end
    
    # 更新偶数位置
    for i in 2:2:size
        left = new_cells[mod1(i-1, size)]
        right = new_cells[mod1(i+1, size)]
        new_cells[i] = rule(cells[i], left, right)
    end

    return new_cells
end

# 示例运行
size = 30
steps = 200
cells = initialize_automaton(size)
println("Initial state: ", cells)

history = []
for step in 1:steps
    cells = update_automaton(cells)
    push!(history, cells)
    # println("Step $step: ", cells)
end

f = Figure()
ax = Axis(f[1,1])
centers_x = 1:steps+1
centers_y = 1:size+1
heatmap!(ax, hcat(history...), colormap = :grays)
f