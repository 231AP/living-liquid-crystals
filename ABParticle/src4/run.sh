#! /bin/bash
#!/bin/bash

# 定义每个变量的值，可以根据需要进行修改

V0=10
ABParticle=1


celllist_size_x=1
celllist_size_y=1
celllist_count_x=128
celllist_count_y=128

celllist_size_x1=4
celllist_size_y1=4
celllist_count_x1=32
celllist_count_y1=32
box_x_length=$((celllist_size_x*celllist_count_x))
box_y_length=$((celllist_size_y*celllist_count_y))
echo $box_x_length

particle_density=4
max_particles_per_cell=200
max_particles_per_cell1=200
neighbor_distance=3
mask0=5
mask1=3
min_distance=0.1
equilibrium_distance=1
force_coefficient=3
kBT=0.2
viscosity_coefficient=1
neighbor_update_threshold=1
neighbor_update_threshold1=1
total_particles=6000
start_time=0
end_time=200
time_step=0.001
tExpo=0.1

dir_name="boxX${box_x_length}_boxY${box_y_lenth}_particles${total_particles}_endTime${end_time}_kBT${kBT}"
echo $dir_name
# 创建并写入 input.dat 文件
cat > input.dat << EOF
# box X length
$box_x_length
# box Y length
$box_y_length
# x方向celllist大小
$celllist_size_x
# y方向celllist大小
$celllist_size_y
# x方向celllist个数
$celllist_count_x
# y方向celllist个数
$celllist_count_y
# 粒子密度
$particle_density
# 每个cell里粒子最大个数
$max_particles_per_cell
# neighbor距离
$neighbor_distance
# mask0
$mask0
# mask1
$mask1
# 每个粒子间的最小距离，用于位运算
$min_distance
# 粒子平衡距离r0
$equilibrium_distance
# 作用力的系数
$force_coefficient
# kBT
$kBT
# 粘滞系数
$viscosity_coefficient
# 近邻表更新临界位移
$neighbor_update_threshold
# 粒子总数
$total_particles
# 开始时间
$start_time
# 结束时间
$end_time
# 时间步长
$time_step
# tExpo
$tExpo
# useField
$V0
# ABParticle
$ABParticle

# x方向celllist大小
$celllist_size_x1
# y方向celllist大小
$celllist_size_y1
# x方向celllist个数
$celllist_count_x1
# y方向celllist个数
$celllist_count_y1
# 每个cell里粒子最大个数
$max_particles_per_cell1
# 近邻表更新临界位移
$neighbor_update_threshold1

EOF

echo "input.dat 文件已生成。"


# 初始目录名称
#dir_name="data"
count=0

#循环查找不存在的目录名称
while [ -d "$dir_name" ]; do
    count=$((count + 1))
    dir_name="data$count"
done

# mkdir $dir_name
# cp kernel.cu $dir_name
# cp input.dat $dir_name
# cd $dir_name

