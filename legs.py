import re
import math

input_file = "legs.txt"
output_data_list = []

def parse_legs_file(filename):
    data = []
    try:
        with open(filename, "r", encoding='utf-8') as f:
            content = f.read()
    except UnicodeDecodeError:
        with open(filename, "r", encoding='gbk') as f:
            content = f.read()
            
    # Split content into blocks based on the header "【装配体】"
    blocks = content.split("【装配体】")
    
    for block in blocks:
        if not block.strip():
            continue
            
        # Extract Y and Z from COM section
        idx_com = block.find("质心")
        if idx_com == -1:
            continue
            
        sub_com = block[idx_com:]
        idx_end_com = sub_com.find("惯性主轴")
        if idx_end_com != -1:
            sub_com = sub_com[:idx_end_com]
            
        y_match = re.search(r"Y\s*=\s*([-\d\.]+)", sub_com)
        z_match = re.search(r"Z\s*=\s*([-\d\.]+)", sub_com)
        
        if not (y_match and z_match):
            continue
            
        y_val = float(y_match.group(1))
        z_val = float(z_match.group(1))
        
        # Extract Ixx from Inertia Tensor section (Taken at Output Coordinate System)
        ixx_match = re.search(r"惯性张量:.*?\n\s*由输出座标系决定.*?\n.*?Ixx\s*=\s*([\-\d\.]+)", block, re.DOTALL)
        
        if not ixx_match:
            continue
        
        ixx_val = float(ixx_match.group(1))
        
        data.append((y_val, z_val, ixx_val))
        
    return data

parsed_data = parse_legs_file(input_file)
print(f"Parsed {len(parsed_data)} data points.")
print("% 格式: [腿长(m), 质心偏离角度(rad), 质心到髋关节距离(m), 转动惯量(kg·m²)]")

start_leg_length = 0.1
increment = 0.01

for i, (y_val, z_val, ixx_val) in enumerate(parsed_data):
    leg_length = start_leg_length + i * increment
    y_m = y_val / 1000.0
    z_m = z_val / 1000.0
    ixx_kgm2 = ixx_val / 1e9
    # 计算重心到髋关节距离
    cg_dist = (y_m**2 + z_m**2) ** 0.5
    # 计算重心偏移Y轴负方向的角度（弧度）
    angle_rad = math.atan2(-z_m, -y_m)
    output_data_list.append([leg_length, angle_rad, cg_dist, ixx_kgm2])

for row in output_data_list:
    print("{};".format(", ".join(map(str, row))))