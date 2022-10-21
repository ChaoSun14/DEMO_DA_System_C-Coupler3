import sys
import os
import re
 
root_path = "./results/"
#grid_info = "10km"
#proc = "1000_10"
grid_info = sys.argv[1]
proc = sys.argv[2]
file_name = "/ccpl_output_time.out"
#pattern = re.compile(r"(\d{1-2}\.\d{6})")
pattern = re.compile(r"[(](.*?)[)]")

total_t1 = 0.0
total_t2 = 0.0
total_t3 = 0.0
number1 = 0
number2 = 0
number3 = 0
result = []
nn = 0
with open(root_path+grid_info+"/"+proc+file_name) as f:
    lines = f.readlines()
    line1 = ""
    line2 = ""
    line3 = ""  
    for line in lines:
        if line.find("Ensemble_procedures_inst::run: individual run copy in") != -1:
            t1 = float(pattern.findall(line)[-1])
            if t1 > 0:
                number1 = number1 + 1
                total_t1 = total_t1 + t1
        if line.find("Ensemble_procedures_inst::run: individual run copy out") != -1:
            t2 = float(pattern.findall(line)[-1])
            if t2 > 0:
                total_t1 = total_t1 + t2
        if line.find("TIME in External_procedures_inst before run do ocn_da_demo_individual run") != -1:
            t3 = float(pattern.findall(line)[-1])
            if t3 > 0: total_t1 = total_t1 + t3
        if line.find("TIME in External_procedures_inst after run do ocn_da_demo_individual run") != -1:
            t4 = float(pattern.findall(line)[-1])
            if t4 > 0: total_t1 = total_t1 + t4
        if line.find("TIME in Ensemble_procedures_inst::run: ensemble run before run do ocn_da_demo_ensmean run") != -1:
            t5 = float(pattern.findall(line)[-1])
            if t5 > 0: 
                total_t2 = total_t2 + t5
        if line.find("TIME in Ensemble_procedures_inst::run: ensemble run after run do ocn_da_demo_ensmean run") != -1:
            t6 = float(pattern.findall(line)[-1])
            if t6 > 0: 
                total_t2 = total_t2 + t6
        if line.find("TIME in External_procedures_inst before run do ocn_da_demo_ensmean run") != -1:
            t7 = float(pattern.findall(line)[-1])
            if t7 > 0: 
                number2 = number2 +1
                total_t2 = total_t2 + t7
        if line.find("TIME in External_procedures_inst after run do ocn_da_demo_ensmean run") != -1:
            t8 = float(pattern.findall(line)[-1])
            if t8 > 0: 
                total_t1 = total_t2 + t8
        if line.find("TIME in Ensemble_procedures_inst::run: ensemble run before run do ocn_da_demo_ensgather run") != -1:
            t9 = float(pattern.findall(line)[-1])
            if t9 > 0: 
                total_t3 = total_t3 + t9
        if line.find("TIME in Ensemble_procedures_inst::run: ensemble run after run do ocn_da_demo_ensgather run") != -1:
            t10 = float(pattern.findall(line)[-1])
            if t10 > 0: 
                total_t3 = total_t3 + t10
        if line.find("TIME in External_procedures_inst before run do ocn_da_demo_ensgather run") != -1:
            t11 = float(pattern.findall(line)[-1])
            if t11 > 0: total_t3 = total_t3 + t11
        if line.find("TIME in External_procedures_inst after run do ocn_da_demo_ensgather run") != -1:
            t12 = float(pattern.findall(line)[-1])
            if t12 > 0: 
                number3 = number3 +1
                total_t3 = total_t3 + t12
#print("number", number1, number2, number3)
total_t1 = total_t1/number1
total_t2 = total_t2/number2
total_t3 = total_t3/number3
#result.append(number1)
result.append('%.6f'%total_t1)
#result.append(number2)
result.append('%.6f'%total_t2)
#result.append(number3)
result.append('%.6f'%total_t3)
print("CCPL_run", result)
