import sys
import os
import re
 
root_path = "./results/"
#grid_info="10km"
#proc = "1000_10"
grid_info = sys.argv[1]
proc = sys.argv[2]
file_name = "/model_output_time.out"
pattern = re.compile(r"(\d{1-2}\.\d{6})")
#pattern = re.compile(r"(\d+(.\d+)?)")

total_da_t01 = 0.0
total_da_t02 = 0.0
total_da_t03 = 0.0
total_da_t04 = 0.0
total_damean_t01 = 0.0
total_damean_t02 = 0.0
total_damean_t03 = 0.0
total_damean_t04 = 0.0
total_dagather_t01 = 0.0
total_dagather_t02 = 0.0
total_dagather_t03 = 0.0
total_dagather_t04 = 0.0
number_da_01 = 0
number_da_02 = 0
number_damean_01 = 0
number_damean_02 = 0
number_dagather_01 = 0
number_dagather_02 = 0
total_t1 = 0.0
total_t2 = 0.0
total_t3 = 0.0
number1 = 0
number2 = 0
number3 = 0
with open(root_path+grid_info+"/"+proc+file_name) as f:
    lines = f.readlines()
    line1 = ""
    line2 = ""
    line3 = ""	
    for line in lines:
        if line.find("[CCPL <da_demo>] ccpl time all") != -1:
            line1 = line.split(':')
            t1 = float(line1[-1])
            number_da_01 = number_da_01 + 1
            total_da_t01 = total_da_t01 + t1
        if line.find("[CCPL <da_demo>] ccpl time grids") != -1:
            line1 = line.split(':')
            t2 = float(line1[-1])
            total_da_t02 = total_da_t02 + t2
        if line.find("[CCPL <da_demo>] ccpl time fields") != -1:
            line1 = line.split(':')
            t3 = float(line1[-1])
            total_da_t03 = total_da_t03 + t3
        if line.find("[CCPL <da_demo>] total time") != -1:
            line1 = line.split(':')
            t4 = float(line1[-1])
            number_da_02 = number_da_02 + 1
            total_da_t04 = total_da_t04 + t4
        if line.find("individual total time") != -1:
            line1 = line.split(':')
            t1 = float(line1[-1])
            number1 = number1 + 1
            total_t1 = total_t1 + t1
        if line.find("[CCPL <da_ensmean_demo>] ccpl time all") != -1:
            line2 = line.split(':')
            t1 = float(line2[-1])
            number_damean_01 = number_damean_01 + 1
            total_damean_t01 = total_damean_t01 + t1
        if line.find("[CCPL <da_ensmean_demo>] ccpl time grids") != -1:
            line2 = line.split(':')
            t2 = float(line2[-1])
            total_damean_t02 = total_damean_t02 + t2
        if line.find("[CCPL <da_ensmean_demo>] ccpl time fields") != -1:
            line2 = line.split(':')
            t3 = float(line2[-1])
            total_damean_t03 = total_damean_t03 + t3
        if line.find("[CCPL <da_ensmean_demo>] total time") != -1:
            line2 = line.split(':')
            t4 = float(line2[-1])
            number_damean_02 = number_damean_02 + 1
            total_damean_t04 = total_damean_t04 + t4
        if line.find("ensmean total time") != -1:
        	line2 = line.split(':')
        	t2 = float(line2[-1])
        	number2 = number2 + 1
        	total_t2 = total_t2 + t2
        if line.find("[CCPL <da_ensgather_demo>] ccpl time all") != -1:
            line3 = line.split(':')
            t1 = float(line3[-1])
            number_dagather_01 = number_dagather_01 + 1
            total_dagather_t01 = total_dagather_t01 + t1
        if line.find("[CCPL <da_ensgather_demo>] ccpl time grids") != -1:
            line3 = line.split(':')
            t2 = float(line3[-1])
            total_dagather_t02 = total_dagather_t02 + t2
        if line.find("[CCPL <da_ensgather_demo>] ccpl time fields") != -1:
            line3 = line.split(':')
            t3 = float(line3[-1])
            total_dagather_t03 = total_dagather_t03 + t3
        if line.find("[CCPL <da_ensgather_demo>] total time") != -1:
            line3 = line.split(':')
            t4 = float(line3[-1])
            number_dagather_02 = number_dagather_02 + 1
            total_dagather_t04 = total_dagather_t04 + t4
        if line.find("ensgather total time") != -1:
        	line3 = line.split(':')
        	t3 = float(line3[-1])
        	number3 = number3 + 1
        	total_t3 = total_t3 + t3
total_da_t01 = total_da_t01/number_da_01
total_da_t02 = total_da_t02/number_da_01
total_da_t03 = total_da_t03/number_da_01
total_da_t04 = total_da_t04/number_da_02
total_damean_t01 = total_damean_t01/number_damean_01
total_damean_t02 = total_damean_t02/number_damean_01
total_damean_t03 = total_damean_t03/number_damean_01
total_damean_t04 = total_damean_t04/number_damean_02
total_dagather_t01 = total_dagather_t01/number_dagather_01
total_dagather_t02 = total_dagather_t02/number_dagather_01
total_dagather_t03 = total_dagather_t03/number_dagather_01
total_dagather_t04 = total_dagather_t04/number_dagather_02
total_t1 = total_t1/number1
total_t2 = total_t2/number2
total_t3 = total_t3/number3

result_da = []
#result_da.append(number_da_01)
#result_da.append('%.6f'%total_da_t01)
#result_da.append('%.6f'%total_da_t02)
#result_da.append('%.6f'%total_da_t03)
#result_da.append(number_da_02)
#result_da.append('%.6f'%total_da_t04)
#result_da.append(number1)
result_da.append('%.6f'%(total_t1-total_da_t04+total_da_t01))

result_damean = []
#result_damean.append(number_damean_01)
#result_damean.append('%.6f'%total_damean_t01)
#result_damean.append('%.6f'%total_damean_t02)
#result_damean.append('%.6f'%total_damean_t03)
#result_damean.append(number_damean_02)
#result_damean.append('%.6f'%total_damean_t04)
#result_damean.append(number2)
#result_damean.append('%.6f'%(total_t2-total_damean_t04+total_damean_t01))
result_da.append('%.6f'%(total_t2-total_damean_t04+total_damean_t01))

result_dagather = []
#result_dagather.append(number_dagather_01)
#result_dagather.append('%.6f'%total_dagather_t01)
#result_dagather.append('%.6f'%total_dagather_t02)
#result_dagather.append('%.6f'%total_dagather_t03)
#result_dagather.append(number_dagather_02)
#result_dagather.append('%.6f'%total_dagather_t04)
#result_dagather.append(number3)
#result_dagather.append('%.6f'%(total_t3-total_dagather_t04+total_dagather_t01))
result_da.append('%.6f'%(total_t3-total_dagather_t04+total_dagather_t01))
#print("Model individual",result_da)
#print("Model ensmean   ",result_damean)
#print("Model ensgather ",result_dagather)
print("CCPL_intial",result_da)
