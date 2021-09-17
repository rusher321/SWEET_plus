# SWEEET
1.SWEET_sample_weight_calculating.py

$ python3 1.SWEET_sample_weight_calculating.py -f ./example/expression.txt -s ./example/weight.txt

Argument | Variable | Description | Default value
------------ | ------------- | ------------- | -------------
-h | help | Get help for any of the commands | NA
-f | file | expression file | ./example/expression.txt
-k | k | Balance parameter | 0.1
-s | save | name of save file | ./example/weight.txt

2.SWEET_edge_weight_calculating.py

$ python3 2.SWEET_edge_weight_calculating.py -f ./example/expression.txt -w ./example/weight.txt -p ./example/patient.txt -g ./example/gene.txt -s ./example

Argument | Variable | Description | Default value
------------ | ------------- | ------------- | -------------
-h | help | Get help for any of the commands | NA
-f | file_e | expression file | ./example/expression.txt
-w | file_w | weight file | ./example/weight.txt
-p | file_p | patient file | ./example/patient.txt
-g | file_g | gene file | ./example/gene.txt
-s | save_path | save path | ./example

3.SWEET_calculating_mean_std.py

$ python3 3.SWEET_calculating_mean_std.py -p ./example/patient.txt -l  ./example -s ./example/mean_std.txt

Argument | Variable | Description | Default value
------------ | ------------- | ------------- | -------------
-h | help | Get help for any of the commands | NA
-p | file_p | patient file | ./example/patient.txt
-l | file_l | Patient edge weight file path | ./example
-s | save | name of save file | ./example/mean_std.txt
