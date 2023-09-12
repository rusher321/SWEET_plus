python3 1.SWEET_sample_weight_calculating.py -f ./example/mf_bacterime.txt -f2 ./example/mf_vrome.txt -s ./test2/weight.spearman -c spearman
python3 2.SWEET_edge_score_calculating.py -f ./example/mf_bacterime.txt -f2 ./example/mf_vrome.txt -w ./test2/weight.spearman.mean.txt -p ./example/patient_mf.txt -g ./example/species.txt -c spearman -s ./test2/
python3 3.SWEET_calculating_mean_std_zscore.py -p ./example/patient_mf.txt -l  ./test2 -s ./test2/mean_std.txt -z True
