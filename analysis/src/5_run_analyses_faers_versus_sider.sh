jupyter nbconvert --to script run_analyses_faers_versus_sider.ipynb
jupyter nbconvert --to script analysis_functions.ipynb
chmod u+x run_analyses_faers_versus_sider.py
chmod u+x analysis_functions.py
./run_analyses_faers_versus_sider.py --faers_dir 20200110_faers_unbound_margin_pred_005_PRR2 --sider_dir 20200110_sider_unbound_margin_pred --dest_dir unbound_margin_pred_faers_vs_sider
#./run_analyses_faers_versus_sider.py --faers_dir 20200110_faers_cutoff6_pred_005_PRR2 --sider_dir 20200110_sider_cutoff6_pred --dest_dir cutoff_pred_faers_vs_sider
