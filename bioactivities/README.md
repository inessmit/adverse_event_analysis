A local installation of the target prediction tool PIDGIN was used <sup id="a1">[1](#f1)</sup><sup id="a2">[,2](#f2)</sup>.

Used PIDGIN environment as described in documentation of PIDGIN [1] for executing target prediction.

Rest of notebooks in 'release environment' as detailed in the top-level README, except 03_prepare_pidgin_input with rdkit environment (version 2019.09.2.0), which was created with:
```
conda create -c rdkit -n my-rdkit-env rdkit python=3.6 xlrd nb_conda_kernels pymysql
```

<b id="f1">1</b> https://pidginv3.readthedocs.io/en/latest/ [↩](#a1)  
<b id="f2">2</b> Mervin LH, Bulusu KC, Kalash L, et al. Orthologue chemical space and its influence on target prediction. Bioinformatics. 2018;34(1):72‐79. doi: [10.1093/bioinformatics/btx525](https://doi.org/10.1093/bioinformatics/btx525) [↩](#a2)